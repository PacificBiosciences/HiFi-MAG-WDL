version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "completeness_aware_binning/completeness_aware_binning.wdl" as CompletenessAwareBinning
import "coverage/coverage.wdl" as Coverage
import "binning/binning.wdl" as Binning
import "checkm2/checkm2.wdl" as CheckM2
import "gtdbtk/gtdbtk.wdl" as GTDBTk
import "mag/mag.wdl" as MAG

workflow metagenomics {
	input {
		String sample_id
		File? bam
		File? fastq

		# Complete-aware binning
		File checkm2_ref_db
		Int min_contig_length = 500000
		Int min_contig_completeness = 93

		# Binning
		Int metabat2_min_contig_size = 30000
		String semibin2_model_flag = "--environment=global"
		String dastool_search_engine = "diamond"
		Float dastool_score_threshold = 0.2

		Int min_mag_completeness = 70
		Int max_mag_contamination = 10
		Int max_contigs = 20

		# GTDB-Tk reference data
		File gtdbtk_data_tar_gz

		# Backend configuration
		String backend
		String? zones
		String? aws_spot_queue_arn
		String? aws_on_demand_queue_arn
		String? container_registry

		Boolean preemptible
	}

	call BackendConfiguration.backend_configuration {
		input:
			backend = backend,
			zones = zones,
			aws_spot_queue_arn = aws_spot_queue_arn,
			aws_on_demand_queue_arn = aws_on_demand_queue_arn,
			container_registry = container_registry
	}

	RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

	Boolean run_bam_to_fastq = if (defined(bam)) then true else false

	if (run_bam_to_fastq) {
		call bam_to_fastq {
			input:
				sample_id = sample_id,
				bam = select_first([bam]),
				runtime_attributes = default_runtime_attributes
		}
	}

	call hifiasm_meta {
		input:
			sample_id = sample_id,
			fastq = select_first([bam_to_fastq.converted_fastq, fastq]),
			runtime_attributes = default_runtime_attributes
	}

	call CompletenessAwareBinning.completeness_aware_binning {
		input:
			sample_id = sample_id,
			contigs_fasta = hifiasm_meta.primary_contig_fasta,
			min_contig_length = min_contig_length,
			min_contig_completeness = min_contig_completeness,
			checkm2_ref_db = checkm2_ref_db,
			default_runtime_attributes = default_runtime_attributes
	}

	call Coverage.coverage {
		input:
			sample_id = sample_id,
			contigs_fasta = hifiasm_meta.primary_contig_fasta,
			hifi_reads_fasta = hifiasm_meta.hifi_reads_fasta,
			bins_contigs_key_txt = completeness_aware_binning.bins_contigs_key_txt,
			default_runtime_attributes = default_runtime_attributes
	}

	call Binning.binning {
		input:
			sample_id = sample_id,
			incomplete_contigs_fasta = completeness_aware_binning.incomplete_contigs_fasta,
			filtered_contig_depth_txt = coverage.filtered_contig_depth_txt,
			sorted_bam = coverage.sorted_bam.data,
			metabat2_min_contig_size = metabat2_min_contig_size,
			semibin2_model_flag = semibin2_model_flag,
			dastool_search_engine = dastool_search_engine,
			dastool_score_threshold = dastool_score_threshold,
			default_runtime_attributes = default_runtime_attributes
	}

	call CheckM2.checkm2 {
		input:
			sample_id = sample_id,
			checkm2_ref_db = checkm2_ref_db,
			filtered_contig_depth_txt = coverage.filtered_contig_depth_txt,
			long_bin_fastas = completeness_aware_binning.long_bin_fastas,
			dastool_bins = binning.dastool_bins,
			min_mag_completeness = min_mag_completeness,
			max_mag_contamination = max_mag_contamination,
			max_contigs = max_contigs,
			default_runtime_attributes = default_runtime_attributes
	}

	if (checkm2.passed_bin_count_nonempty) {
		call GTDBTk.gtdbtk {
			input:
				sample_id = sample_id,
				gtdb_batch_txt = checkm2.gtdb_batch_txt,
				gtdbtk_data_tar_gz = gtdbtk_data_tar_gz,
				derep_bins = checkm2.derep_bins,
				default_runtime_attributes = default_runtime_attributes
		}

		call MAG.mag {
			input:
				sample_id = sample_id,
				gtdbtk_summary_txt = gtdbtk.gtdbtk_summary_txt,
				filtered_quality_report_tsv = checkm2.filtered_quality_report_tsv,
				derep_bins = checkm2.derep_bins,
				min_mag_completeness = min_mag_completeness,
				max_mag_contamination = max_mag_contamination,
				default_runtime_attributes = default_runtime_attributes
		}
	}

	output {
		# Preprocessing
		File? converted_fastq = bam_to_fastq.converted_fastq
		File primary_contig_gfa = hifiasm_meta.primary_contig_gfa
		File primary_contig_fasta = hifiasm_meta.primary_contig_fasta
		File hifi_reads_fasta = hifiasm_meta.hifi_reads_fasta

		# Completeness-aware binning output
		File bins_contigs_key_txt = completeness_aware_binning.bins_contigs_key_txt
		Array[File] long_bin_fastas = completeness_aware_binning.long_bin_fastas
		File incomplete_contigs_fasta = completeness_aware_binning.incomplete_contigs_fasta
		File? contig_quality_report_tsv = completeness_aware_binning.contig_quality_report_tsv
		File? passed_bins_txt = completeness_aware_binning.passed_bins_txt
		File? scatterplot_pdf = completeness_aware_binning.scatterplot_pdf
		File? histogram_pdf = completeness_aware_binning.histogram_pdf

		# Coverage output
		IndexData sorted_bam = coverage.sorted_bam
		File filtered_contig_depth_txt = coverage.filtered_contig_depth_txt

		# Binning output
		Array[File] metabat2_reconstructed_bins_fastas = binning.metabat2_reconstructed_bins_fastas
		File metabat2_bin_sets_tsv = binning.metabat2_bin_sets_tsv
		File semibin2_bins_tsv = binning.semibin2_bins_tsv
		Array[File] semibin2_reconstructed_bins_fastas = binning.semibin2_reconstructed_bins_fastas
		File semibin2_bin_sets_tsv = binning.semibin2_bin_sets_tsv
		Array[File] dastool_bins = binning.dastool_bins

		# CheckM2 output
		Array[File] derep_bins = checkm2.derep_bins
		File bin_quality_report_tsv = checkm2.bin_quality_report_tsv
		File gtdb_batch_txt = checkm2.gtdb_batch_txt
		File passed_bin_count_txt = checkm2.passed_bin_count_txt
		File filtered_quality_report_tsv = checkm2.filtered_quality_report_tsv

		# GTDB-Tk output
		File? gtdbtk_align_tar_gz = gtdbtk.gtdbtk_align_tar_gz
		File? gtdbtk_classify_tar_gz = gtdbtk.gtdbtk_classify_tar_gz
		File? gtdbtk_identify_tar_gz = gtdbtk.gtdbtk_identify_tar_gz
		File? gtdbtk_summary_txt = gtdbtk.gtdbtk_summary_txt

		# MAG summary and plot output
		File? mag_summary_txt = mag.mag_summary_txt
		Array[File]? filtered_mags_fastas = mag.filtered_mags_fastas
		File? dastool_bins_plot_pdf = mag.dastool_bins_plot_pdf
		File? contigs_quality_plot_pdf = mag.contigs_quality_plot_pdf
		File? genome_size_depths_plot_df = mag.genome_size_depths_plot_df
	}

	parameter_meta {
		sample_id: {help: "Sample ID"}
		bam: {help: "Optional sample BAM; one of [bam, fastq] must be provided as input"}
		fastq: {help: "Optional sample in FASTQ format; one of [bam, fastq] must be provided as input"}
		run_bam_to_fastq: {help: "Optional step to convert sample BAM to FASTQ format if BAM is provided"}
		# Completeness-aware binning
		checkm2_ref_db: {help: "CheckM2 DIAMOND reference database Uniref100/KO"}
		min_contig_length: {help: "Minimum size of a contig to consider for completeness scores; default value is set to 500kb. This value should not be increased"}
		min_contig_completeness: {help: "Minimum completeness score (from CheckM2) to mark a contig as complete and place it in a distinct bin; default value is set to 93%. This value should not be lower than 90%"}
		# Binning
		metabat2_min_contig_size: {help: "The minimum size of contig to be included in binning for MetaBAT2; default value is set to 30000"}
		semibin2_model_flag: {help: "The trained model to be used in SemiBin2; default value is set to 'global'"}
		dastool_search_engine: {help: "The engine for single copy gene searching used in DAS Tool; default is set to 'diamond'"}
		dastool_score_threshold: {help: "Score threshold until selection algorithm will keep selecting bins [0 to 1] used in DAS Tool; default value is set to 0.2 (20%)"}
		# Quality filters for MAGs
		min_mag_completeness: {help: "Minimum completeness score for a genome bin; default value is set to 70%"}
		max_mag_contamination: {help: "Maximum contamination threshold for a genome bin; default value is set to 10%"}
		max_contigs: {help: "The maximum number of contigs allowed in a genome bin; default value is set to 20"}
		# GTDBT-k
		gtdbtk_data_tar_gz: {help: "A TAR GZ file of GTDB-Tk (Genome Database Taxonomy toolkit) reference data, release207_v2 used for assigning taxonomic classifications to bacterial and archaeal genomes"}
		# Backend configuration
		backend: {help: "Backend where the workflow will be executed ['GCP', 'Azure', 'AWS']"}
		zones: {help: "Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'"}
		aws_spot_queue_arn: {help: "Queue ARN for the spot batch queue; required if backend is set to 'AWS'"}
		aws_on_demand_queue_arn: {help: "Queue ARN for the on demand batch queue; required if backend is set to 'AWS'"}
		container_registry: {help: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used"}
		preemptible: {help: "Where possible, run tasks preemptibly"}
	}
}

task bam_to_fastq {
	input {
		String sample_id
		File bam

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		samtools --version

		samtools fastq \
			~{bam} \
		| bgzip -c > "~{sample_id}.fastq.gz"
	>>>

	output {
		File converted_fastq = "~{sample_id}.fastq.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools:5e8307c"
		cpu: 2
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}

task hifiasm_meta {
	input {
		String sample_id
		File fastq

		RuntimeAttributes runtime_attributes
	}

	Int threads = 32
	Int mem_gb = threads * 4
	Int disk_size = ceil(size(fastq, "GB") * 4 + 20)

	command <<<
		set -euo pipefail

		hifiasm_meta --version

		hifiasm_meta \
			-t ~{threads} \
			-o ~{sample_id} \
			~{fastq}

		awk '/^S/{print ">"$2;print $3}' "~{sample_id}.p_ctg.gfa" > "~{sample_id}.p_ctg.fa"

		gzip -d <~{fastq} | sed -n '1~4s/^@/>/p;2~4p' > "~{sample_id}.fa"
	>>>

	output {
		File primary_contig_gfa = "~{sample_id}.p_ctg.gfa"
		File primary_contig_fasta = "~{sample_id}.p_ctg.fa"
		File hifi_reads_fasta = "~{sample_id}.fa"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/hifiasm-meta:0.3"
		cpu: threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
