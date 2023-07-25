version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "completeness_aware_binning/completeness_aware_binning.wdl" as CompletenessAwareBinning
import "coverage/coverage.wdl" as Coverage
import "binning/binning.wdl" as Binning


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
			reads_fasta = hifiasm_meta.reads_fasta,
			bins_contigs_key = completeness_aware_binning.bins_contigs_key,
			default_runtime_attributes = default_runtime_attributes
	}

	call Binning.binning {
		input:
			sample_id = sample_id,
			incomplete_contigs_fasta = completeness_aware_binning.incomplete_contigs,
			filtered_contig_depth_txt = coverage.filtered_contig_depth_txt,
			sorted_bam = coverage.sorted_bam,
			metabat2_min_contig_size = metabat2_min_contig_size,
			semibin2_model_flag = semibin2_model_flag,
			dastool_search_engine = dastool_search_engine,
			dastool_score_threshold = dastool_score_threshold,
			default_runtime_attributes = default_runtime_attributes
	}

	output {
		# BAM to FASTQ output
		File? converted_fastq = bam_to_fastq.converted_fastq

		# Assembly output
		File primary_contig_graph = hifiasm_meta.primary_contig_graph
		File primary_contig_fasta = hifiasm_meta.primary_contig_fasta
		File reads_fasta = hifiasm_meta.reads_fasta

		# Completeness-aware binning output
		File bins_contigs_key = completeness_aware_binning.bins_contigs_key
		File incomplete_contigs = completeness_aware_binning.incomplete_contigs
		File? report = completeness_aware_binning.report
		File? passed_bins = completeness_aware_binning.passed_bins
		File? scatterplot = completeness_aware_binning.scatterplot
		File? histogram = completeness_aware_binning.histogram

		# Coverage output
		File sorted_bam = coverage.sorted_bam
		File sorted_bam_index = coverage.sorted_bam_index
		File filtered_contig_depth_txt = coverage.filtered_contig_depth_txt

		# Binning output
		Array[File] metabat2_reconstructed_bins_fastas = binning.metabat2_reconstructed_bins_fastas
		File metabat2_bin_sets_tsv = binning.metabat2_bin_sets_tsv
		File semibin2_bins_tsv = binning.semibin2_bins_tsv
		Array[File] semibin2_reconstructed_bins_fastas = binning.semibin2_reconstructed_bins_fastas
		File semibin2_bin_sets_tsv = binning.semibin2_bin_sets_tsv
		Array[File] dastool_bins = binning.dastool_bins
	}

	parameter_meta {
		sample_id: {help: "Sample ID"}
		# BAM to FASTQ
		bam: {help: "Optional sample BAM to convert to FASTQ format; one of [bam, fastq] must be provided as input"}
		# Assembly
		fastq: {help: "Optional sample in FASTQ format; one of [bam, fastq] must be provided as input"}
		# Completeness-aware binning
		checkm2_ref_db: {help: "CheckM2 DIAMOND reference database Uniref100/KO"}
		min_contig_length: {help: "Minimum size of a contig to consider for completeness scores; default value is set to 500kb. This value should not be increased"}
		min_contig_completeness: {help: "Minimum completeness score (from CheckM2) to mark a contig as complete and place it in a distinct bin; default value is set to 93%. This value should not be lower than 90%"}
		# Binning
		metabat2_min_contig_size: {help: "The minimum size of contig to be included in binning for MetaBAT2; default value is set to 30000"}
		semibin2_model_flag: {help: "The trained model to be used in SemiBin2; default value is set to 'global'"}
		dastool_search_engine: {help: "The engine for single copy gene searching used in DAS Tool; default is set to 'diamond'"}
		dastool_score_threshold: {help: "Score threshold until selection algorithm will keep selecting bins [0 to 1] used in DAS Tool; default value is set to 0.2 (20%)"}
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
	Int mem_gb = threads * 2
	Int disk_size = ceil(size(fastq, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		hifiasm_meta --version

		hifiasm_meta \
			-t ~{threads} \
			-o ~{sample_id} \
			~{fastq}

		awk '/^S/{print ">"$2;print $3}' "~{sample_id}.p_ctg.gfa" > "~{sample_id}.p_ctg.fa"
	>>>

	output {
		File primary_contig_graph = "~{sample_id}.p_ctg.gfa"
		File primary_contig_fasta = "~{sample_id}.p_ctg.fa"
		File reads_fasta = "~{sample_id}.rescue.fa"
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
