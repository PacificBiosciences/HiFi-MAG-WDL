version 1.0

import "../wdl-common/wdl/structs.wdl"
import "bin_long_contigs/bin_long_contigs.wdl" as BinLongContigs
import "bin_incomplete_contigs/bin_incomplete_contigs.wdl" as BinIncompleteContigs
import "common/predict_bin_quality.wdl" as PredictBinQuality

workflow bin_contigs {
	input {
		String sample_id
		File hifi_reads_fastq
		File assembled_contigs_fa
		File assembled_contigs_fa_gz

		File checkm2_ref_db
		Int min_contig_length = 500000
		Int min_contig_completeness = 93

		Int metabat2_min_contig_size = 30000
		String semibin2_model = "global"
		String dastool_search_engine = "diamond"
		Float dastool_score_threshold = 0.2

		Int min_mag_completeness = 70
		Int max_mag_contamination = 10
		Int max_contigs = 20

		RuntimeAttributes default_runtime_attributes
	}

	call BinLongContigs.bin_long_contigs {
		input:
			sample_id = sample_id,
			assembled_contigs_fa = assembled_contigs_fa,
			checkm2_ref_db = checkm2_ref_db,
			min_contig_length = min_contig_length,
			min_contig_completeness = min_contig_completeness,
			default_runtime_attributes = default_runtime_attributes
	}

	call BinIncompleteContigs.bin_incomplete_contigs {
		input:
			sample_id = sample_id,
			assembled_contigs_fa_gz = assembled_contigs_fa_gz,
			hifi_reads_fastq = hifi_reads_fastq,
			incomplete_contigs_fa = bin_long_contigs.incomplete_contigs_fa,
			passing_long_contig_bin_map = bin_long_contigs.passing_long_contig_bin_map,
			metabat2_min_contig_size = metabat2_min_contig_size,
			semibin2_model = semibin2_model,
			dastool_search_engine = dastool_search_engine,
			dastool_score_threshold = dastool_score_threshold,
			default_runtime_attributes = default_runtime_attributes
	}

	Array[File] derep_bin_fas_tars = [
		bin_long_contigs.filtered_long_bin_fas_tar,
		bin_incomplete_contigs.merged_incomplete_bin_fas_tar
	]

	call PredictBinQuality.predict_bin_quality {
		input:
			prefix = "~{sample_id}.dereplicated_bins",
			checkm2_ref_db = checkm2_ref_db,
			bin_fas_tars = derep_bin_fas_tars,
			runtime_attributes = default_runtime_attributes
	}

	call filter_dereplicated_bins {
		input:
			sample_id = sample_id,
			bin_quality_report_tsv = predict_bin_quality.bin_quality_report_tsv,
			contig_depth_txt = bin_incomplete_contigs.contig_depth_txt,
			dereplicated_bin_fas_tars = derep_bin_fas_tars,
			min_mag_completeness = min_mag_completeness,
			max_mag_contamination = max_mag_contamination,
			max_contigs = max_contigs,
			runtime_attributes = default_runtime_attributes
	}

	output {
		# bin_long_contigs output
		File long_contig_bin_map = bin_long_contigs.long_contig_bin_map
		File filtered_long_bin_fas_tar = bin_long_contigs.filtered_long_bin_fas_tar
		File incomplete_contigs_fa = bin_long_contigs.incomplete_contigs_fa

		File? long_contig_bin_quality_report_tsv = bin_long_contigs.long_contig_bin_quality_report_tsv
		File? filtered_long_contig_bin_map = bin_long_contigs.filtered_long_contig_bin_map
		File? scatterplot_pdf = bin_long_contigs.scatterplot_pdf
		File? histogram_pdf = bin_long_contigs.histogram_pdf

		File passing_long_contig_bin_map = bin_long_contigs.passing_long_contig_bin_map

		# bin_incomplete_contigs output
		# Coverage
		IndexData aligned_sorted_bam = bin_incomplete_contigs.aligned_sorted_bam
		File contig_depth_txt = bin_incomplete_contigs.contig_depth_txt

		# Incomplete contig binning
		File metabat2_bin_fas_tar = bin_incomplete_contigs.metabat2_bin_fas_tar
		File metabat2_contig_bin_map = bin_incomplete_contigs.metabat2_contig_bin_map

		File semibin2_bins_tsv = bin_incomplete_contigs.semibin2_bins_tsv
		File semibin2_bin_fas_tar = bin_incomplete_contigs.semibin2_bin_fas_tar
		File semibin2_contig_bin_map = bin_incomplete_contigs.semibin2_contig_bin_map

		File merged_incomplete_bin_fas_tar = bin_incomplete_contigs.merged_incomplete_bin_fas_tar

		Array[File] dereplicated_bin_fas_tars = derep_bin_fas_tars

		# Bin quality
		File bin_quality_report_tsv = predict_bin_quality.bin_quality_report_tsv
		File gtdb_batch_txt = filter_dereplicated_bins.gtdb_batch_txt
		File passed_bin_count_txt = filter_dereplicated_bins.passed_bin_count_txt
		File filtered_quality_report_tsv = filter_dereplicated_bins.filtered_quality_report_tsv
	}
}

task filter_dereplicated_bins {
	input {
		String sample_id

		File bin_quality_report_tsv
		File contig_depth_txt
		Array[File] dereplicated_bin_fas_tars

		Int min_mag_completeness
		Int max_mag_contamination
		Int max_contigs

		RuntimeAttributes runtime_attributes
	}

	String bin_quality_report_tsv_basename = basename(bin_quality_report_tsv, ".tsv")
	Int disk_size = ceil(size(dereplicated_bin_fas_tars, "GB") * 3 + 20)

	command <<<
		set -euo pipefail

		mkdir bins
		while read -r bin_fas_tar || [[ -n "${bin_fas_tar}" ]]; do
			tar -zxvf "${bin_fas_tar}" \
				-C bins \
				--strip-components 1
		done < ~{write_lines(dereplicated_bin_fas_tars)}

		python /opt/scripts/Filter-Checkm2-Bins.py \
			--input_tsv ~{bin_quality_report_tsv} \
			--bin_dir "$(pwd)/bins" \
			--depth_file ~{contig_depth_txt} \
			--min_completeness ~{min_mag_completeness} \
			--max_contamination ~{max_mag_contamination} \
			--max_contigs ~{max_contigs} \
			--gtdb_outfile "~{sample_id}.GTDBTk_batch_file.txt" \
			--target_outfile "~{sample_id}.BinCount.txt" \
			--updated_tsv "~{bin_quality_report_tsv_basename}.filtered.tsv"
	>>>

	output {
		File gtdb_batch_txt = "~{sample_id}.GTDBTk_batch_file.txt"
		File passed_bin_count_txt = "~{sample_id}.BinCount.txt"
		File filtered_quality_report_tsv = "~{bin_quality_report_tsv_basename}.filtered.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python@sha256:23d8a139685afad5b0c763ded6db6f02dd3884ad150a142bc1cd17a5d21910d5"
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
