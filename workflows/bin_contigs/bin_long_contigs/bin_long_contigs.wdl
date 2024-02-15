version 1.0

import "../../wdl-common/wdl/structs.wdl"
import "../common/predict_bin_quality.wdl" as PredictBinQuality

workflow bin_long_contigs {
	input {
		String sample_id
		File assembled_contigs_fa

		File checkm2_ref_db
		Int min_contig_length
		Int min_contig_completeness

		RuntimeAttributes default_runtime_attributes
	}

	call long_contigs_to_bins {
		input:
			sample_id = sample_id,
			assembled_contigs_fa = assembled_contigs_fa,
			min_contig_length = min_contig_length,
			runtime_attributes = default_runtime_attributes
	}

	if (read_int(long_contigs_to_bins.long_bin_fas_count) > 0) {
		call PredictBinQuality.predict_bin_quality {
			input:
				prefix = "~{sample_id}.long_contigs",
				checkm2_ref_db = checkm2_ref_db,
				bin_fas_tars = [long_contigs_to_bins.long_bin_fas_tar],
				runtime_attributes = default_runtime_attributes
		}

		# Filter long contigs by completeness percent
		call filter_complete_contigs {
			input:
				sample_id = sample_id,
				assembled_contigs_fa = assembled_contigs_fa,
				bin_quality_report_tsv = predict_bin_quality.bin_quality_report_tsv,
				contig_bin_map = long_contigs_to_bins.long_contig_bin_map,
				min_contig_completeness = min_contig_completeness,
				runtime_attributes = default_runtime_attributes
			}
	}

	File passing_contig_bin_map = select_first([filter_complete_contigs.filtered_long_contig_bin_map, long_contigs_to_bins.long_contig_bin_map])

	call make_incomplete_contigs {
		input:
			sample_id = sample_id,
			assembled_contigs_fa = assembled_contigs_fa,
			passing_contig_bin_map = passing_contig_bin_map,
			long_bin_fas_tar = long_contigs_to_bins.long_bin_fas_tar,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File long_contig_bin_map = long_contigs_to_bins.long_contig_bin_map

		File? long_contig_bin_quality_report_tsv = predict_bin_quality.bin_quality_report_tsv
		File? filtered_long_contig_bin_map = filter_complete_contigs.filtered_long_contig_bin_map
		File? scatterplot_pdf = filter_complete_contigs.scatterplot_pdf
		File? histogram_pdf = filter_complete_contigs.histogram_pdf

		File passing_long_contig_bin_map = passing_contig_bin_map

		File filtered_long_bin_fas_tar = make_incomplete_contigs.filtered_long_bin_fas_tar
		File incomplete_contigs_fa = make_incomplete_contigs.incomplete_contigs_fa
	}
}

task long_contigs_to_bins {
	input {
		String sample_id
		File assembled_contigs_fa

		Int min_contig_length

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(assembled_contigs_fa, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Fasta-Make-Long-Seq-Bins.py \
			--input_fasta ~{assembled_contigs_fa} \
			--bins_contigs "~{sample_id}.long_contig_bin_map.tsv" \
			--length ~{min_contig_length} \
			--outdir long_bin_fas \
			--prefix ~{sample_id}

		tar -zcvf ~{sample_id}.long_bin_fas.tar.gz long_bin_fas/

		find long_bin_fas/ -type f -name "*.fa" | wc -l > ~{sample_id}.long_bin_fas_count.txt
	>>>

	output {
		File long_contig_bin_map = "~{sample_id}.long_contig_bin_map.tsv"
		File long_bin_fas_tar = "~{sample_id}.long_bin_fas.tar.gz"
		File long_bin_fas_count = "~{sample_id}.long_bin_fas_count.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python@sha256:e76da216b7adcc498c35e64220a838a3a17d4e51420e59cfc18aa2064e2ef1f7"
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

task filter_complete_contigs {
	input {
		String sample_id
		File assembled_contigs_fa

		File bin_quality_report_tsv
		File contig_bin_map

		Int min_contig_completeness

		RuntimeAttributes runtime_attributes
	}

	String contig_bin_map_basename = basename(contig_bin_map, ".tsv")
	Int disk_size = ceil(size(assembled_contigs_fa, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Filter-Complete-Contigs.py \
			--input_fasta ~{assembled_contigs_fa} \
			--checkm ~{bin_quality_report_tsv} \
			--bins_contigs ~{contig_bin_map} \
			--min_completeness ~{min_contig_completeness} \
			--passed_bins "~{contig_bin_map_basename}.filtered.tsv" \
			--plot_scatter "~{sample_id}.completeness_vs_size_scatter.pdf" \
			--plot_histo "~{sample_id}.completeness_histo.pdf"
	>>>

	output {
		File filtered_long_contig_bin_map = "~{contig_bin_map_basename}.filtered.tsv"
		File scatterplot_pdf = "~{sample_id}.completeness_vs_size_scatter.pdf"
		File histogram_pdf = "~{sample_id}.completeness_histo.pdf"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python@sha256:e76da216b7adcc498c35e64220a838a3a17d4e51420e59cfc18aa2064e2ef1f7"
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

task make_incomplete_contigs {
	input {
		String sample_id
		File assembled_contigs_fa

		File passing_contig_bin_map
		File long_bin_fas_tar

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(assembled_contigs_fa, "GB") + size(long_bin_fas_tar, "GB") * 3 + 20)

	command <<<
		set -euo pipefail

		mkdir bins
		tar -zxvf ~{long_bin_fas_tar} \
			-C bins \
			--strip-components 1

		python /opt/scripts/Make-Incomplete-Contigs.py \
			--input_fasta ~{assembled_contigs_fa} \
			--output_fasta "~{sample_id}.incomplete_contigs.fa" \
			--passed_bins ~{passing_contig_bin_map} \
			--fastadir "$(pwd)/bins" \
			--outdir filtered_long_bin_fas

		tar -zcvf ~{sample_id}.filtered_long_bin_fas.tar.gz filtered_long_bin_fas/
	>>>

	output {
		File filtered_long_bin_fas_tar = "~{sample_id}.filtered_long_bin_fas.tar.gz"
		File incomplete_contigs_fa = "~{sample_id}.incomplete_contigs.fa"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python@sha256:e76da216b7adcc498c35e64220a838a3a17d4e51420e59cfc18aa2064e2ef1f7"
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
