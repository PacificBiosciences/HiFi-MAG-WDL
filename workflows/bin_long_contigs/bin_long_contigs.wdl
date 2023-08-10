version 1.0

import "../wdl-common/wdl/structs.wdl"
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

	if (length(long_contigs_to_bins.long_bin_fas) > 0) {
		call PredictBinQuality.predict_bin_quality {
			input:
				prefix = "~{sample_id}.long_contigs",
				checkm2_ref_db = checkm2_ref_db,
				bin_fas = long_contigs_to_bins.long_bin_fas,
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
			long_bin_fas = long_contigs_to_bins.long_bin_fas,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File long_contig_bin_map = long_contigs_to_bins.long_contig_bin_map

		File? long_contig_bin_quality_report_tsv = predict_bin_quality.bin_quality_report_tsv
		File? filtered_long_contig_bin_map = filter_complete_contigs.filtered_long_contig_bin_map
		File? scatterplot_pdf = filter_complete_contigs.scatterplot_pdf
		File? histogram_pdf = filter_complete_contigs.histogram_pdf

		File passing_long_contig_bin_map = passing_contig_bin_map

		Array[File] filtered_long_bin_fas = make_incomplete_contigs.filtered_long_bin_fas
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
	>>>

	output {
		File long_contig_bin_map = "~{sample_id}.long_contig_bin_map.tsv"
		Array[File] long_bin_fas = glob("long_bin_fas/~{sample_id}.*.fa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python@sha256:c7e594d86c35d2c3b2cd8fabf51d9274d74347464433c4f3e55e5306be7bd1ea"
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
		docker: "~{runtime_attributes.container_registry}/python@sha256:c7e594d86c35d2c3b2cd8fabf51d9274d74347464433c4f3e55e5306be7bd1ea"
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
		Array[File] long_bin_fas

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(assembled_contigs_fa, "GB") + size(long_bin_fas[0], "GB") * length(long_bin_fas) * 2 + 20)

	command <<<
		set -euo pipefail

		long_bin_fa_dir=$(dirname ~{long_bin_fas[0]})

		python /opt/scripts/Make-Incomplete-Contigs.py \
			--input_fasta ~{assembled_contigs_fa} \
			--output_fasta "~{sample_id}.incomplete_contigs.fa" \
			--passed_bins ~{passing_contig_bin_map} \
			--fastadir "${long_bin_fa_dir}" \
			--outdir filtered_long_bin_fa_dir
	>>>

	output {
		Array[File] filtered_long_bin_fas = glob("filtered_long_bin_fa_dir/*.fa")
		File incomplete_contigs_fa = "~{sample_id}.incomplete_contigs.fa"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python@sha256:c7e594d86c35d2c3b2cd8fabf51d9274d74347464433c4f3e55e5306be7bd1ea"
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
