version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow mag {
	input {
		String sample_id

		File gtdbk_summary_txt
		File filtered_quality_report_tsv
		Array[File] derep_bins

		Int min_mag_completeness
		Int max_mag_contamination

		RuntimeAttributes default_runtime_attributes
	}

	call mag_summary {
		input:
			sample_id = sample_id,
			gtdbk_summary_txt = gtdbk_summary_txt,
			filtered_quality_report_tsv = filtered_quality_report_tsv,
			runtime_attributes = default_runtime_attributes
	}

	call mag_copy {
		input:
			mag_summary_txt = mag_summary.mag_summary_txt,
			derep_bins = derep_bins,
			runtime_attributes = default_runtime_attributes
	}

	call mag_plots {
		input:
			sample_id = sample_id,
			filtered_quality_report_tsv = filtered_quality_report_tsv,
			mag_summary_txt = mag_summary.mag_summary_txt,
			min_mag_completeness = min_mag_completeness,
			max_mag_contamination = max_mag_contamination,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File mag_summary_txt = mag_summary.mag_summary_txt
		Array[File] filtered_mags_fastas = mag_copy.filtered_mags_fastas
		File dastool_bins_plot_pdf = mag_plots.dastool_bins_plot_pdf
		File contigs_quality_plot_pdf = mag_plots.contigs_quality_plot_pdf
		File genome_size_depths_plot_df = mag_plots.genome_size_depths_plot_df
	}
}

task mag_summary {
	input {
		String sample_id

		File gtdbk_summary_txt
		File filtered_quality_report_tsv

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(gtdbk_summary_txt, "GB") + size(filtered_quality_report_tsv, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/MAG-Summary.py \
			--gtdb_summary ~{gtdbk_summary_txt} \
			--checmk2_summary ~{filtered_quality_report_tsv} \
			--outfile "~{sample_id}.HiFi_MAG.summary.txt"
	>>>

	output {
		File mag_summary_txt = "~{sample_id}.HiFi_MAG.summary.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python:5e8307c"
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

task mag_copy {
	input {
		File mag_summary_txt
		Array[File] derep_bins

		RuntimeAttributes runtime_attributes
	}
 
	Int disk_size = ceil((size(mag_summary_txt, "GB") + (size(derep_bins[0], "GB") * length(derep_bins))) * 2 + 20)

	command <<<
		set -euo pipefail

		derep_bins_dir=$(dirname ~{derep_bins[0]})

		mkdir filtered_mags_out_dir

		python /opt/scripts/Copy-Final-MAGs.py \
			--mag_summary ~{mag_summary_txt} \
			--magdir "${derep_bins_dir}" \
			--outdir filtered_mags_out_dir
	>>>

	output {
		Array[File] filtered_mags_fastas = glob("filtered_mags_out_dir/*.fa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python:5e8307c"
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

task mag_plots {
	input {
		String sample_id

		File filtered_quality_report_tsv
		File mag_summary_txt

		Int min_mag_completeness
		Int max_mag_contamination

		RuntimeAttributes runtime_attributes
	}
 
	Int disk_size = ceil((size(filtered_quality_report_tsv, "GB") + size(mag_summary_txt, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Plot-Figures.py \
			--checkm_eval ~{filtered_quality_report_tsv} \
			--mag_eval ~{mag_summary_txt} \
			--label ~{sample_id} \
			--completeness ~{min_mag_completeness} \
			--contamination ~{max_mag_contamination} \
			--output1 "~{sample_id}.All-DASTool-Bins.pdf" \
			--output2 "~{sample_id}.Completeness-Contamination-Contigs.pdf" \
			--output3 "~{sample_id}.GenomeSizes-Depths.pdf"
	>>>

	output {
		File dastool_bins_plot_pdf = "~{sample_id}.All-DASTool-Bins.pdf"
		File contigs_quality_plot_pdf = "~{sample_id}.Completeness-Contamination-Contigs.pdf"
		File genome_size_depths_plot_df = "~{sample_id}.GenomeSizes-Depths.pdf"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python:5e8307c"
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
