version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow assign_taxonomy {
	input {
		String sample_id

		File gtdb_batch_txt
		File gtdbtk_data_tar_gz
		Array[File] dereplicated_bin_fas_tars

		File filtered_quality_report_tsv

		Int min_mag_completeness
		Int max_mag_contamination


		RuntimeAttributes default_runtime_attributes
	}

	call assign_taxonomy_gtdbtk {
		input:
			sample_id = sample_id,
			gtdb_batch_txt = gtdb_batch_txt,
			gtdbtk_data_tar_gz = gtdbtk_data_tar_gz,
			dereplicated_bin_fas_tars = dereplicated_bin_fas_tars,
			runtime_attributes = default_runtime_attributes
	}

	call mag_summary {
		input:
			sample_id = sample_id,
			gtdbtk_summary_txt = assign_taxonomy_gtdbtk.gtdbtk_summary_txt,
			filtered_quality_report_tsv = filtered_quality_report_tsv,
			runtime_attributes = default_runtime_attributes
	}

	call mag_copy {
		input:
			mag_summary_txt = mag_summary.mag_summary_txt,
			dereplicated_bin_fas_tars = dereplicated_bin_fas_tars,
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
		File gtdbtk_summary_txt = assign_taxonomy_gtdbtk.gtdbtk_summary_txt
		File gtdbk_output_tar_gz = assign_taxonomy_gtdbtk.gtdbk_output_tar_gz
		File mag_summary_txt = mag_summary.mag_summary_txt
		Array[File] filtered_mags_fas = mag_copy.filtered_mags_fas
		File dastool_bins_plot_pdf = mag_plots.dastool_bins_plot_pdf
		File contigs_quality_plot_pdf = mag_plots.contigs_quality_plot_pdf
		File genome_size_depths_plot_df = mag_plots.genome_size_depths_plot_df
	}
}

task assign_taxonomy_gtdbtk {
	input {
		String sample_id

		File gtdb_batch_txt
		File gtdbtk_data_tar_gz
		Array[File] dereplicated_bin_fas_tars

		RuntimeAttributes runtime_attributes
	}

	Int threads = 48
	Int mem_gb = threads * 2
	Int disk_size = ceil((size(gtdbtk_data_tar_gz, "GB") + (size(dereplicated_bin_fas_tars, "GB") )) * 3 + 20)

	command <<<
		set -euo pipefail

		# Must set $GTDBTK_DATA_PATH variable to use gtdbtk command
		mkdir gtdbtk_data
		tar \
			-xzvf ~{gtdbtk_data_tar_gz} \
			-C gtdbtk_data \
			--strip-components 1
		GTDBTK_DATA_PATH="$(pwd)/gtdbtk_data"
		export GTDBTK_DATA_PATH

		gtdbtk --version

		mkdir bins
		while read -r bin_fas_tar || [[ -n "${bin_fas_tar}" ]]; do
			tar -zxvf "${bin_fas_tar}" \
				-C bins \
				--strip-components 1
		done < ~{write_lines(dereplicated_bin_fas_tars)}

		mkdir ~{sample_id}_gtdbtk tmp_dir

		gtdbtk classify_wf \
			--batchfile ~{gtdb_batch_txt} \
			--out_dir ~{sample_id}_gtdbtk \
			--extension fa \
			--prefix ~{sample_id} \
			--cpus ~{threads} \
			--tmpdir tmp_dir

		python /opt/scripts/GTDBTk-Organize.py \
			--input_dir ~{sample_id}_gtdbtk/classify \
			--outfile "~{sample_id}.GTDBTk_Summary.txt"

		tar -zcvf "~{sample_id}_gtdbtk.tar.gz" ~{sample_id}_gtdbtk
	>>>

	output {
		File gtdbtk_summary_txt = "~{sample_id}.GTDBTk_Summary.txt"
		File gtdbk_output_tar_gz = "~{sample_id}_gtdbtk.tar.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/gtdbtk@sha256:c2b2129af53057fdfc03a0837750e4c5b0562708dde3091747e099f9d30640b7"
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

task mag_summary {
	input {
		String sample_id

		File gtdbtk_summary_txt
		File filtered_quality_report_tsv

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(gtdbtk_summary_txt, "GB") + size(filtered_quality_report_tsv, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/MAG-Summary.py \
			--gtdb_summary ~{gtdbtk_summary_txt} \
			--checmk2_summary ~{filtered_quality_report_tsv} \
			--outfile "~{sample_id}.HiFi_MAG.summary.txt"
	>>>

	output {
		File mag_summary_txt = "~{sample_id}.HiFi_MAG.summary.txt"
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

task mag_copy {
	input {
		File mag_summary_txt
		Array[File] dereplicated_bin_fas_tars

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(mag_summary_txt, "GB") + (size(dereplicated_bin_fas_tars, "GB"))) * 3 + 20)

	command <<<
		set -euo pipefail

		mkdir bins
		while read -r bin_fas_tar || [[ -n "${bin_fas_tar}" ]]; do
			tar -zxvf "${bin_fas_tar}" \
				-C bins \
				--strip-components 1
		done < ~{write_lines(dereplicated_bin_fas_tars)}

		mkdir filtered_mags_out_dir

		python /opt/scripts/Copy-Final-MAGs.py \
			--mag_summary ~{mag_summary_txt} \
			--magdir "$(pwd)/bins" \
			--outdir filtered_mags_out_dir
	>>>

	output {
		Array[File] filtered_mags_fas = glob("filtered_mags_out_dir/*.fa")
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
