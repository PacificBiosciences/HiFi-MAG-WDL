version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow assign_summarize_taxonomy {
	input {
		String sample_id

		File gtdb_batch_txt
		File gtdbtk_data_tar_gz
		Array[File] dereplicated_bins

		File filtered_quality_report_tsv

		Int min_mag_completeness
		Int max_mag_contamination


		RuntimeAttributes default_runtime_attributes
	}

	call assign_taxonomy {
		input:
			sample_id = sample_id,
			gtdb_batch_txt = gtdb_batch_txt,
			gtdbtk_data_tar_gz = gtdbtk_data_tar_gz,
			dereplicated_bins= dereplicated_bins,
			runtime_attributes = default_runtime_attributes
	}

	call mag_summary {
		input:
			sample_id = sample_id,
			gtdbtk_summary_txt = assign_taxonomy.gtdbtk_summary_txt,
			filtered_quality_report_tsv = filtered_quality_report_tsv,
			runtime_attributes = default_runtime_attributes
	}

	call mag_copy {
		input:
			mag_summary_txt = mag_summary.mag_summary_txt,
			dereplicated_bins = dereplicated_bins,
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
		File gtdbtk_summary_txt = assign_taxonomy.gtdbtk_summary_txt
		File gtdbk_output_tar_gz = assign_taxonomy.gtdbk_output_tar_gz
		File mag_summary_txt = mag_summary.mag_summary_txt
		Array[File] filtered_mags_fastas = mag_copy.filtered_mags_fastas
		File dastool_bins_plot_pdf = mag_plots.dastool_bins_plot_pdf
		File contigs_quality_plot_pdf = mag_plots.contigs_quality_plot_pdf
		File genome_size_depths_plot_df = mag_plots.genome_size_depths_plot_df
	}
}

task assign_taxonomy {
	input {
		String sample_id

		File gtdb_batch_txt
		File gtdbtk_data_tar_gz
		Array[File] dereplicated_bins

		RuntimeAttributes runtime_attributes
	}

	Int threads = 48
	Int mem_gb = threads * 2
	Int disk_size = ceil((size(gtdbtk_data_tar_gz, "GB") + (size(dereplicated_bins[0], "GB") * length(dereplicated_bins))) * 2 + 20)

	command <<<
		set -euo pipefail

		# Must set $GTDBTK_DATA_PATH variable to use gtdbtk command
		tar -xzvf ~{gtdbtk_data_tar_gz}
		GTDBTK_DATA_PATH="$(pwd)/$(tar -tzf ~{gtdbtk_data_tar_gz} | head -1 | cut -d '/' -f 1)"
		export GTDBTK_DATA_PATH

		gtdbtk --version

		# Ensure all bins are in the bin_dir; this will match the structure of the gtdb_batch_txt
		mkdir bin_dir
		while read -r bin || [[ -n "${bin}" ]]; do
			ln -s "${bin}" "$(pwd)/bin_dir"
		done < ~{write_lines(dereplicated_bins)}

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
		docker: "~{runtime_attributes.container_registry}/gtdbtk:2.1.1"
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
		Array[File] dereplicated_bins

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(mag_summary_txt, "GB") + (size(dereplicated_bins[0], "GB") * length(dereplicated_bins))) * 2 + 20)

	command <<<
		set -euo pipefail

		# Ensure all bins are in the bin_dir
		mkdir bin_dir
		while read -r bin || [[ -n "${bin}" ]]; do
			ln -s "${bin}" "$(pwd)/bin_dir"
		done < ~{write_lines(dereplicated_bins)}

		mkdir filtered_mags_out_dir

		python /opt/scripts/Copy-Final-MAGs.py \
			--mag_summary ~{mag_summary_txt} \
			--magdir bin_dir \
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