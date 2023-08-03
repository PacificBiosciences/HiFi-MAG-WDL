version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow gtdbtk {
	input {
		String sample_id

		File gtdb_batch_txt
		File gtdbtk_data_tar_gz
		Array[File] derep_bins

		RuntimeAttributes default_runtime_attributes
	}

	call gtdbtk_analysis {
		input:
			sample_id = sample_id,
			gtdb_batch_txt = gtdb_batch_txt,
			gtdbtk_data_tar_gz = gtdbtk_data_tar_gz,
			derep_bins= derep_bins,
			runtime_attributes = default_runtime_attributes
	}

	call gtdbtk_cleanup {
		input:
			sample_id = sample_id,
			gtdbtk_classify_tar_gz = gtdbtk_analysis.gtdbtk_classify_tar_gz,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File gtdbtk_align_tar_gz = gtdbtk_analysis.gtdbtk_align_tar_gz
		File gtdbtk_classify_tar_gz = gtdbtk_analysis.gtdbtk_classify_tar_gz
		File gtdbtk_identify_tar_gz = gtdbtk_analysis.gtdbtk_identify_tar_gz
		File gtdbk_summary_txt = gtdbtk_cleanup.gtdbk_summary_txt
	}
}

task gtdbtk_analysis {
	input {
		String sample_id
		
		File gtdb_batch_txt
		File gtdbtk_data_tar_gz
		Array[File] derep_bins

		RuntimeAttributes runtime_attributes
	}

	Int threads = 48
	Int mem_gb = threads * 2
	Int disk_size = ceil((size(gtdbtk_data_tar_gz, "GB") + (size(derep_bins[0], "GB") * length(derep_bins))) * 2 + 20)

	command <<<
		set -euo pipefail

		mkdir gtdbtk_out_dir
		mkdir tmp_dir

		tar -xzvf ~{gtdbtk_data_tar_gz}

		# Must set $GTDBTK_DATA_PATH variable to use gtdbtk command
		GTDBTK_DATA_PATH="$(pwd)/release207_v2"
		export GTDBTK_DATA_PATH

		gtdbtk --version

		gtdbtk classify_wf \
			--batchfile ~{gtdb_batch_txt} \
			--out_dir gtdbtk_out_dir \
			--extension fa \
			--prefix ~{sample_id} \
			--cpus ~{threads} \
			--tmpdir tmp_dir

		tar -C gtdbtk_out_dir -czvf "~{sample_id}.align.tar.gz" align
		tar -C gtdbtk_out_dir -czvf "~{sample_id}.classify.tar.gz" classify
		tar -C gtdbtk_out_dir -czvf "~{sample_id}.identify.tar.gz" identify
	>>>

	output {
		File gtdbtk_align_tar_gz = "~{sample_id}.align.tar.gz"
		File gtdbtk_classify_tar_gz = "~{sample_id}.classify.tar.gz"
		File gtdbtk_identify_tar_gz = "~{sample_id}.identify.tar.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/gtdbtk:5e8307c"
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

task gtdbtk_cleanup {
	input {
		String sample_id
		
		File gtdbtk_classify_tar_gz

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(gtdbtk_classify_tar_gz, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		classify_dir=$(tar -tzf ~{gtdbtk_classify_tar_gz} | head -1 | cut -d '/' -f 1)
		tar -xzvf ~{gtdbtk_classify_tar_gz}

		python /opt/scripts/GTDBTk-Organize.py \
			--input_dir "${classify_dir}" \
			--outfile "~{sample_id}.GTDBTk_Summary.txt"
	>>>

	output {
		File gtdbk_summary_txt = "~{sample_id}.GTDBTk_Summary.txt"
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
