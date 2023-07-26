version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow gtdbtk {
	input {
		String sample_id

		File gtdb_batch_txt
		String gtdbtk_data_path

		RuntimeAttributes default_runtime_attributes
	}

	call gtdbtk_analysis {
		input:
			sample_id = sample_id,
			gtdb_batch_txt = gtdb_batch_txt,
			gtdbtk_data_path = gtdbtk_data_path,
			runtime_attributes = default_runtime_attributes
	}

	call gtdbtk_cleanup {
		input:
			sample_id = sample_id,
			classify = gtdbtk_analysis.classify,
			runtime_attributes = default_runtime_attributes
	}

	output {
		Array[File] all = gtdbtk_analysis.all
		File gtdbk_summary_txt = gtdbtk_cleanup.gtdbk_summary_txt
	}
}

task gtdbtk_analysis {
	input {
		String sample_id
		
		File gtdb_batch_txt
		String gtdbtk_data_path

		RuntimeAttributes runtime_attributes
	}

	Int threads = 48
	Int mem_gb = threads * 2
	Int disk_size = ceil(size(gtdb_batch_txt, "GB") * 2 + 60)

	command <<<
		set -euo pipefail

		gtdbtk --version

		mkdir gtdbtk_out_dir

		# Must set $GTDBTK_DATA_PATH variable to use gtdbtk command
		GTDBTK_DATA_PATH=~{gtdbtk_data_path}

		gtdbtk classify_wf \
			--batchfile ~{gtdb_batch_txt} \
			--out_dir gtdbtk_out_dir \
			-x fa \
			--prefix ~{sample_id} \
			--cpus ~{threads}
	>>>

	output {
		# TODO
		#Array[File] align = glob("gtdbtk_out_dir/")
		#Array[File] classify = glob("gtdbtk_out_dir/")
		#Array[File] identify = glob("gtdbtk_out_dir/")
		Array[File] all = glob("gtdbtk_out_dir/*")
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
		
		Array[File] classify

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(classify, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		classify_dir=$(dirname ~{classify[0]})

		python /opt/scripts/GTDBTk-Organize.py \
			--input_dir "$classify_dir" \
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
