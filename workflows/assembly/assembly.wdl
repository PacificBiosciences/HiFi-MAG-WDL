version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow assembly {
	input {
		String sample
		File fastq

		RuntimeAttributes default_runtime_attributes
	}

	call hifiasm_meta {
		input:
			sample = sample,
			fastq = fastq,
			runtime_attributes = default_runtime_attributes
	}

	output {
		 
	}

	parameter_meta {
		sample: {help: "Sample name"}
		fastq: {help: "Sample"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task hifiasm_meta {
	input {
		String sample
		File fastq

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(fastq, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		hifiasm_meta \
			-t32 \
			-oasm \
			~{fastq}
	>>>

	output {
		 
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/hifiasm-meta@sha256:c513c314af775a879ebb642a9c959a3b13d92a82ccda2dd514abf2139404a481"
		cpu: 4
		memory: "32 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
