version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow qc {
	input {
		Array[File] contigs
		String path

		RuntimeAttributes default_runtime_attributes
	}

	call checkm2 {
		input:
			contigs = contigs,
			path = path,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File key = checkm2.key
		File complete = checkm2.complete
	}

	parameter_meta {
		contigs: {help: "Contigs"} #TODO
		path: {help: "Path to install CheckM2"} #TODO
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

#TODO - Setup CheckM2; put in docker image?
task checkm2 {
	input {
		Array[File] contigs
		String path

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24 #TODO
	Int disk_size = ceil(size(contigs, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		checkm2 database \
			--download \
			--path ~{path}
	>>>

	output {
		File key = "uniref100.KO.1.dmnd"
		File complete = "CheckM2.complete.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/" #TODO
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
