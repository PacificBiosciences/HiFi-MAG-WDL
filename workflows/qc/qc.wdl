version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow qc {
	input {
		Array[File] contigs

		RuntimeAttributes default_runtime_attributes
	}

	call checkm2 {
		input:
			contigs = contigs,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File key = checkm2.key
		File complete = checkm2.complete
	}

	parameter_meta {
		
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task checkm2 {
	input {
		Array[File] contigs

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int disk_size = ceil(size(contigs, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		checkm2 database \
			--download \
			--path {params.installdir} &> {log} && touch {output.complete}
	>>>

	output {
		File key = "~{gvcf_basename}.expanded.g.vcf.gz"
		File complete = "~{gvcf_basename}.expanded.g.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/"
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
