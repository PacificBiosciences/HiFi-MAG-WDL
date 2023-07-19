version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow bam_to_fastq {
	input {
		String sample
		File bam

		RuntimeAttributes default_runtime_attributes
	}

	call convert_bam_to_fastq {
		input:
			sample = sample,
			bam = bam,
			runtime_attributes = default_runtime_attributes
	}

	output {
		 File converted_fastq = convert_bam_to_fastq.converted_fastq
	}

	parameter_meta {
		sample: {help: "Sample name"}
		bam: {help: "Sample BAM to convert"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task convert_bam_to_fastq {
	input {
		String sample
		File bam

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		samtools fastq \
			~{bam} \
		| bgzip -c >"~{sample}.fastq.gz"
	>>>

	output {
		 File converted_fastq = "~{sample}.fastq.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools@sha256:d3f5cfb7be7a175fe0471152528e8175ad2b57c348bacc8b97be818f31241837"
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
