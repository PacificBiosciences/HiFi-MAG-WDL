version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow depth {
	input {
		String sample

		File contigs
		File reads

		RuntimeAttributes default_runtime_attributes
	}

	call minimap_to_bam {
		input:
			sample = sample,
			contigs = contigs,
			reads = reads,
			runtime_attributes = default_runtime_attributes
	}

	call bam_depth {
		input:
			sample = sample,
			bam = minimap_to_bam.bam,
			runtime_attributes = default_runtime_attributes
	}

	call convert_depth {
		input:
			sample = sample,
			depth = bam_depth.depth,
			passed = passed,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File bam = minimap_to_bam.bam
		File bam_index = minimap_to_bam.bam_index
		File complete = minimap_to_bam.complete
		File filtered_depth = convert_depth.filtered_depth
	}

	parameter_meta {
		sample: {help: "Sample name"}
		contigs: {help: "Contigs"} #TODO
		reads: {help: "Fasta file containing sample reads"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task minimap_to_bam {
	input {
		String sample
		
		File reads
		File contigs

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int disk_size = ceil(size(reads, "GB") * 2 + size(contigs, "GB") + 20)

	command <<<
		set -euo pipefail

		minimap2 \
			-a \
			-k 19 \
			-w 10 \
			-I 10G \
			-g 5000 \
			-r 2000 \
			--lj-min-ratio 0.5 \
			-A 2 \
			-B 5 \
			-O 5,56 \
			-E 4,1 \
			-z 400,50 \
			--sam-hit-only \
			-t ~{threads} \
			~{contigs} \
			~{reads} \
		| samtools sort \
			-@ ~{threads} \
			-o "~{sample}.bam"

		samtools index \
			-@ ~{threads} \
			"~{sample}.bam"
	>>>

	output {
		File bam = "~{sample}.bam"
		File bam_index = "~{sample}.bam.bai"
		File complete = "~{sample}.index.completed.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools@sha256:d3f5cfb7be7a175fe0471152528e8175ad2b57c348bacc8b97be818f31241837"
		cpu: threads
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

task bam_depth {
	input {
		String sample
		
		File bam

		RuntimeAttributes runtime_attributes
	}
 
	Int disk_size = ceil(size(bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		jgi_summarize_bam_contig_depths \
			--outputDepth "~{sample}.JGI.depth.txt" \
			~{bam}
	>>>

	output {
		File depth = "~{sample}.JGI.depth.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/metabat@sha256:d68d401803c4e99d85a4c93e48eb223275171f08b49d3314ff245a2d68264651"
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

task convert_depth {
	input {
		String sample
		
		File depth
		File passed

		RuntimeAttributes runtime_attributes
	}
 
	Int disk_size = ceil(size(depth, "GB") * 2 + size(passed, "GB") + 20)

	command <<<
		set -euo pipefail

		Convert-JGI-Coverages.py \
			-i ~{depth} \
			-p ~{passed} \
			-o1 "~{sample}.JGI.filtered.depth.txt"
	>>>

	output {
		File filtered_depth = "~{sample}.JGI.filtered.depth.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python@sha256:f4bf8adfa4987cf61d037440ec6f7fc1cff80c33adbe620620579f1e1f9c08f1"
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
