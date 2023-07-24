version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow coverage {
	input {
		String sample_id

		File contigs_fasta
		File reads_fasta
		File key

		RuntimeAttributes default_runtime_attributes
	}

	call minimap_to_bam {
		input:
			sample_id = sample_id,
			contigs_fasta = contigs_fasta,
			reads_fasta = reads_fasta,
			runtime_attributes = default_runtime_attributes
	}

	call bam_depth {
		input:
			sample_id = sample_id,
			sorted_bam = minimap_to_bam.sorted_bam,
			runtime_attributes = default_runtime_attributes
	}

	call convert_depth {
		input:
			sample_id = sample_id,
			depth = bam_depth.depth,
			key = key,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File sorted_bam = minimap_to_bam.sorted_bam
		File sorted_bam_index = minimap_to_bam.sorted_bam_index
		File filtered_depth = convert_depth.filtered_depth
	}
}

task minimap_to_bam {
	input {
		String sample_id
		
		File contigs_fasta
		File reads_fasta

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = threads * 2
	Int disk_size = ceil(size(contigs_fasta, "GB") * 2 + size(reads_fasta, "GB") + 20)

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
			~{contigs_fasta} \
			~{reads_fasta} \
		| samtools sort \
			-@ ~{threads} \
			-o "~{sample_id}.sorted.bam"

		samtools index \
			-@ ~{threads} \
			"~{sample_id}.sorted.bam"
	>>>

	output {
		File sorted_bam = "~{sample_id}.sorted.bam"
		File sorted_bam_index = "~{sample_id}.sorted.bam.bai"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools@sha256:d3f5cfb7be7a175fe0471152528e8175ad2b57c348bacc8b97be818f31241837"
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

task bam_depth {
	input {
		String sample_id
		
		File sorted_bam

		RuntimeAttributes runtime_attributes
	}
 
	Int disk_size = ceil(size(sorted_bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		jgi_summarize_bam_contig_depths \
			--outputDepth "~{sample_id}.JGI.depth.txt" \
			~{sorted_bam}
	>>>

	output {
		File depth = "~{sample_id}.JGI.depth.txt"
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
		String sample_id
		
		File depth
		File key

		RuntimeAttributes runtime_attributes
	}
 
	Int disk_size = ceil(size(depth, "GB") * 2 + size(key, "GB") + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Convert-JGI-Coverages.py \
			-i ~{depth} \
			-p ~{key} \
			-o1 "~{sample_id}.JGI.filtered.depth.txt"
	>>>

	output {
		File filtered_depth = "~{sample_id}.JGI.filtered.depth.txt"
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
