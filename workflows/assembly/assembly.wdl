version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow assembly {
	input {
		String sample_id
		File fastq

		RuntimeAttributes default_runtime_attributes
	}

	call hifiasm_meta {
		input:
			sample_id = sample_id,
			fastq = fastq,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File primary_contig_graph = hifiasm_meta.primary_contig_graph
		File primary_contig_fasta = hifiasm_meta.primary_contig_fasta
		File reads_fasta = hifiasm_meta.reads_fasta
	}

	parameter_meta {
		sample_id: {help: "Sample ID"}
		fastq: {help: "Sample in fastq format"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task hifiasm_meta {
	input {
		String sample_id
		File fastq

		RuntimeAttributes runtime_attributes
	}

	Int threads = 32
	Int mem_gb = threads * 2
	Int disk_size = ceil(size(fastq, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		hifiasm_meta \
			-t ~{threads} \
			-o ~{sample_id} \
			~{fastq}

		awk '/^S/{print ">"$2;print $3}' "~{sample_id}.p_ctg.gfa" > "~{sample_id}.p_ctg.fa"
	>>>

	output {
		File primary_contig_graph = "~{sample_id}.p_ctg.gfa"
		File primary_contig_fasta = "~{sample_id}.p_ctg.fa"
		File reads_fasta = "~{sample_id}.rescue.fa"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/hifiasm-meta@sha256:c513c314af775a879ebb642a9c959a3b13d92a82ccda2dd514abf2139404a481"
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
