version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow assemble_metagenomes {
	input {
		String sample_id
		File? hifi_reads_bam
		File? hifi_reads_fastq

		RuntimeAttributes default_runtime_attributes
	}

	if (! defined(hifi_reads_fastq)) {
		call bam_to_fastq {
			input:
				hifi_reads_bam = select_first([hifi_reads_bam]),
				runtime_attributes = default_runtime_attributes
		}
	}

	call assemble_reads {
		input:
			sample_id = sample_id,
			fastq = select_first([bam_to_fastq.converted_fastq, hifi_reads_fastq]),
			runtime_attributes = default_runtime_attributes
	}

	output {
		File? converted_fastq = bam_to_fastq.converted_fastq

		File assembled_contigs_gfa = assemble_reads.assembled_contigs_gfa
		File assembled_contigs_fa = assemble_reads.assembled_contigs_fa
		File assembled_contigs_fa_gz = assemble_reads.assembled_contigs_fa_gz
	}
}

task bam_to_fastq {
	input {
		File hifi_reads_bam

		RuntimeAttributes runtime_attributes
	}

	String bam_basename = basename(hifi_reads_bam, ".bam")
	Int threads = 2
	Int disk_size = ceil(size(hifi_reads_bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		samtools --version

		samtools fastq \
			-@ ~{threads - 1} \
			~{hifi_reads_bam} \
		| bgzip \
			--stdout \
		> "~{bam_basename}.fastq.gz"
	>>>

	output {
		File converted_fastq = "~{bam_basename}.fastq.gz"
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

task assemble_reads {
	input {
		String sample_id
		File fastq

		RuntimeAttributes runtime_attributes
	}

	Int threads = 48
	Int mem_gb = ceil(size(fastq, "GB") / 0.07 + 20)
	Int disk_size = ceil(size(fastq, "GB") * 8 + 20)

	command <<<
		set -euo pipefail

		hifiasm_meta --version
		gfatools version

		hifiasm_meta \
			-t ~{threads} \
			-o ~{sample_id}.asm \
			~{fastq}

		gfatools gfa2fa \
			~{sample_id}.asm.p_ctg.gfa \
		> ~{sample_id}.asm.p_ctg.fa

		bgzip \
			--threads ~{threads} \
			--stdout \
			~{sample_id}.asm.p_ctg.fa \
		> ~{sample_id}.asm.p_ctg.fa.gz
	>>>

	output {
		File assembled_contigs_gfa = "~{sample_id}.asm.p_ctg.gfa"
		File assembled_contigs_fa = "~{sample_id}.asm.p_ctg.fa"
		File assembled_contigs_fa_gz = "~{sample_id}.asm.p_ctg.fa.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/hifiasm-meta@sha256:13d4378116d137f7b3bfe2af4555866bfe2a52085e48a4c006d804757b7a5cd7"
		cpu: threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: if (size(fastq, "GB") < 10) then 0 else runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
