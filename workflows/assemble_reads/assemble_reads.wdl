version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow assemble_reads {
	input {
		String sample_id
		File hifi_reads

		RuntimeAttributes default_runtime_attributes
	}

	# Detect the input format of the hifi_reads; if BAM, first convert to FASTQ
	String hifi_reads_extension = sub(basename(hifi_reads), basename(hifi_reads, ".bam"), "")

	if (hifi_reads_extension == ".bam") {
		call bam_to_fastq {
			input:
				bam = hifi_reads,
				runtime_attributes = default_runtime_attributes
		}
	}

	call assemble_metagenomes {
		input:
			sample_id = sample_id,
			fastq = select_first([bam_to_fastq.fastq, hifi_reads]),
			runtime_attributes = default_runtime_attributes
	}

	output {
		File? fastq = bam_to_fastq.fastq

		File primary_contig_gfa = assemble_metagenomes.primary_contig_gfa
		File primary_contig_fasta = assemble_metagenomes.primary_contig_fasta
		File primary_contig_fasta_gz = assemble_metagenomes.primary_contig_fasta_gz
	}
}

task bam_to_fastq {
	input {
		File bam

		RuntimeAttributes runtime_attributes
	}

	String bam_basename = basename(bam, ".bam")
	Int threads = 2
	Int disk_size = ceil(size(bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		samtools --version

		samtools fastq \
			-@ ~{threads - 1} \
			~{bam} \
		| bgzip \
			--stdout \
		> "~{bam_basename}.fastq.gz"
	>>>

	output {
		File fastq = "~{bam_basename}.fastq.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools:5e8307c"
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

task assemble_metagenomes {
	input {
		String sample_id
		File fastq

		RuntimeAttributes runtime_attributes
	}

	Int threads = 32
	Int mem_gb = threads * 4
	Int disk_size = ceil(size(fastq, "GB") * 4 + 20)

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
		File primary_contig_gfa = "~{sample_id}.asm.p_ctg.gfa"
		File primary_contig_fasta = "~{sample_id}.asm.p_ctg.fa"
		File primary_contig_fasta_gz = "~{sample_id}.asm.p_ctg.fa.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/hifiasm-meta:0.3.1"
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
