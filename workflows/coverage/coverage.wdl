version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow coverage {
	input {
		String sample_id

		File contigs_fasta_gz
		File hifi_reads_fastq
		File bins_contigs_key_txt

		RuntimeAttributes default_runtime_attributes
	}

	call align_hifiasm {
		input:
			sample_id = sample_id,
			contigs_fasta_gz = contigs_fasta_gz,
			hifi_reads_fastq = hifi_reads_fastq,
			runtime_attributes = default_runtime_attributes
	}

	call jgi_bam_depth {
		input:
			sample_id = sample_id,
			aligned_sorted_bam = align_hifiasm.aligned_sorted_bam,
			aligned_sorted_bam_index = align_hifiasm.aligned_sorted_bam_index,
			runtime_attributes = default_runtime_attributes
	}

	call convert_jgi_bamdepth {
		input:
			sample_id = sample_id,
			contig_depth_txt = jgi_bam_depth.contig_depth_txt,
			bins_contigs_key_txt = bins_contigs_key_txt,
			runtime_attributes = default_runtime_attributes
	}

	output {
		IndexData aligned_sorted_bam = {
			"data": align_hifiasm.aligned_sorted_bam,
			"data_index": align_hifiasm.aligned_sorted_bam_index
		}
		File filtered_contig_depth_txt = convert_jgi_bamdepth.filtered_contig_depth_txt
	}
}

task align_hifiasm {
	input {
		String sample_id

		File contigs_fasta_gz
		File hifi_reads_fastq

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = threads * 2
	Int disk_size = ceil((size(contigs_fasta_gz, "GB") + size(hifi_reads_fastq, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		minimap2 --version
		samtools --version

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
			-t ~{threads / 2} \
			~{contigs_fasta_gz} \
			~{hifi_reads_fastq} \
		| samtools sort \
			-@ ~{threads / 2 - 1} \
			-o "~{sample_id}.aligned.sorted.bam"

		samtools index \
			-@ ~{threads - 1} \
			"~{sample_id}.aligned.sorted.bam"
	>>>

	output {
		File aligned_sorted_bam = "~{sample_id}.aligned.sorted.bam"
		File aligned_sorted_bam_index = "~{sample_id}.aligned.sorted.bam.bai"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools:5e8307c"
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

task jgi_bam_depth {
	input {
		String sample_id

		File aligned_sorted_bam
		File aligned_sorted_bam_index

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(aligned_sorted_bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		metabat --help |& grep "version"

		jgi_summarize_bam_contig_depths \
			--outputDepth "~{sample_id}.JGI.depth.txt" \
			~{aligned_sorted_bam}
	>>>

	output {
		File contig_depth_txt = "~{sample_id}.JGI.depth.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/metabat:5e8307c"
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

task convert_jgi_bamdepth {
	input {
		String sample_id

		File contig_depth_txt
		File bins_contigs_key_txt

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(contig_depth_txt, "GB") + size(bins_contigs_key_txt, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Convert-JGI-Coverages.py \
			--in_jgi ~{contig_depth_txt} \
			--passed_bins ~{bins_contigs_key_txt} \
			--out_jgi "~{sample_id}.JGI.filtered.depth.txt"
	>>>

	output {
		File filtered_contig_depth_txt = "~{sample_id}.JGI.filtered.depth.txt"
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
