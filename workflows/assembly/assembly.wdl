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
		File primary_contig_graph = hifiasm_meta.primary_contig_graph
		File primary_contig_fasta = hifiasm_meta.primary_contig_fasta
		File reads_fasta = hifiasm_meta.reads_fasta
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
			-o ~{sample} \
			~{fastq}

		awk '/^S/{print ">"$2;print $3}' "~{sample}.p_ctg.gfa" > "~{sample}.p_ctg.fa"
	>>>

	output {
		File bins_tsv = "~{sample}.bins.tsv"
		File cmd = "{sample}.cmd"
		File primary_contig_graph = "~{sample}.p_ctg.gfa"
		File primary_contig_graph_noseq = "~{sample}.p_ctg.noseq.gfa"
		File primary_contig_fasta = "~{sample}.p_ctg.fa"
		File reads_fasta_all = "~{sample}.rescue.all.fa"
		File reads_fasta = "~{sample}.rescue.fa"

		Array[File] raw_unitig_graphs = glob("~{sample}.r_utg*.gfa")
		Array[File] unitig_graphs = glob("~{sample}.p_utg*.gfa")
		Array[File] alternate_contig_graphs = glob("~{sample}.a_ctg*.gfa")
		Array[File] binary_files = glob("*.bin")
		Array[File] tangle_unitig_graphs = glob("~{sample}.GresolveTangle_before_*.utg.gfa")
		Array[File] cln_unitig_graphs = glob("~{sample}.GHEYbefore_circle_cln_*.utg.gfa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/hifiasm-meta@sha256:c513c314af775a879ebb642a9c959a3b13d92a82ccda2dd514abf2139404a481"
		cpu: 2
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
