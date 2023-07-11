version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow completeness {
	input {
		String sample
		Array[File] contigs

		#TODO - conditional for these values?
		Int min_length = 500000
		Int min_completeness = 93

		String outdir

		RuntimeAttributes default_runtime_attributes
	}

	#TODO - use call or conditional?
	#Checkpoint 1: Fork samples into two groups: those with long contigs (> min length), and those with no contigs.
	call long_contigs_to_bins {
		input:
			sample = sample,
			contigs = contigs,
			min_length = min_length,
			outdir = outdir,
			runtime_attributes = default_runtime_attributes
	}

	call close_long_bin_fork {
		input:
			sample = sample,
			contigs = contigs,
			key = long_contigs_to_bins.key,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File key = long_contigs_to_bins.key
	}

	parameter_meta {
		sample: {help: "Sample name"}
		contigs: {help: "Contigs"} #TODO
		min_length: {help: "Minimum size of a contig to consider for completeness scores; default value is set to 500kb. This value should not be increased"}
		min_completeness: {help: "Minimum completeness score (from CheckM2) to mark a contig as complete and place it in a distinct bin; default value is set to 93. This value should not be lower than 90"}
		outdir: {help: ""} #TODO
		key: {help: ""} #TODO
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task long_contigs_to_bins {
	input {
		String sample
		Array[File] contigs

		Int min_length

		String outdir

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contigs, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		Fasta-Make-Long-Seq-Bins.py \
			-i ~{sep=' ' contigs} \ #TODO - check if need to avoid using sep due to 'argument list too long' error
			-o ~{outdir} \
			-b "~{sample}.bin_key.txt" \
			-l ~{min_length}
	>>>

	output {
		File key = "~{sample}.bin_key.txt"
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

task close_long_bin_fork {
	input {
		String sample
		Array[File] contigs
		
		File key

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contigs, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		#TODO - add this function to get targets of forked workflow:
		#https://github.com/PacificBiosciences/pb-metagenomics-tools/blob/5e8307ccaf9b5f8b5b40421a4a39a7a9165af01b/HiFi-MAG-Pipeline/Snakefile-hifimags.smk#L66-L73

		Make-Incomplete-Contigs.py \
			-i ~{sep=' ' contigs} \ #TODO - check if need to avoid using sep due to 'argument list too long' error
			-f "~{sample}.incomplete_contigs.fasta" \ #TODO - is this an array?
			-p {input} \ #TODO - function
			-d {params.fastadir} \ #TODO
			-o {output.outdir} && cp {input} {output.complete} #TODO
	>>>

	output {
		File incomplete_contigs = "~{sample}.incomplete_contigs.fasta"
		File complete = "~{sample}.LongBinCompleted.txt"
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




