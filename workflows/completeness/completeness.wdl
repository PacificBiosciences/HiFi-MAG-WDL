version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow completeness {
	input {
		String sample
		Array[File] contigs

		#TODO - conditional for these values?
		Int min_length = 500000
		#Int min_completeness = 93

		RuntimeAttributes default_runtime_attributes
	}

	#TODO - use call or conditional?
	#Checkpoint 1: Fork samples into two groups: those with long contigs (> min length), and those with no contigs.
	## Write fasta seqs to individual fasta files if longer than length threshhold
	call long_contigs_to_bins {
		input:
			sample = sample,
			contigs = contigs,
			min_length = min_length,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File key = long_contigs_to_bins.key
	}

	parameter_meta {
		sample: {help: "Sample name"}
		contigs: {help: "Contigs"} #TODO
		min_length: {help: "Minimum size of a contig to consider for completeness scores; default value is set to 500kb. This value should not be increased"}
		min_completeness: {help: "Minimum completeness score (from CheckM2) to mark a contig as complete and place it in a distinct bin; default value is set to 93%. This value should not be lower than 90%"}
		key: {help: ""} #TODO
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task long_contigs_to_bins {
	input {
		String sample
		Array[File] contigs

		Int min_length

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contigs, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		Fasta-Make-Long-Seq-Bins.py \
			-i ~{sep=' ' contigs} \ #TODO - check if need to avoid using sep due to 'argument list too long' error
			-b "~{sample}.bin_key.txt" \
			-l ~{min_length}
	>>>

	output {
		File key = "~{sample}.bin_key.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python@sha256:e4d921e252c3c19fe64097aa619c369c50cc862768d5fcb5e19d2877c55cfdd2"
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
