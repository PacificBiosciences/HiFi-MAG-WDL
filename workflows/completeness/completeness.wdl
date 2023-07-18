version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow completeness {
	input {
		String sample
		File contig

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
			contig = contig,
			min_length = min_length,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File key = long_contigs_to_bins.key
	}

	parameter_meta {
		sample: {help: "Sample name"}
		contig: {help: "Contigs"} #TODO
		min_length: {help: "Minimum size of a contig to consider for completeness scores; default value is set to 500kb. This value should not be increased"}
		min_completeness: {help: "Minimum completeness score (from CheckM2) to mark a contig as complete and place it in a distinct bin; default value is set to 93%. This value should not be lower than 90%"}
		key: {help: ""} #TODO
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task long_contigs_to_bins {
	input {
		String sample
		File contig

		Int min_length

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contig, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Fasta-Make-Long-Seq-Bins.py \
			-i ~{contig} \
			-b "~{sample}.bin_key.txt" \
			-l ~{min_length}
			-o ./ #TODO
	>>>

	output {
		File key = "~{sample}.bin_key.txt"
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
