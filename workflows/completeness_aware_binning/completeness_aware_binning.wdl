version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow completeness_aware_binning {
	input {
		String sample
		File contig
		File db

		#TODO - conditional for these values?
		Int min_length = 500000
		Int min_completeness = 93

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

	# Checkpoint 1 aggregator; close the fork; outputs '/1-long-contigs/SAMPLE/SAMPLE.incomplete_contigs.fasta'
	## Identify long complete contigs from checkm2 scores
	call make_incomplete_contigs {
		input:
			sample = sample,
			contig = contig,
			key = long_contigs_to_bins.key,
			long_bin_fastas = long_contigs_to_bins.long_bin_fastas,
			runtime_attributes = default_runtime_attributes
	}

	if (long_contigs_to_bins.bin_key_nonempty) {
		call checkm2_contig_analysis {
		input:
			db = db,
			long_bin_fastas = long_contigs_to_bins.long_bin_fastas,
			runtime_attributes = default_runtime_attributes
		}

		call filter_complete_contigs {
		input:
			sample = sample,
			contig = contig,
			report = checkm2_contig_analysis.report,
			key = long_contigs_to_bins.key,
			min_length = min_length,
			min_completeness = min_completeness,
			runtime_attributes = default_runtime_attributes
		}
	}

	output {
		File key = long_contigs_to_bins.key
		File incomplete_contigs = make_incomplete_contigs.incomplete_contigs

		File? report = checkm2_contig_analysis.report
		File? passed = filter_complete_contigs.passed
		File? p1 = filter_complete_contigs.p1
		File? p2 = filter_complete_contigs.p2
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
			-l ~{min_length} \
			-o ./ #TODO

		# Check if any long contigs (>=500kb) were identified
		if [[ -s "~{sample}.bin_key.txt" ]]; then
			echo "true" > bin_key_nonempty.txt
		else
			echo "false" > bin_key_nonempty.txt
		fi
	>>>

	output {
		File key = "~{sample}.bin_key.txt"
		Boolean bin_key_nonempty = read_boolean("bin_key_nonempty.txt")
		Array[File] long_bin_fastas = glob("complete.*.fa")
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

task make_incomplete_contigs {
	input {
		String sample
		File contig
		
		File key
		Array[File] long_bin_fastas

		RuntimeAttributes runtime_attributes
	}

	String long_bin_fastas_dir = sub(long_bin_fastas[0], "complete.1.fa", "")
	Int disk_size = ceil(size(contig, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		# TODO - need to rename long bin fastas here

		python /opt/scripts/Make-Incomplete-Contigs.py \
			-i ~{contig} \
			-f "~{sample}.incomplete_contigs.fasta" \
			-p ~{key} \
			-d ~{long_bin_fastas_dir} \
			-o ./ \ # TODO - this step is just copying the complete AKA long bin fastas to an output directory
	>>>

	output {
		File incomplete_contigs = "~{sample}.incomplete_contigs.fasta"
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

task checkm2_contig_analysis {
	input {
		File db

		Array[File] long_bin_fastas

		RuntimeAttributes runtime_attributes
	}

	String long_bin_fastas_dir = sub(long_bin_fastas[0], "complete.1.fa", "")
	Int threads = 24
	Int disk_size = ceil(size(db, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		checkm2 predict \
			-i ~{long_bin_fastas_dir} \
			-o ./ \ # TODO
			-x fa \
			-t ~{threads} \
			--force \
			--database_path ~{db} \
			--remove_intermediates
	>>>

	output {
		File report = "quality_report.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/" #TODO
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

task filter_complete_contigs {
	input {
		String sample
		File contig

		File report
		File key

		Int min_length
		Int min_completeness

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contig, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Filter-Complete-Contigs.py \
			-i ~{contig} \
			-c ~{report} \
			-b ~{key} \
			-l ~{min_length} \
			-m ~{min_completeness} \
			-p "~{sample}.passed_bins.txt" \
			-p1 "~{sample}.completeness_vs_size_scatter.pdf" \
			-p2 "~{sample}.completeness_histo.pdf"
	>>>

	output {
		File passed = "~{sample}.passed_bins.txt"
		File p1 = "~{sample}.completeness_vs_size_scatter.pdf"
		File p2 = "~{sample}.completeness_histo.pdf"
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
