version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow completeness_aware_binning {
	input {
		String sample_id
		File contigs_fasta
		File db

		#TODO - conditional for these values?
		Int min_length = 500000
		Int min_completeness = 93

		RuntimeAttributes default_runtime_attributes
	}

	call long_contigs_to_bins {
		input:
			sample_id = sample_id,
			contigs_fasta = contigs_fasta,
			min_length = min_length,
			runtime_attributes = default_runtime_attributes
	}

	call make_incomplete_contigs {
		input:
			sample_id = sample_id,
			contigs_fasta = contigs_fasta,
			key = long_contigs_to_bins.key,
			long_bin_fastas = long_contigs_to_bins.long_bin_fastas,
			runtime_attributes = default_runtime_attributes
	}

	#if (long_contigs_to_bins.bin_key_nonempty) {
	#	call checkm2_contig_analysis {
	#	input:
	#		db = db,
	#		long_bin_fastas = long_contigs_to_bins.long_bin_fastas,
	#		runtime_attributes = default_runtime_attributes
	#	}

	#	call filter_complete_contigs {
	#	input:
	#		sample_id = sample_id,
	#		contigs_fasta = contigs_fasta,
	#		report = checkm2_contig_analysis.report,
	#		key = long_contigs_to_bins.key,
	#		min_length = min_length,
	#		min_completeness = min_completeness,
	#		runtime_attributes = default_runtime_attributes
	#	}
	#}

	output {
		File key = long_contigs_to_bins.key
		File incomplete_contigs = make_incomplete_contigs.incomplete_contigs

		#File? report = checkm2_contig_analysis.report
		#File? passed_bins = filter_complete_contigs.passed_bins
		#File? scatterplot = filter_complete_contigs.scatterplot
		#File? histogram = filter_complete_contigs.histogram
	}
}

task long_contigs_to_bins {
	input {
		String sample_id
		File contigs_fasta

		Int min_length

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contigs_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Fasta-Make-Long-Seq-Bins.py \
			-i ~{contigs_fasta} \
			-b "~{sample_id}.bin_key.txt" \
			-l ~{min_length} \
			-o ./

		# Check if any long contigs (>500kb) were identified
		if [[ -s "~{sample_id}.bin_key.txt" ]]; then
			echo "true" > bin_key_nonempty.txt
		else
			echo "false" > bin_key_nonempty.txt
		fi
	>>>

	output {
		File key = "~{sample_id}.bin_key.txt"
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
		String sample_id
		File contigs_fasta
		
		File key
		Array[File] long_bin_fastas

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contigs_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		long_bin_fastas_dir=$(dirname ~{long_bin_fastas[0]})

		# TODO - need to rename long bin fastas here
		while IFS= read -r fasta; do
			mv "$(echo "$fasta" | awk '{print $1}' | awk '{print "'"$long_bin_fastas_dir"'/"$1".fa"}')" "$(echo "$fasta" | awk '{print $2}' | awk '{print "'"$long_bin_fastas_dir"'/"$1".fa"}')"
		done < ~{key}

		python /opt/scripts/Make-Incomplete-Contigs.py \
			-i ~{contigs_fasta} \
			-f "~{sample_id}.incomplete_contigs.fasta" \
			-p ~{key} \
			-d "${long_bin_fastas_dir}" \
			-o ./ # TODO - this step is just copying the complete AKA long bin fastas to an output directory
	>>>

	output {
		File incomplete_contigs = "~{sample_id}.incomplete_contigs.fasta"
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

	Int threads = 24
	Int mem_gb = threads * 4
	Int disk_size = ceil(size(db, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		long_bin_fastas_dir=$(dirname ~{long_bin_fastas[0]})

		checkm2 predict \
			-i "$long_bin_fastas_dir" \
			-o ./ \
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

task filter_complete_contigs {
	input {
		String sample_id
		File contigs_fasta

		File report
		File key

		Int min_length
		Int min_completeness

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contigs_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Filter-Complete-Contigs.py \
			-i ~{contigs_fasta} \
			-c ~{report} \
			-b ~{key} \
			-l ~{min_length} \
			-m ~{min_completeness} \
			-p "~{sample_id}.passed_bins.txt" \
			-p1 "~{sample_id}.completeness_vs_size_scatter.pdf" \
			-p2 "~{sample_id}.completeness_histo.pdf"
	>>>

	output {
		File passed_bins = "~{sample_id}.passed_bins.txt"
		File scatterplot = "~{sample_id}.completeness_vs_size_scatter.pdf"
		File histogram = "~{sample_id}.completeness_histo.pdf"
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
