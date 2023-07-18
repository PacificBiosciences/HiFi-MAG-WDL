version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow checkm2 {
	input {
		String sample

		File db
		File dbrep_bins

		File depth

		Int min_completeness = 70
		Int max_contamination = 10
		Int max_contigs = 20

		RuntimeAttributes default_runtime_attributes
	}

	#CheckM2 analysis and checkpoint, fork to GTDB-Tk & summary (if bins found) or end analysis (no bins found)
	call checkm2_bin_analysis {
		input:
			sample = sample,
			db = db,
			dbrep_bins = dbrep_bins,
			runtime_attributes = default_runtime_attributes
	}

	#Checkpoint 2 - Ensure there are bins after CheckM2, before running GTDB-Tk and the summary
	call assess_checkm2_bins {
		input:
			sample = sample,
			qv = checkm2_bin_analysis.qv,
			depth = depth,
			min_completeness = min_completeness,
			max_contamination = max_contamination,
			max_contigs = max_contigs,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File qv = checkm2_bin_analysis.qv
		File gtdb = assess_checkm2_bins.gtdb
		File target = assess_checkm2_bins.target
		File tsv = assess_checkm2_bins.tsv
	}

	parameter_meta {
		sample: {help: "Sample name"}
		db: {help: "DIAMOND database CheckM2 relies on for rapid annotation"}
		dbrep_bins: {help: "Dereplicated bin set that consists of both merged bin set and all long complete contigs"}
		depth: {help: "Depth file from BAM files"}
		# The quality filters that are applied to MAGs
		min_completeness: {help: "The minimum percent completeness for a genome bin; default value is set to 70"}
		max_contamination: {help: "The maximum percent contamination for a genome bin; default value is set to 10"}
		max_contigs: {help: "The maximum number of contigs allowed in a genome bin; default value is set to 20"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task checkm2_bin_analysis {
	input {
		String sample
		
		File db
		File dbrep_bins

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int disk_size = ceil(size(db, "GB") * 2 + size(dbrep_bins, "GB") + 20)

	command <<<
		set -euo pipefail

		checkm2 predict \
			-i ~{dbrep_bins} \
			-o "quality_report.tsv" \
			-x fa \
			-t ~{threads} \
			--force
			--remove_intermediates \
			--database_path ~{db} \
			--tmpdir {params.tmp} #TODO
	>>>

	output {
		File qv = "quality_report.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/checkm2@" #TODO
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

task assess_checkm2_bins {
	input {
		String sample
		
		File qv
		File depth

		Int min_completeness
		Int max_contamination
		Int max_contigs

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(qv, "GB") * 2 + size(depth, "GB") + 20)

	command <<<
		set -euo pipefail

		Filter-Checkm2-Bins.py \
			-i ~{qv} \
			-b {params.bin_dir} \ #TODO
			-d ~{depth} \\
			-c1 ~{min_completeness} \
			-c2 ~{max_contamination} \
			-c3 ~{max_contigs} \
			-o "~{sample}.GTDBTk_batch_file.txt" \
			-t "~{sample}.BinCount.txt" \
			-u "~{sample}.quality_report.tsv"
	>>>

	output {
		File gtdb = "~{sample}.GTDBTk_batch_file.txt"
		File target = "~{sample}.BinCount.txt"
		File tsv = "~{sample}.quality_report.tsv"
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
