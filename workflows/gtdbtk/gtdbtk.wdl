version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow gtdbtk {
	input {
		String sample

		File gtdb
		File target
		File gtdbtk_data

		RuntimeAttributes default_runtime_attributes
	}

	# Checkpoint 2: Fork 2 - Bins passed filters; GTDBTkAnalysis -> GTDBTkCleanup -> MAGSummary -> MAGCopy -> MAGPlots
	call gtdbtk_analysis {
		input:
			sample = sample,
			gtdb = gtdb,
			gtdbtk_data = gtdbtk_data,
			target = target, #TODO - remove?
			runtime_attributes = default_runtime_attributes
	}

	## Grab all summary files from GTDB-Tk and merge
	call gtdbtk_cleanup {
		input:
			sample = sample,
			dir_classify = dir_classify,
			complete = gtdbtk_analysis.complete,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File 
	}

	parameter_meta {
		sample: {help: "Sample name"}
		gtdb: {help: "File describing genomes - tab separated in 2 or 3 columns (FASTA file, genome ID, translation table [optional])"
		gtdbtk_data: {help: "The path to the GTDB-Tk database (Genome Database Taxonomy Toolkit)"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task gtdbtk_analysis {
	input {
		String sample
		
		File gtdb
		File target
		File gtdbtk_data

		RuntimeAttributes runtime_attributes
	}

	Int threads = 48
	Int disk_size = ceil(size(gtdb, "GB") * 2 + size(target, "GB") + size(gtdbtk_data, "GB") + 20)

	command <<<
		set -euo pipefail

		GTDBTK_DATA_PATH=~{gtdbtk_data}

		gtdbtk classify_wf \
			--batchfile ~{gtdb} \
			--out_dir ./ \ #TODO
			-x fa \
			--prefix ~{sample} \
			--cpus ~{threads}
	>>>

	output {
		# TODO
		Array[File] align
		Array[File] classify
		Array[File] identify
		File complete = "~{sample}.Complete.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/gtdbtk@sha256:5c92b776841612e07f1c2375a3e950effcf47c0c79d4e2377ade06cafa717fa7"
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

task gtdbtk_cleanup {
	input {
		String sample
		
		File complete

		RuntimeAttributes runtime_attributes
	}

	# TODO - combine with task above to avoid reading in ~{classify} twice?
	String dir_classify = dirname(classify) 
	Int disk_size = ceil(size(complete, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		GTDBTk-Organize.py \
			-i ~{dir_classify} \
			-o "~{sample}.GTDBTk_Summary.txt"
	>>>

	output {
		File gtdbk = "~{sample}.GTDBTk_Summary.txt"
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
