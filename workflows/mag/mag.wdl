version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow mag {
	input {
		String sample

		File gtdbk
		File qv

		Int min_completeness = 70
		Int max_contamination = 10

		RuntimeAttributes default_runtime_attributes
	}

	call mag_summary {
		input:
			sample = sample,
			gtdbk = gtdbk,
			qv = qv,
			runtime_attributes = default_runtime_attributes
	}

	call mag_copy {
		input:
			sample = sample,
			summary = mag_summary.summary,
			runtime_attributes = default_runtime_attributes
	}

	call mag_plots {
		input:
			sample = sample,
			qv = qv,
			summary = mag_summary.summary,
			min_completeness = min_completeness,
			max_contamination = max_contamination,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File summary = mag_summary.summary
		Array[File] plots = mag_plots.plots
	}

	parameter_meta {
		sample: {help: "Sample name"}
		gtdbk: {help: "Summary text for taxnomic classifcations (Genome Database Taxonomy Toolkit)"}
		qv: {help: "The quality report from CheckM2 containing the completeness and contamination information for each bin"}
		# The quality filters that are applied to MAGs
		min_completeness: {help: "The minimum percent completeness for a genome bin; default value is set to 70"}
		max_contamination: {help: "The maximum percent contamination for a genome bin; default value is set to 10"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task mag_summary {
	input {
		String sample
		
		File gtdbk
		File qv

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(gtdbk, "GB") * 2 + size(qv, "GB") + 20)

	command <<<
		set -euo pipefail

		MAG-Summary.py \
			-g ~{gtdbtk} \
			-c ~{qv} \
			-o "~{sample}.HiFi_MAG.summary.txt"
	>>>

	output {
		File summary = "~{sample}.HiFi_MAG.summary.txt"
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

task mag_copy {
	input {
		String sample
		
		File summary

		RuntimeAttributes runtime_attributes
	}
 
	Int disk_size = ceil(size(summary, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		Copy-Final-MAGs.py \
			-i ~{summary} \
			-m {input.mag_dir} \ #TODO - deplicated-bins directory
			-o {output} #TODO - output directory
	>>>

	output {

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

task mag_plots {
	input {
		String sample
		
		File qv
		File summary

		Int min_completeness
		Int max_contamination

		RuntimeAttributes runtime_attributes
	}
 
	Int disk_size = ceil(size(qv, "GB") * 2 + size(summary, "GB") + 20)

	command <<<
		set -euo pipefail

		Plot-Figures.py \
			-i1 ~{qv} \
			-i2 ~{summary} \
			-l ~{sample} \
			-c1 ~{min_completeness} \
			-c2 ~{max_contamination} \
			-o1 "~{sample}.All-DASTool-Bins.pdf" \
			-o2 "~{sample}.Completeness-Contamination-Contigs.pdf" \
			-o3 "~{sample}.GenomeSizes-Depths.pdf"
	>>>

	output {
		Array[File] plots = ["~{sample}.All-DASTool-Bins.pdf", "~{sample}.Completeness-Contamination-Contigs.pdf", "~{sample}.GenomeSizes-Depths.pdf"]
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

