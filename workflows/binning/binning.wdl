version 1.0

# Binning of incomplete contigs (contigs <500kb)

import "../wdl-common/wdl/structs.wdl"

workflow binning {
	input {
		String sample

		File incomplete_contigs
		File filtered_depth

		Int min_contig_size = 30000
		String model_flag = 'model: "--environment=global"'

		String search_engine = "diamond"
		Int score_threshold = 0.2

		RuntimeAttributes default_runtime_attributes
	}

	call metabat2_analysis {
		input:
			sample = sample,
			incomplete_contigs = incomplete_contigs,
			filtered_depth = filtered_depth,
			min_contig_size = min_contig_size,
			runtime_attributes = default_runtime_attributes
	}

	call semibin2_analysis {
		input:
			sample = sample,
			incomplete_contigs = incomplete_contigs,
			bam = bam,
			model_flag = model_flag,
			runtime_attributes = default_runtime_attributes
	}

	call dastool_input as dastool_input_metabat2 {
		input:
			sample = sample,
			binning_algorithm = "metabat2",
			runtime_attributes = default_runtime_attributes
	}

	call dastool_input as dastool_input_semibin2 {
		input:
			sample = sample,
			binning_algorithm = "semibin2",
			runtime_attributes = default_runtime_attributes
	}

	call dastool_analysis {
		input:
			sample = sample,
			incomplete_contigs = incomplete_contigs,
			metabat = dastool_input_metabat2.bin_sets,
			semibin = dastool_input_semibin2.bin_sets,
			search_engine = search_engine,
			score_threshold = score_threshold,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File 
	}

	parameter_meta {
		sample: {help: "Sample name"}
		incomplete_contigs: {help: ""} #TODO
		filtered_depth: {help: "Filtered depth file"} #TODO
		# Metabat2 params
		min_contig_size: {help: "The minimum size of contig to be included in binning for Metabat2; default value is set to 30000"}
		# Semibin2 params
		model_flag: {help: "The trained model to be used in Semibin2; default value is set to 'global'"}
		# DAS tool params
		search_engine: {help: "The engine for single copy gene searching; default is set to 'diamond'"}
		score_threshold: {help: "Score threshold until selection algorithm will keep selecting bins [0 to 1]; default value is set to 0.2 (20%)"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task metabat2_analysis {
	input {
		String sample
		
		File incomplete_contigs
		File filtered_depth

		Int min_contig_size

		RuntimeAttributes runtime_attributes
	}

	Int threads = 12
	Int disk_size = ceil(size(incomplete_contigs, "GB") * 2 + size(filtered_depth, "GB") + 20)

	command <<<
		set -euo pipefail

		metabat2 \
			-v \
			-i ~{incomplete_contigs} \
			-a ~{filtered_depth} \
			-o ~{sample} \
			-t ~{threads} \
			-m ~{min_contig_size}
	>>>

	output {
		File read_1 = "~{sample}.1.fa"
		File read_2 = "~{sample}.2.fa"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/metabat@sha256:d68d401803c4e99d85a4c93e48eb223275171f08b49d3314ff245a2d68264651"
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

task semibin2_analysis {
	input {
		String sample
		
		File incomplete_contigs
		File bam

		String model_flag

		RuntimeAttributes runtime_attributes
	}
 
 	Int threads = 48
	Int disk_size = ceil(size(incomplete_contigs, "GB") * 2 + size(bam, "GB") + 20)

	command <<<
		set -euo pipefail

		SemiBin \
			single_easy_bin \
			-i ~{incomplete_contigs} \
			-b ~{bam} \
			-o {output.outdir} \
			--self-supervised \
			--sequencing-type=long_reads \
			--compression=none \
			-t ~{threads} \
			--tag-output \
			semibin2 ~{model_flag} \
			--verbose \
			--tmpdir={params.tmp} #TODO
	>>>

	output {
		File bins = "bins_info.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/semibin@sha256:ddea641ec796eb6259aa6e1298c64e49f70e49b2057c2d00a626aa525c9d8cb4"
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

task dastool_input {
	input {
		String sample
		String binning_algorithm

		# TODO
		File 

		RuntimeAttributes runtime_attributes
	}
 
	Int disk_size = ceil(size(depth, "GB") * 2 + size(passed, "GB") + 20)

	command <<<
		set -euo pipefail

		Fasta_to_Contig2Bin.sh \
			-i {params.indir} \ #TODO
			-e fa 1> "~{sample}.~{binning_algorithm}.tsv"
	>>>

	output {
		File bin_sets = "~{sample}.~{binning_algorithm}.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/dastool@sha256:ff71f31b7603be06c738e12ed1a91d1e448d0b96a53e073188b1ef870dd8da1c"
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

task dastool_analysis {
	input {
		String sample

		File incomplete_contigs
		File metabat
		File semibin

		String search_engine
		Int score_threshold

		RuntimeAttributes runtime_attributes
	}
 	
 	Int threads = 24
	Int disk_size = ceil(size(depth, "GB") * 2 + size(passed, "GB") + 20)

	command <<<
		set -euo pipefail

		DAS_Tool \
			-i ~{metabat},~{semibin} \
			-c ~{incomplete_contigs} \
			-l metabat2,semibin2 \
			-o {params.outlabel} \ #TODO
			--search_engine ~{search_engine} \
			--write_bins \
			-t ~{threads} \
			--score_threshold ~{score_threshold} \
			--debug
	>>>

	output {
		File complete = "~{sample}.Complete.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/dastool@sha256:ff71f31b7603be06c738e12ed1a91d1e448d0b96a53e073188b1ef870dd8da1c"
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
