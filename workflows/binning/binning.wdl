version 1.0

# Binning of incomplete contigs (contigs <500kb)

import "../wdl-common/wdl/structs.wdl"

workflow binning {
	input {
		String sample_id

		File incomplete_contigs_fasta
		File filtered_depth
		File sorted_bam

		Int min_contig_size = 30000
		String model_flag = 'model: "--environment=global"'

		String search_engine = "diamond"
		Float score_threshold = 0.2

		RuntimeAttributes default_runtime_attributes
	}

	call metabat2_analysis {
		input:
			sample_id = sample_id,
			incomplete_contigs_fasta = incomplete_contigs_fasta,
			filtered_depth = filtered_depth,
			min_contig_size = min_contig_size,
			runtime_attributes = default_runtime_attributes
	}

	call semibin2_analysis {
		input:
			incomplete_contigs_fasta = incomplete_contigs_fasta,
			sorted_bam = sorted_bam,
			model_flag = model_flag,
			runtime_attributes = default_runtime_attributes
	}

	call dastool_input as dastool_input_metabat2 {
		input:
			sample_id = sample_id,
			binning_algorithm = "metabat2",
			reconstructed_bins_fastas = metabat2_analysis.discovered_bins_fastas,
			runtime_attributes = default_runtime_attributes
	}

	call dastool_input as dastool_input_semibin2 {
		input:
			sample_id = sample_id,
			binning_algorithm = "semibin2",
			reconstructed_bins_fastas = semibin2_analysis.reconstructed_bins_fastas,
			runtime_attributes = default_runtime_attributes
	}

	call dastool_analysis {
		input:
			sample_id = sample_id,
			incomplete_contigs_fasta = incomplete_contigs_fasta,
			metabat2_bin_sets_tsv = dastool_input_metabat2.bin_sets_tsv,
			semibin2_bin_sets_tsv = dastool_input_semibin2.bin_sets_tsv,
			search_engine = search_engine,
			score_threshold = score_threshold,
			runtime_attributes = default_runtime_attributes
	}

	output {
		# Metabat2
		Array[File] metabat2_reconstructed_bins_fastas = metabat2_analysis.discovered_bins_fastas
		File metabat2_bin_sets_tsv = dastool_input_semibin2.bin_sets_tsv

		# Semibin2
		File semibin2_bins_tsv = semibin2_analysis.bins_tsv
		Array[File] semibin2_reconstructed_bins_fastas = semibin2_analysis.reconstructed_bins_fastas
		File semibin2_bin_sets_tsv = dastool_input_semibin2.bin_sets_tsv

		# DAS Tool
		Array[File] dastool_bins = dastool_analysis.dastool_bins
	}

	parameter_meta {
		sample_id: {help: "Sample ID"}
		incomplete_contigs_fasta: {help: "Incomplete contigs identified (<93% complete or <500kb) in fasta format"}
		filtered_depth: {help: "A tab-delimited table of the average coverage depth and variance for the incomplete contigs"}
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
		String sample_id
		
		File incomplete_contigs_fasta
		File filtered_depth

		Int min_contig_size

		RuntimeAttributes runtime_attributes
	}

	Int threads = 12
	Int mem_gb = threads * 4
	Int disk_size = ceil(size(incomplete_contigs_fasta, "GB") * 2 + size(filtered_depth, "GB") + 20)

	command <<<
		set -euo pipefail

		metabat2 \
			-v \
			-i ~{incomplete_contigs_fasta} \
			-a ~{filtered_depth} \
			-o ~{sample_id} \
			-t ~{threads} \
			-m ~{min_contig_size}
	>>>

	output {
		Array[File] discovered_bins_fastas = glob("~{sample_id}.*.fa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/metabat@sha256:d68d401803c4e99d85a4c93e48eb223275171f08b49d3314ff245a2d68264651"
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

task semibin2_analysis {
	input {
		File incomplete_contigs_fasta
		File sorted_bam

		String model_flag

		RuntimeAttributes runtime_attributes
	}
 
 	Int threads = 48
 	Int mem_gb = threads * 4
	Int disk_size = ceil(size(incomplete_contigs_fasta, "GB") * 2 + size(sorted_bam, "GB") + 20)

	command <<<
		set -euo pipefail

		SemiBin \
			single_easy_bin \
			-i ~{incomplete_contigs_fasta} \
			-b ~{sorted_bam} \
			-o ./ \
			--self-supervised \
			--sequencing-type=long_reads \
			--compression=none \
			-t ~{threads} \
			--tag-output \
			semibin2 ~{model_flag} \
			--verbose
	>>>

	output {
		File bins_tsv = "bins_info.tsv"
		Array[File] reconstructed_bins_fastas = glob("output_bins/semibin2_*.fa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/semibin@sha256:ddea641ec796eb6259aa6e1298c64e49f70e49b2057c2d00a626aa525c9d8cb4"
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

task dastool_input {
	input {
		String sample_id
		String binning_algorithm

		Array[File] reconstructed_bins_fastas

		RuntimeAttributes runtime_attributes
	}
 	
	Int disk_size = ceil(size(reconstructed_bins_fastas, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		reconstructed_bins_fastas_dir=$(dirname ~{reconstructed_bins_fastas[0]})

		Fasta_to_Contig2Bin.sh \
			-i "${reconstructed_bins_fastas_dir}" \
			-e fa 1> "~{sample_id}.~{binning_algorithm}.tsv"
	>>>

	output {
		File bin_sets_tsv = "~{sample_id}.~{binning_algorithm}.tsv"
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
		String sample_id

		File incomplete_contigs_fasta
		File metabat2_bin_sets_tsv
		File semibin2_bin_sets_tsv

		String search_engine
		Float score_threshold

		RuntimeAttributes runtime_attributes
	}
 	
 	Int threads = 24
 	Int mem_gb = threads * 4
	Int disk_size = ceil(size(incomplete_contigs_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		DAS_Tool \
			-i ~{metabat2_bin_sets_tsv},~{semibin2_bin_sets_tsv} \
			-c ~{incomplete_contigs_fasta} \
			-l metabat2,semibin2 \
			-o ~{sample_id} \
			--search_engine ~{search_engine} \
			--write_bins \
			-t ~{threads} \
			--score_threshold ~{score_threshold} \
			--debug
	>>>

	output {
		Array[File] dastool_bins = glob("~{sample_id}_DASTool_bins/*.fa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/dastool@sha256:ff71f31b7603be06c738e12ed1a91d1e448d0b96a53e073188b1ef870dd8da1c"
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
