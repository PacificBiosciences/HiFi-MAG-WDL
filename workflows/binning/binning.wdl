version 1.0

# Binning of incomplete contigs (contigs <500kb)

import "../wdl-common/wdl/structs.wdl"

workflow binning {
	input {
		String sample_id

		File incomplete_contigs_fasta
		File filtered_contig_depth_txt
		File sorted_bam

		Int metabat2_min_contig_size
		String semibin2_model_flag
		String dastool_search_engine
		Float dastool_score_threshold

		RuntimeAttributes default_runtime_attributes
	}

	call metabat2_analysis {
		input:
			sample_id = sample_id,
			incomplete_contigs_fasta = incomplete_contigs_fasta,
			filtered_contig_depth_txt = filtered_contig_depth_txt,
			metabat2_min_contig_size = metabat2_min_contig_size,
			runtime_attributes = default_runtime_attributes
	}

	call semibin2_analysis {
		input:
			incomplete_contigs_fasta = incomplete_contigs_fasta,
			sorted_bam = sorted_bam,
			semibin2_model_flag = semibin2_model_flag,
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
			dastool_search_engine = dastool_search_engine,
			dastool_score_threshold = dastool_score_threshold,
			runtime_attributes = default_runtime_attributes
	}

	output {
		# MetaBAT2
		Array[File] metabat2_reconstructed_bins_fastas = metabat2_analysis.discovered_bins_fastas
		File metabat2_bin_sets_tsv = dastool_input_semibin2.bin_sets_tsv

		# SemiBin2
		File semibin2_bins_tsv = semibin2_analysis.bins_tsv
		Array[File] semibin2_reconstructed_bins_fastas = semibin2_analysis.reconstructed_bins_fastas
		File semibin2_bin_sets_tsv = dastool_input_semibin2.bin_sets_tsv

		# DAS Tool
		Array[File] dastool_bins = dastool_analysis.dastool_bins
	}
}

task metabat2_analysis {
	input {
		String sample_id
		
		File incomplete_contigs_fasta
		File filtered_contig_depth_txt

		Int metabat2_min_contig_size

		RuntimeAttributes runtime_attributes
	}

	Int threads = 12
	Int mem_gb = threads * 4
	Int disk_size = ceil(size(incomplete_contigs_fasta, "GB") * 2 + size(filtered_contig_depth_txt, "GB") + 20)

	command <<<
		set -euo pipefail

		# TODO - get metabat version. It's in the --help message

		metabat2 \
			--verbose \
			--inFile ~{incomplete_contigs_fasta} \
			--abdFile ~{filtered_contig_depth_txt} \
			--outFile ~{sample_id} \
			--numThreads ~{threads} \
			--minContig ~{metabat2_min_contig_size}
	>>>

	output {
		Array[File] discovered_bins_fastas = glob("~{sample_id}.*.fa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/metabat:5e8307c"
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

		String semibin2_model_flag

		RuntimeAttributes runtime_attributes
	}
 
 	Int threads = 48
 	Int mem_gb = threads * 4
	Int disk_size = ceil(size(incomplete_contigs_fasta, "GB") * 2 + size(sorted_bam, "GB") + 20)

	command <<<
		set -euo pipefail

		SemiBin --version

		mkdir semibin2_out_dir

		SemiBin \
			single_easy_bin \
			--input-fasta ~{incomplete_contigs_fasta} \
			--input-bam ~{sorted_bam} \
			--output semibin2_out_dir \
			--self-supervised \
			--sequencing-type=long_reads \
			--compression=none \
			--threads ~{threads} \
			--tag-output \
			semibin2 ~{semibin2_model_flag} \
			--verbose
	>>>

	output {
		File bins_tsv = "semibin2_out_dir/bins_info.tsv"
		Array[File] reconstructed_bins_fastas = glob("semibin2_out_dir/output_bins/semibin2_*.fa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/semibin:5e8307c"
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

		DAS_Tool --version

		reconstructed_bins_fastas_dir=$(dirname ~{reconstructed_bins_fastas[0]})

		Fasta_to_Contig2Bin.sh \
			--input_folder "${reconstructed_bins_fastas_dir}" \
			--extension fa 1> "~{sample_id}.~{binning_algorithm}.tsv"
	>>>

	output {
		File bin_sets_tsv = "~{sample_id}.~{binning_algorithm}.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/dastool:5e8307c"
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

		String dastool_search_engine
		Float dastool_score_threshold

		RuntimeAttributes runtime_attributes
	}
 	
 	Int threads = 24
 	Int mem_gb = threads * 4
	Int disk_size = ceil(size(incomplete_contigs_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		DAS_Tool --version

		DAS_Tool \
			--bins ~{metabat2_bin_sets_tsv},~{semibin2_bin_sets_tsv} \
			--contigs ~{incomplete_contigs_fasta} \
			--labels metabat2,semibin2 \
			--outputbasename ~{sample_id} \
			--search_engine ~{dastool_search_engine} \
			--write_bins \
			--threads ~{threads} \
			--score_threshold ~{dastool_score_threshold} \
			--debug
	>>>

	output {
		Array[File] dastool_bins = glob("~{sample_id}_DASTool_bins/*.fa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/dastool:5e8307c"
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
