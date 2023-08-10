version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow bin_incomplete_contigs {
	input {
		String sample_id
		File assembled_contigs_fa_gz
		File hifi_reads_fastq

		File incomplete_contigs_fa
		File passing_long_contig_bin_map

		Int metabat2_min_contig_size
		String semibin2_model
		String dastool_search_engine
		Float dastool_score_threshold

		RuntimeAttributes default_runtime_attributes
	}

	# Calculate coverage
	call align_reads_to_assembled_contigs {
		input:
			sample_id = sample_id,
			assembled_contigs = assembled_contigs_fa_gz,
			hifi_reads_fastq = hifi_reads_fastq,
			runtime_attributes = default_runtime_attributes
	}

	# This gets depth for all contigs, not just incomplete ones
	call summarize_contig_depth {
		input:
			sample_id = sample_id,
			aligned_sorted_bam = align_reads_to_assembled_contigs.aligned_sorted_bam,
			aligned_sorted_bam_index = align_reads_to_assembled_contigs.aligned_sorted_bam_index,
			runtime_attributes = default_runtime_attributes
	}

	# Filter the contig depth to only include incomplete contigs
	call filter_contig_depth {
		input:
			contig_depth_txt = summarize_contig_depth.contig_depth_txt,
			long_contig_bin_map = passing_long_contig_bin_map,
			runtime_attributes = default_runtime_attributes
	}

	# Bin incomplete contigs (contigs <500kb)
	call bin_incomplete_contigs_metabat2 {
		input:
			sample_id = sample_id,
			incomplete_contigs_fa = incomplete_contigs_fa,
			filtered_contig_depth_txt = filter_contig_depth.filtered_contig_depth_txt,
			metabat2_min_contig_size = metabat2_min_contig_size,
			runtime_attributes = default_runtime_attributes
	}

	call bin_incomplete_contigs_semibin2 {
		input:
			sample_id = sample_id,
			incomplete_contigs_fa = incomplete_contigs_fa,
			aligned_sorted_bam = align_reads_to_assembled_contigs.aligned_sorted_bam,
			semibin2_model = semibin2_model,
			runtime_attributes = default_runtime_attributes
	}

	call map_contig_to_bin as map_contig_to_bin_metabat2 {
		input:
			sample_id = sample_id,
			binning_algorithm = "metabat2",
			bin_fas = bin_incomplete_contigs_metabat2.bin_fas,
			runtime_attributes = default_runtime_attributes
	}

	call map_contig_to_bin as map_contig_to_bin_semibin2 {
		input:
			sample_id = sample_id,
			binning_algorithm = "semibin2",
			bin_fas = bin_incomplete_contigs_semibin2.bin_fas,
			runtime_attributes = default_runtime_attributes
	}

	call merge_incomplete_bins {
		input:
			sample_id = sample_id,
			incomplete_contigs_fa = incomplete_contigs_fa,
			metabat2_contig_bin_map = map_contig_to_bin_metabat2.contig_bin_map,
			semibin2_contig_bin_map = map_contig_to_bin_semibin2.contig_bin_map,
			dastool_search_engine = dastool_search_engine,
			dastool_score_threshold = dastool_score_threshold,
			runtime_attributes = default_runtime_attributes
	}

	output {
		# Coverage
		IndexData aligned_sorted_bam = {
			"data": align_reads_to_assembled_contigs.aligned_sorted_bam,
			"data_index": align_reads_to_assembled_contigs.aligned_sorted_bam_index
		}

		File contig_depth_txt = summarize_contig_depth.contig_depth_txt

		# Incomplete contig binning
		Array[File] metabat2_bin_fas = bin_incomplete_contigs_metabat2.bin_fas
		File metabat2_contig_bin_map = map_contig_to_bin_metabat2.contig_bin_map

		File semibin2_bins_tsv = bin_incomplete_contigs_semibin2.bins_tsv
		Array[File] semibin2_bin_fas = bin_incomplete_contigs_semibin2.bin_fas
		File semibin2_contig_bin_map = map_contig_to_bin_semibin2.contig_bin_map

		Array[File] merged_incomplete_bin_fas = merge_incomplete_bins.merged_incomplete_bin_fas
	}
}

task align_reads_to_assembled_contigs {
	input {
		String sample_id

		File assembled_contigs
		File hifi_reads_fastq

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = threads * 2
	Int disk_size = ceil((size(assembled_contigs, "GB") + size(hifi_reads_fastq, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		minimap2 --version
		samtools --version

		minimap2 \
			-a \
			-k 19 \
			-w 10 \
			-I 10G \
			-g 5000 \
			-r 2000 \
			--lj-min-ratio 0.5 \
			-A 2 \
			-B 5 \
			-O 5,56 \
			-E 4,1 \
			-z 400,50 \
			--sam-hit-only \
			-t ~{threads / 2} \
			~{assembled_contigs} \
			~{hifi_reads_fastq} \
		| samtools sort \
			-@ ~{threads / 2 - 1} \
			-o "~{sample_id}.aligned.sorted.bam"

		samtools index \
			-@ ~{threads - 1} \
			"~{sample_id}.aligned.sorted.bam"
	>>>

	output {
		File aligned_sorted_bam = "~{sample_id}.aligned.sorted.bam"
		File aligned_sorted_bam_index = "~{sample_id}.aligned.sorted.bam.bai"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools:5e8307c"
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

task summarize_contig_depth {
	input {
		String sample_id

		File aligned_sorted_bam
		File aligned_sorted_bam_index

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(aligned_sorted_bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		metabat --help |& grep "version"

		jgi_summarize_bam_contig_depths \
			--outputDepth "~{sample_id}.JGI.depth.txt" \
			~{aligned_sorted_bam}
	>>>

	output {
		File contig_depth_txt = "~{sample_id}.JGI.depth.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/metabat:5e8307c"
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

# Filter out passing bins from the contig depth JGI
task filter_contig_depth {
	input {
		File contig_depth_txt
		File long_contig_bin_map

		RuntimeAttributes runtime_attributes
	}

	String contig_depth_txt_basename = basename(contig_depth_txt, ".txt")
	Int disk_size = ceil((size(contig_depth_txt, "GB") + size(long_contig_bin_map, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Convert-JGI-Coverages.py \
			--in_jgi ~{contig_depth_txt} \
			--passed_bins ~{long_contig_bin_map} \
			--out_jgi "~{contig_depth_txt_basename}.filtered.txt"
	>>>

	output {
		File filtered_contig_depth_txt = "~{contig_depth_txt_basename}.filtered.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python@sha256:c7e594d86c35d2c3b2cd8fabf51d9274d74347464433c4f3e55e5306be7bd1ea"
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

task bin_incomplete_contigs_metabat2 {
	input {
		String sample_id

		File incomplete_contigs_fa
		File filtered_contig_depth_txt

		Int metabat2_min_contig_size

		RuntimeAttributes runtime_attributes
	}

	Int threads = 12
	Int mem_gb = threads * 4
	Int disk_size = ceil(size(incomplete_contigs_fa, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		metabat --help |& grep "version"

		metabat2 \
			--verbose \
			--inFile ~{incomplete_contigs_fa} \
			--abdFile ~{filtered_contig_depth_txt} \
			--outFile ~{sample_id}.metabat2 \
			--numThreads ~{threads} \
			--minContig ~{metabat2_min_contig_size}
	>>>

	output {
		Array[File] bin_fas = glob("~{sample_id}.metabat2.*.fa")
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

task bin_incomplete_contigs_semibin2 {
	input {
		String sample_id
		File incomplete_contigs_fa
		File aligned_sorted_bam

		String semibin2_model

		RuntimeAttributes runtime_attributes
	}

	Int threads = 48
	Int mem_gb = threads * 4
	Int disk_size = ceil((size(incomplete_contigs_fa, "GB") + size(aligned_sorted_bam, "GB")) * 2 + 20)

	String semibin2_model_flag = if (semibin2_model == "TRAIN") then "" else "--environment ~{semibin2_model}"

	command <<<
		set -euo pipefail

		SemiBin --version

		SemiBin \
			single_easy_bin \
			--input-fasta ~{incomplete_contigs_fa} \
			--input-bam ~{aligned_sorted_bam} \
			--output semibin2_out_dir \
			--self-supervised \
			--sequencing-type long_reads \
			--compression none \
			--threads ~{threads} \
			--tag-output ~{sample_id}.semibin2 \
			~{semibin2_model_flag} \
			--verbose

		mv semibin2_out_dir/bins_info.tsv "semibin2_out_dir/~{sample_id}.bins_info.tsv"
	>>>

	output {
		File bins_tsv = "semibin2_out_dir/~{sample_id}.bins_info.tsv"
		Array[File] bin_fas = glob("semibin2_out_dir/output_bins/~{sample_id}.semibin2_*.fa")
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

task map_contig_to_bin {
	input {
		String sample_id
		String binning_algorithm

		Array[File] bin_fas

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(bin_fas[0], "GB") * length(bin_fas) * 2 + 20)

	command <<<
		set -euo pipefail

		DAS_Tool --version

		bin_fa_dir=$(dirname ~{bin_fas[0]})

		Fasta_to_Contig2Bin.sh \
			--input_folder "${bin_fa_dir}" \
			--extension fa \
		1> "~{sample_id}.~{binning_algorithm}.contig_bin_map.tsv"
	>>>

	output {
		File contig_bin_map = "~{sample_id}.~{binning_algorithm}.contig_bin_map.tsv"
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

task merge_incomplete_bins {
	input {
		String sample_id

		File incomplete_contigs_fa
		File metabat2_contig_bin_map
		File semibin2_contig_bin_map

		String dastool_search_engine
		Float dastool_score_threshold

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = threads * 4
	Int disk_size = ceil(size(incomplete_contigs_fa, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		DAS_Tool --version

		DAS_Tool \
			--contigs ~{incomplete_contigs_fa} \
			--bins ~{metabat2_contig_bin_map},~{semibin2_contig_bin_map} \
			--labels metabat2,semibin2 \
			--outputbasename ~{sample_id} \
			--search_engine ~{dastool_search_engine} \
			--write_bins \
			--threads ~{threads} \
			--score_threshold ~{dastool_score_threshold} \
			--debug
	>>>

	output {
		Array[File] merged_incomplete_bin_fas = glob("~{sample_id}_DASTool_bins/*.fa")
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
