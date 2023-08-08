version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow bin_reads {
	input {
		String sample_id
		File contigs_fasta
		File contigs_fasta_gz
		File hifi_reads_fastq

		File checkm2_ref_db
		Int min_contig_length
		Int min_contig_completeness

		Int metabat2_min_contig_size
		String semibin2_model
		String dastool_search_engine
		Float dastool_score_threshold

		Int min_mag_completeness
		Int max_mag_contamination
		Int max_contigs

		RuntimeAttributes default_runtime_attributes
	}

	# Completeness-aware binning
	call long_contigs_to_bins {
		input:
			sample_id = sample_id,
			contigs_fasta = contigs_fasta,
			min_contig_length = min_contig_length,
			runtime_attributes = default_runtime_attributes
	}

	call make_incomplete_contigs {
		input:
			sample_id = sample_id,
			contigs_fasta = contigs_fasta,
			bins_contigs_key_txt = long_contigs_to_bins.bins_contigs_key_txt,
			renamed_long_bin_fastas = long_contigs_to_bins.renamed_long_bin_fastas,
			runtime_attributes = default_runtime_attributes
	}

	if (long_contigs_to_bins.bin_key_nonempty) {
		call checkm2_contig_analysis {
		input:
			sample_id = sample_id,
			checkm2_ref_db = checkm2_ref_db,
			long_bin_fastas = long_contigs_to_bins.long_bin_fastas,
			runtime_attributes = default_runtime_attributes
		}

		call filter_complete_contigs {
		input:
			sample_id = sample_id,
			contigs_fasta = contigs_fasta,
			contig_quality_report_tsv = checkm2_contig_analysis.contig_quality_report_tsv,
			bins_contigs_key_txt = long_contigs_to_bins.bins_contigs_key_txt,
			min_contig_length = min_contig_length,
			min_contig_completeness = min_contig_completeness,
			runtime_attributes = default_runtime_attributes
		}
	}

	# Calculate coverage
	call align_hifiasm {
		input:
			sample_id = sample_id,
			contigs_fasta_gz = contigs_fasta_gz,
			hifi_reads_fastq = hifi_reads_fastq,
			runtime_attributes = default_runtime_attributes
	}

	call jgi_bam_depth {
		input:
			sample_id = sample_id,
			aligned_sorted_bam = align_hifiasm.aligned_sorted_bam,
			aligned_sorted_bam_index = align_hifiasm.aligned_sorted_bam_index,
			runtime_attributes = default_runtime_attributes
	}

	call convert_jgi_bamdepth {
		input:
			sample_id = sample_id,
			contig_depth_txt = jgi_bam_depth.contig_depth_txt,
			bins_contigs_key_txt = long_contigs_to_bins.bins_contigs_key_txt,
			runtime_attributes = default_runtime_attributes
	}

	# Bin incomplete contigs (contigs <500kb)
	call metabat2_analysis {
		input:
			sample_id = sample_id,
			incomplete_contigs_fasta = make_incomplete_contigs.incomplete_contigs_fasta,
			filtered_contig_depth_txt = convert_jgi_bamdepth.filtered_contig_depth_txt,
			metabat2_min_contig_size = metabat2_min_contig_size,
			runtime_attributes = default_runtime_attributes
	}

	call dastool_input as dastool_input_metabat2 {
		input:
			sample_id = sample_id,
			binning_algorithm = "metabat2",
			reconstructed_bins_fastas = metabat2_analysis.discovered_bins_fastas,
			runtime_attributes = default_runtime_attributes
	}

	call semibin2_analysis {
		input:
			sample_id = sample_id,
			incomplete_contigs_fasta = make_incomplete_contigs.incomplete_contigs_fasta,
			aligned_sorted_bam = align_hifiasm.aligned_sorted_bam,
			semibin2_model = semibin2_model,
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
			incomplete_contigs_fasta = make_incomplete_contigs.incomplete_contigs_fasta,
			metabat2_bin_sets_tsv = dastool_input_metabat2.bin_sets_tsv,
			semibin2_bin_sets_tsv = dastool_input_semibin2.bin_sets_tsv,
			dastool_search_engine = dastool_search_engine,
			dastool_score_threshold = dastool_score_threshold,
			runtime_attributes = default_runtime_attributes
	}

	# Assess bin quality
	call checkm2_bin_analysis {
		input:
			sample_id = sample_id,
			checkm2_ref_db = checkm2_ref_db,
			long_bin_fastas = long_contigs_to_bins.long_bin_fastas,
			dastool_bins = dastool_analysis.dastool_bins,
			runtime_attributes = default_runtime_attributes
	}

	call assess_checkm2_bins {
		input:
			sample_id = sample_id,
			bin_quality_report_tsv = checkm2_bin_analysis.bin_quality_report_tsv,
			filtered_contig_depth_txt = convert_jgi_bamdepth.filtered_contig_depth_txt,
			derep_bins = checkm2_bin_analysis.derep_bins,
			min_mag_completeness = min_mag_completeness,
			max_mag_contamination = max_mag_contamination,
			max_contigs = max_contigs,
			runtime_attributes = default_runtime_attributes
	}

	output {
		# Completeness-aware binning
		File bins_contigs_key_txt = long_contigs_to_bins.bins_contigs_key_txt
		Array[File] long_bin_fastas = long_contigs_to_bins.long_bin_fastas
		File incomplete_contigs_fasta = make_incomplete_contigs.incomplete_contigs_fasta

		File? contig_quality_report_tsv = checkm2_contig_analysis.contig_quality_report_tsv
		File? passed_bins_txt = filter_complete_contigs.passed_bins_txt
		File? scatterplot_pdf = filter_complete_contigs.scatterplot_pdf
		File? histogram_pdf = filter_complete_contigs.histogram_pdf

		# Coverage
		IndexData aligned_sorted_bam = {
			"data": align_hifiasm.aligned_sorted_bam,
			"data_index": align_hifiasm.aligned_sorted_bam_index
		}
		File filtered_contig_depth_txt = convert_jgi_bamdepth.filtered_contig_depth_txt

		# Incomplete contig binning
		Array[File] metabat2_reconstructed_bins_fastas = metabat2_analysis.discovered_bins_fastas
		File metabat2_bin_sets_tsv = dastool_input_semibin2.bin_sets_tsv
		File semibin2_bins_tsv = semibin2_analysis.bins_tsv
		Array[File] semibin2_reconstructed_bins_fastas = semibin2_analysis.reconstructed_bins_fastas
		File semibin2_bin_sets_tsv = dastool_input_semibin2.bin_sets_tsv
		Array[File] dastool_bins = dastool_analysis.dastool_bins

		# Bin quality
		Array[File] derep_bins = checkm2_bin_analysis.derep_bins
		File bin_quality_report_tsv = checkm2_bin_analysis.bin_quality_report_tsv
		File gtdb_batch_txt = assess_checkm2_bins.gtdb_batch_txt
		File passed_bin_count_txt = assess_checkm2_bins.passed_bin_count_txt
		File filtered_quality_report_tsv = assess_checkm2_bins.filtered_quality_report_tsv
		Boolean bin_count_nonempty = assess_checkm2_bins.bin_count_nonempty
	}
}

task long_contigs_to_bins {
	input {
		String sample_id
		File contigs_fasta

		Int min_contig_length

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contigs_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		mkdir long_bin_fastas_out_dir

		python /opt/scripts/Fasta-Make-Long-Seq-Bins.py \
			--input_fasta ~{contigs_fasta} \
			--bins_contigs "~{sample_id}.bin_key.txt" \
			--length ~{min_contig_length} \
			--outdir long_bin_fastas_out_dir

		# Check if any long contigs (>500kb) were identified
		if [[ -s "~{sample_id}.bin_key.txt" ]]; then
			echo "true" > bin_key_nonempty.txt
		else
			echo "false" > bin_key_nonempty.txt
		fi

		# Rename long bin fastas from "complete.x.fa" to "${contig}.fa" in order to use Make-Incomplete-Contigs.py script in next task
		mkdir long_bin_fastas_renamed_out_dir

		while IFS= read -r bins_contigs_line || [[ -n "${bins_contigs_line}" ]]; do
			original_fasta="long_bin_fastas_out_dir/$(echo "${bins_contigs_line}" | cut -f 1).fa"
			renamed_fasta="long_bin_fastas_renamed_out_dir/$(echo "${bins_contigs_line}" | cut -f 2).fa"
			cp "${original_fasta}" "${renamed_fasta}"
		done < ~{sample_id}.bin_key.txt
	>>>

	output {
		File bins_contigs_key_txt = "~{sample_id}.bin_key.txt"
		Boolean bin_key_nonempty = read_boolean("bin_key_nonempty.txt")
		Array[File] long_bin_fastas = glob("long_bin_fastas_out_dir/*.fa")
		Array[File] renamed_long_bin_fastas = glob("long_bin_fastas_renamed_out_dir/*.fa")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python:5e8307c"
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

		File bins_contigs_key_txt
		Array[File] renamed_long_bin_fastas

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(contigs_fasta, "GB") + (size(renamed_long_bin_fastas[0], "GB") * length(renamed_long_bin_fastas))) * 2 + 20)

	command <<<
		set -euo pipefail

		renamed_long_bin_fastas_dir=$(dirname ~{renamed_long_bin_fastas[0]})

		mkdir long_bin_fastas_copy_out_dir

		python /opt/scripts/Make-Incomplete-Contigs.py \
			--input_fasta ~{contigs_fasta} \
			--output_fasta "~{sample_id}.incomplete_contigs.fasta" \
			--passed_bins ~{bins_contigs_key_txt} \
			--fastadir "${renamed_long_bin_fastas_dir}" \
			--outdir long_bin_fastas_copy_out_dir
	>>>

	output {
		File incomplete_contigs_fasta = "~{sample_id}.incomplete_contigs.fasta"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python:5e8307c"
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
		String sample_id
		File checkm2_ref_db

		Array[File] long_bin_fastas

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = threads * 4
	Int disk_size = ceil((size(checkm2_ref_db, "GB") + (size(long_bin_fastas[0], "GB") * length(long_bin_fastas))) * 2 + 20)

	command <<<
		set -euo pipefail

		checkm2 --version

		long_bin_fastas_dir=$(dirname ~{long_bin_fastas[0]})

		mkdir checkm2_out_dir
		mkdir tmp_dir

		checkm2 predict \
			--input "${long_bin_fastas_dir}" \
			--output-directory checkm2_out_dir \
			--extension fa \
			--threads ~{threads} \
			--database_path ~{checkm2_ref_db} \
			--remove_intermediates \
			--tmpdir tmp_dir

		mv checkm2_out_dir/quality_report.tsv checkm2_out_dir/"~{sample_id}.contig.quality_report.tsv"
	>>>

	output {
		File contig_quality_report_tsv = "checkm2_out_dir/~{sample_id}.contig.quality_report.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/checkm2:5e8307c"
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

		File contig_quality_report_tsv
		File bins_contigs_key_txt

		Int min_contig_length
		Int min_contig_completeness

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(contigs_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Filter-Complete-Contigs.py \
			--input_fasta ~{contigs_fasta} \
			--checkm ~{contig_quality_report_tsv} \
			--bins_contigs ~{bins_contigs_key_txt} \
			--length ~{min_contig_length} \
			--min_completeness ~{min_contig_completeness} \
			--passed_bins "~{sample_id}.passed_bins.txt" \
			--plot_scatter "~{sample_id}.completeness_vs_size_scatter.pdf" \
			--plot_histo "~{sample_id}.completeness_histo.pdf"
	>>>

	output {
		File passed_bins_txt = "~{sample_id}.passed_bins.txt"
		File scatterplot_pdf = "~{sample_id}.completeness_vs_size_scatter.pdf"
		File histogram_pdf = "~{sample_id}.completeness_histo.pdf"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python:5e8307c"
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

task align_hifiasm {
	input {
		String sample_id

		File contigs_fasta_gz
		File hifi_reads_fastq

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = threads * 2
	Int disk_size = ceil((size(contigs_fasta_gz, "GB") + size(hifi_reads_fastq, "GB")) * 2 + 20)

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
			~{contigs_fasta_gz} \
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

task jgi_bam_depth {
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

task convert_jgi_bamdepth {
	input {
		String sample_id

		File contig_depth_txt
		File bins_contigs_key_txt

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(contig_depth_txt, "GB") + size(bins_contigs_key_txt, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/Convert-JGI-Coverages.py \
			--in_jgi ~{contig_depth_txt} \
			--passed_bins ~{bins_contigs_key_txt} \
			--out_jgi "~{sample_id}.JGI.filtered.depth.txt"
	>>>

	output {
		File filtered_contig_depth_txt = "~{sample_id}.JGI.filtered.depth.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python:5e8307c"
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
	Int disk_size = ceil(size(incomplete_contigs_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		metabat --help |& grep "version"

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
		String sample_id
		File incomplete_contigs_fasta
		File aligned_sorted_bam

		String semibin2_model

		RuntimeAttributes runtime_attributes
	}

	Int threads = 48
	Int mem_gb = threads * 4
	Int disk_size = ceil((size(incomplete_contigs_fasta, "GB") + size(aligned_sorted_bam, "GB")) * 2 + 20)

	# If the model is set to the empty string, a new model will be built
	String semibin2_model_flag = if (semibin2_model == "") then "" else "--environment=~{semibin2_model}"

	command <<<
		set -euo pipefail

		SemiBin --version

		mkdir semibin2_out_dir

		SemiBin \
			single_easy_bin \
			semibin2 \
			--input-fasta ~{incomplete_contigs_fasta} \
			--input-bam ~{aligned_sorted_bam} \
			--output semibin2_out_dir \
			--self-supervised \
			--sequencing-type=long_reads \
			--compression=none \
			--threads ~{threads} \
			--tag-output \
			~{semibin2_model_flag} \
			--verbose

		mv semibin2_out_dir/bins_info.tsv semibin2_out_dir/"~{sample_id}.bins_info.tsv"
	>>>

	output {
		File bins_tsv = "semibin2_out_dir/~{sample_id}.bins_info.tsv"
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

	Int disk_size = ceil(size(reconstructed_bins_fastas[0], "GB") * length(reconstructed_bins_fastas) * 2 + 20)

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

task checkm2_bin_analysis {
	input {
		String sample_id
		File checkm2_ref_db
		Array[File] long_bin_fastas
		Array[File] dastool_bins

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = threads * 4
	Int disk_size = ceil((size(checkm2_ref_db, "GB") + (size(long_bin_fastas[0], "GB") * length(long_bin_fastas)) + (size(dastool_bins[0], "GB") * length(dastool_bins))) * 2 + 20)

	command <<<
		set -euo pipefail

		checkm2 --version

		mkdir derep_bins_dir
		mkdir checkm2_out_dir
		mkdir tmp_dir

		while read -r fasta || [[ -n "${fasta}" ]]; do
			ln -s "${fasta}" "$(pwd)"/derep_bins_dir/
		done < ~{write_lines(long_bin_fastas)}

		while read -r fasta || [[ -n "${fasta}" ]]; do
			ln -s "${fasta}" "$(pwd)"/derep_bins_dir/
		done < ~{write_lines(dastool_bins)}

		checkm2 predict \
			--input derep_bins_dir \
			--output-directory checkm2_out_dir \
			--extension fa \
			--threads ~{threads} \
			--database_path ~{checkm2_ref_db} \
			--remove_intermediates \
			--tmpdir tmp_dir

		mv checkm2_out_dir/quality_report.tsv checkm2_out_dir/"~{sample_id}.bin.quality_report.tsv"
	>>>

	output {
		Array[File] derep_bins = glob("derep_bins_dir/*.fa")
		File bin_quality_report_tsv = "checkm2_out_dir/~{sample_id}.bin.quality_report.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/checkm2:5e8307c"
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

task assess_checkm2_bins {
	input {
		String sample_id

		File bin_quality_report_tsv
		File filtered_contig_depth_txt
		Array[File] derep_bins

		Int min_mag_completeness
		Int max_mag_contamination
		Int max_contigs

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(derep_bins[0], "GB") * length(derep_bins) * 2 + 20)

	command <<<
		set -euo pipefail

		derep_bins_dir=$(dirname ~{derep_bins[0]})

		python /opt/scripts/Filter-Checkm2-Bins.py \
			--input_tsv ~{bin_quality_report_tsv} \
			--bin_dir "${derep_bins_dir}" \
			--depth_file ~{filtered_contig_depth_txt} \
			--min_completeness ~{min_mag_completeness} \
			--max_contamination ~{max_mag_contamination} \
			--max_contigs ~{max_contigs} \
			--gtdb_outfile "~{sample_id}.GTDBTk_batch_file.txt" \
			--target_outfile "~{sample_id}.BinCount.txt" \
			--updated_tsv "~{sample_id}.quality_report.tsv"

		# Check if there are bins after CheckM2, before running GTDB-Tk and the MAG summary
		if [[ -s "~{sample_id}.BinCount.txt" ]]; then
			echo "true" > bin_count_nonempty.txt
		else
			echo "false" > bin_count_nonempty.txt
			echo "No bins passed filtering in CheckM2, see ~{sample_id}.quality_report.tsv for more information"
		fi
	>>>

	output {
		File gtdb_batch_txt = "~{sample_id}.GTDBTk_batch_file.txt"
		File passed_bin_count_txt = "~{sample_id}.BinCount.txt"
		File filtered_quality_report_tsv = "~{sample_id}.quality_report.tsv"
		Boolean bin_count_nonempty = read_boolean("bin_count_nonempty.txt")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/python:5e8307c"
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
