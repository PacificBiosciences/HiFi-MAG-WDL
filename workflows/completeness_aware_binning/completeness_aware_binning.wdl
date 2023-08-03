version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow completeness_aware_binning {
	input {
		String sample_id
		File contigs_fasta
		File checkm2_ref_db

		Int min_contig_length
		Int min_contig_completeness

		RuntimeAttributes default_runtime_attributes
	}

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

	output {
		File bins_contigs_key_txt = long_contigs_to_bins.bins_contigs_key_txt
		Array[File] long_bin_fastas = long_contigs_to_bins.long_bin_fastas
		File incomplete_contigs_fasta = make_incomplete_contigs.incomplete_contigs_fasta

		File? contig_quality_report_tsv = checkm2_contig_analysis.contig_quality_report_tsv
		File? passed_bins_txt = filter_complete_contigs.passed_bins_txt
		File? scatterplot_pdf = filter_complete_contigs.scatterplot_pdf
		File? histogram_pdf = filter_complete_contigs.histogram_pdf
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
			mv "${original_fasta}" "${renamed_fasta}"
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
			--input "$long_bin_fastas_dir" \
			--output-directory checkm2_out_dir \
			--extension fa \
			--threads ~{threads} \
			--database_path ~{checkm2_ref_db} \
			--remove_intermediates \
			--tmpdir tmp_dir
	>>>

	output {
		File contig_quality_report_tsv = "checkm2_out_dir/quality_report.tsv"
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
