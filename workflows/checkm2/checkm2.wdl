version 1.0

import "../wdl-common/wdl/structs.wdl"

workflow checkm2 {
	input {
		String sample_id

		File checkm2_ref_db
		File filtered_contig_depth_txt
		Array[File] derep_bins

		Int min_mag_completeness
		Int max_mag_contamination
		Int max_contigs

		RuntimeAttributes default_runtime_attributes
	}

	call checkm2_bin_analysis {
		input:
			checkm2_ref_db = checkm2_ref_db,
			derep_bins = derep_bins,
			runtime_attributes = default_runtime_attributes
	}

	call assess_checkm2_bins {
		input:
			sample_id = sample_id,
			bin_quality_report_tsv = checkm2_bin_analysis.bin_quality_report_tsv,
			filtered_contig_depth_txt = filtered_contig_depth_txt,
			derep_bins = derep_bins,
			min_mag_completeness = min_mag_completeness,
			max_mag_contamination = max_mag_contamination,
			max_contigs = max_contigs,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File bin_quality_report_tsv = checkm2_bin_analysis.bin_quality_report_tsv
		File gtdb_batch_txt = assess_checkm2_bins.gtdb_batch_txt
		File passed_bin_count_txt = assess_checkm2_bins.passed_bin_count_txt
		File filtered_quality_report_tsv = assess_checkm2_bins.filtered_quality_report_tsv
		Boolean passed_bin_count_nonempty = assess_checkm2_bins.passed_bin_count_nonempty
	}
}

task checkm2_bin_analysis {
	input {
		File checkm2_ref_db
		Array[File] derep_bins

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = threads * 4
	Int disk_size = ceil(size(checkm2_ref_db, "GB") * 2 + size(derep_bins, "GB") + 20)

	command <<<
		set -euo pipefail

		checkm2 --version

		derep_bins_dir=$(dirname ~{derep_bins[0]})

		mkdir checkm2_out_dir
		mkdir tmp_dir

		checkm2 predict \
			--input "$derep_bins_dir" \
			--output-directory checkm2_out_dir \
			--extension fa \
			--threads ~{threads} \
			--database_path ~{checkm2_ref_db} \
			--remove_intermediates \
			--tmpdir tmp_dir
	>>>

	output {
		File bin_quality_report_tsv = "checkm2_out_dir/quality_report.tsv"
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

	Int disk_size = ceil(size(bin_quality_report_tsv, "GB") * 2 + size(filtered_contig_depth_txt, "GB") + 20)

	command <<<
		set -euo pipefail

		derep_bins_dir=$(dirname ~{derep_bins[0]})

		python /opt/scripts/Filter-Checkm2-Bins.py \
			--input_tsv ~{bin_quality_report_tsv} \
			--bin_dir "$derep_bins_dir" \
			--depth_file ~{filtered_contig_depth_txt} \
			--min_completeness ~{min_mag_completeness} \
			--max_contamination ~{max_mag_contamination} \
			--max_contigs ~{max_contigs} \
			--gtdb_outfile "~{sample_id}.GTDBTk_batch_file.txt" \
			--target_outfile "~{sample_id}.BinCount.txt" \
			--updated_tsv "~{sample_id}.quality_report.tsv"

		# Check if there are bins after CheckM2, before running GTDB-Tk and the MAG summary
		if [[ -s "~{sample_id}.BinCount.txt" ]]; then
			echo "true" > passed_bin_count_nonempty.txt
		else
			echo "false" > passed_bin_count_nonempty.txt
			echo "No bins passed filtering in CheckM2, see quality_report.tsv for more information"
		fi
	>>>

	output {
		File gtdb_batch_txt = "~{sample_id}.GTDBTk_batch_file.txt"
		File passed_bin_count_txt = "~{sample_id}.BinCount.txt"
		File filtered_quality_report_tsv = "~{sample_id}.quality_report.tsv"
		Boolean passed_bin_count_nonempty = read_boolean("passed_bin_count_nonempty.txt")
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
