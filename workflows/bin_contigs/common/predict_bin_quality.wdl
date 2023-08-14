version 1.0

import "../../wdl-common/wdl/structs.wdl"

task predict_bin_quality {
	input {
		String prefix
		File checkm2_ref_db

		Array[File] bin_fas

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = threads * 2
	Int disk_size = ceil((size(checkm2_ref_db, "GB") + (size(bin_fas[0], "GB") * length(bin_fas))) * 2 + 20)

	command <<<
		set -euo pipefail

		checkm2 --version

		# Ensure all bins are in the bin_dir
		mkdir tmp_dir bin_dir
		while read -r bin || [[ -n "${bin}" ]]; do
			ln -s "${bin}" "$(pwd)/bin_dir"
		done < ~{write_lines(bin_fas)}

		checkm2 predict \
			--input bin_dir \
			--output-directory checkm2_out_dir \
			--extension fa \
			--threads ~{threads} \
			--database_path ~{checkm2_ref_db} \
			--remove_intermediates \
			--tmpdir tmp_dir

		mv checkm2_out_dir/quality_report.tsv "checkm2_out_dir/~{prefix}.quality_report.tsv"
	>>>

	output {
		File bin_quality_report_tsv = "checkm2_out_dir/~{prefix}.quality_report.tsv"
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
