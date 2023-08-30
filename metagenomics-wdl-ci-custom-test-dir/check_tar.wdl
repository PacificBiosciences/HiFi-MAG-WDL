version 1.0

# Check contents of tar.gz file
# Input type: tar.gz file

task check_tar {
	input {
		File current_run_output
		File validated_output
	}

	Int disk_size = ceil(size(current_run_output, "GB") + size(validated_output, "GB") + 50)

	command <<<
		set -euo pipefail

		err() {
			message=$1

			echo -e "[ERROR] $message" >&2
		}

		validated_output_tar_contents=$(tar -tzf ~{validated_output} | sort)
		current_run_output_tar_contents=$(tar -tzf ~{current_run_output} | sort)

		if [[ "$validated_output_tar_contents" != "$current_run_output_tar_contents" ]]; then
			err "Contents in tar.gz do not match:
				Expected output: [$validated_output_tar_contents]
				Current run output: [$current_run_output_tar_contents]"
			exit 1
		else
			echo "Contents in tar.gz match: [$current_run_output_tar_contents]"
		fi
	>>>

	output {
	}

	runtime {
		docker: "ubuntu:xenial"
		cpu: 1
		memory: "3.75 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 1
	}
}
