version 1.0

# Validate input FASTA files in tar.gz file
# Input type: FASTA file in tar.gz file

task fasta_validator_tar {
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

		validated_fa_filename=$(tar -tzf ~{validated_output} | sort | grep -m 1 '\.fa$' || [[ $? == 1 ]])

		current_run_fa_filename=$(tar -tzf ~{current_run_output} | sort | grep -m 1 '\.fa$' || [[ $? == 1 ]])

		mkdir validated_fas_dir
		validated_fa=$(tar -xvf ~{validated_output} --wildcards --no-anchored "$validated_fa_filename")
		mv "$validated_fa" validated_fas_dir
		validated_fa_basename=$(basename "$validated_fa")

		mkdir current_run_fas_dir
		current_run_fa=$(tar -xvf ~{current_run_output} --wildcards --no-anchored "$current_run_fa_filename")
		mv "$current_run_fa" current_run_fas_dir
		current_run_fa_basename=$(basename "$current_run_fa")

		# Checks both compressed and uncompressed fastas
		if ! py_fasta_validator -f validated_fas_dir/"$validated_fa_basename"; then
			err "Validated FASTA: [$validated_fa_basename] is invalid"
			exit 1
		else
			if ! py_fasta_validator -f current_run_fas_dir/"$current_run_fa_basename"; then
				err "Current run FASTA: [$current_run_fa_basename] is invalid"
				exit 1
			else
				echo "Current run FASTA: [$current_run_fa_basename] is valid"
			fi
		fi
	>>>

	output {
	}

	runtime {
		docker: "dnastack/dnastack-wdl-ci-tools:0.0.1"
		cpu: 1
		memory: "3.75 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 1
	}
}