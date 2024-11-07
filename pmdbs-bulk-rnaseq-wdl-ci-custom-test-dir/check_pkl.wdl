version 1.0

# Validate input PKL files
# Input type: PKL

task check_pkl {
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

		current_run_pkl_magic_number=$(head -c 2 ~{current_run_output} | od -An -N2 -tx1 | grep "80  04" || [[ $? == 1 ]])
		validated_pkl_magic_number=$(head -c 2 ~{validated_output} | od -An -N2 -tx1 | grep "80  04" || [[ $? == 1 ]])

		if [[ -z "$validated_pkl_magic_number" ]]; then
			err "Validated output file [~{basename(validated_output)}] is not a valid pkl file"
			exit 1
		fi
		
		if [[ -z "$current_run_pkl_magic_number" ]]; then
			err "Current run output [~{basename(current_run_output)}] is not a valid pkl file"
			exit 1
		else
			echo "Current run output [~{basename(current_run_output)}] is a valid pkl file"
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
