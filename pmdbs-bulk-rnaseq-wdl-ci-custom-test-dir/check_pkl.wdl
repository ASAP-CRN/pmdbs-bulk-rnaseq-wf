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

		validate_pkl() {
			test_file=$1

			python3.9 -c "import pickle; with open('$test_file', 'rb') as f: data = pickle.load(f)"
		}

		if ! (validate_pkl ~{validated_output}); then
			err "Validated output file [~{basename(validated_output)}] is not a valid pkl file"
			exit 1
		fi
		
		if (validate_pkl ~{current_run_output}); then
			echo "Current run output [~{basename(current_run_output)}] is a valid pkl file"
		else
			err "Current run output [~{basename(current_run_output)}] is not a valid pkl file"
			exit 1
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
