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

		current_run_pkl_file=~{current_run_output}
		validated_pkl_file=~{validated_output}

		python3 - <<EOF
		import pickle

		try:
			with open("$validated_pkl_file", "rb") as f:
				data = pickle.load(f)
			print("Validated PKL file loaded successfully")
		except (pickle.UnpicklingError, EOFError, FileNotFoundError) as e:
			print("Error loading validated PKL file:", e)
			exit(1)

		try:
			with open("$current_run_pkl_file", "rb") as f:
				data = pickle.load(f)
			print("Current run PKL file loaded successfully")
		except (pickle.UnpicklingError, EOFError, FileNotFoundError) as e:
			print("Error loading current run PKL file:", e)
			exit(1)
		EOF

		if [[ $? -ne 0 ]]; then
			err "PKL file failed validation: 
				Expected output: [$validated_pkl_file]
				Current run output: [$current_run_pkl_file]"
			exit 1
		else
			echo "PKL file succeeded validation:
				Expected output: [$validated_pkl_file]
				Current run output: [$current_run_pkl_file]"
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
