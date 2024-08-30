version 1.0

# Map and quantify reads (mapping-based mode)

import "../pmdbs_bulk_rnaseq_analysis_structs.wdl"

workflow pseudo_mapping_quantification {
	input {
		Array[Sample] samples

		Boolean run_index_ref_genome
		ReferenceData reference
		File salmon_genome_dir_tar_gz

		Array[File] trimmed_fastq_R1s
		Array[File] trimmed_fastq_R2s

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	# Task and subworkflow versions
	String sub_workflow_name = "pseudo_mapping_quantification"
	String mapping_quantification_task_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String salmon_raw_data_path = "~{workflow_raw_data_path_prefix}/mapping_quantification/~{mapping_quantification_task_version}"

	scatter (sample_object in samples) {
		String mapping_quantification_output = "~{salmon_raw_data_path}/~{sample_object.sample_id}.mapping_mode.quant.sf"
	}

	# For each sample, outputs an array of true/false: [mapping_quantification_complete]
	call check_output_files_exist {
		input:
			mapping_quantification_output_files = mapping_quantification_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (sample_index in range(length(samples))) {
		Sample sample = samples[sample_index]

		String mapping_quantification_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][0]

		if (run_index_ref_genome) {
			call generate_decoy {
			input:
				primary_assembly_fasta = reference.primary_assembly_fasta,
				transcripts_fasta = reference.transcripts_fasta,
				container_registry = container_registry,
				zones = zones
			}

			call index_ref_genome {
			input:
				gentrome_fasta = generate_decoy.gentrome_fasta,
				decoys_txt = generate_decoy.decoys_txt,
				container_registry = container_registry,
				zones = zones
			}
		}

		File genome_dir_tar_gz = select_first([
			index_ref_genome.salmon_genome_dir_tar_gz,
			salmon_genome_dir_tar_gz
		])

		String salmon_quant_file = "~{salmon_raw_data_path}/~{sample.sample_id}.mapping_mode.quant.sf"
		String salmon_command_info_json = "~{salmon_raw_data_path}/~{sample.sample_id}.mapping_mode.cmd_info.json"
		String salmon_aux_info_tar_gz = "~{salmon_raw_data_path}/~{sample.sample_id}.mapping_mode.aux_info.tar.gz"

		if (mapping_quantification_complete == "false") {
			call mapping_quantification {
				input:
					sample_id = sample.sample_id,
					salmon_genome_dir_tar_gz = genome_dir_tar_gz,
					trimmed_fastq_R1s = trimmed_fastq_R1s,
					trimmed_fastq_R2s = trimmed_fastq_R2s,
					raw_data_path = salmon_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File quant_file_output = select_first([mapping_quantification.quant_file, salmon_quant_file]) #!FileCoercion
		File command_info_json_output = select_first([mapping_quantification.command_info_json, salmon_command_info_json]) #!FileCoercion
		File aux_info_tar_gz_output = select_first([mapping_quantification.aux_info_tar_gz, salmon_aux_info_tar_gz]) #!FileCoercion
	}

	output {
		# Salmon mapping and quantification
		Array[File] quant_file = quant_file_output #!FileCoercion
		Array[File] command_info_json = command_info_json_output #!FileCoercion
		Array[File] aux_info_tar_gz = aux_info_tar_gz_output #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] mapping_quantification_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			mapping_quantification_file=$(echo "${output_files}" | cut -f 1)

			if gsutil -u ~{billing_project} ls "${mapping_quantification_file}"; then
				# If we find all outputs, don't rerun anything
				echo -e "true" >> sample_preprocessing_complete.tsv
			else
				# If we don't find mapping_quantification outputs, run (or rerun)
				echo -e "false" >> sample_preprocessing_complete.tsv
			fi
		done < ~{write_lines(mapping_quantification_output_files)}
	>>>

	output {
		Array[Array[String]] sample_preprocessing_complete = read_tsv("sample_preprocessing_complete.tsv")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
		zones: zones
	}
}

task generate_decoy {
	input {
		File primary_assembly_fasta
		File transcripts_fasta
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([primary_assembly_fasta, transcripts_fasta], "GB") *2 + 100)

	command <<<
		set -euo pipefail

		grep "^>" <(gunzip -c ~{primary_assembly_fasta}) | cut -d " " -f 1 > decoys.txt
		sed -i.bak -e 's/>//g' decoys.txt
		cat ~{transcripts_fasta} ~{primary_assembly_fasta} > gentrome.fa.gz
	>>>

	output {
		File decoys_txt = "decoys.txt"
		File gentrome_fasta = "gentrome.fa.gz"
	}

	runtime {
		docker: "~{container_registry}/salmon:1.10.1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}

task index_ref_genome {
	input {
		File gentrome_fasta
		File decoys_txt
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([gentrome_fasta, decoys_txt], "GB") *2 + 100)

	command <<<
		set -euo pipefail

		salmon index \
			--transcripts ~{gentrome_fasta} \
			--index salmon_index \
			--kmerLen 31 \
			--gencode \
			--decoys ~{decoys_txt} \
			--threads ~{threads}

		tar -czvf salmon_genome_dir.tar.gz salmon_index
	>>>

	output {
		File salmon_genome_dir_tar_gz = "salmon_genome_dir.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/salmon:1.10.1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}

task mapping_quantification {
	input {
		String sample_id

		File salmon_genome_dir_tar_gz

		Array[File] trimmed_fastq_R1s
		Array[File] trimmed_fastq_R2s

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(salmon_genome_dir_tar_gz, "GB") + size(flatten([trimmed_fastq_R1s, trimmed_fastq_R2s])) * 2 + 50)

	command <<<
		set -euo pipefail

		tar -xzvf ~{salmon_genome_dir_tar_gz}

		salmon quant \
			--index salmon_index \
			--libType A \
			--mates1 ~{sep=' ' trimmed_fastq_R1s} \
			--mates2 ~{sep=' ' trimmed_fastq_R2s} \
			--output salmon_quant \
			--validateMappings \
			--threads ~{threads}

		cd salmon_quant
		mv quant.sf ~{sample_id}.mapping_mode.quant.sf
		mv cmd_info.json ~{sample_id}.mapping_mode.cmd_info.json
		tar -czvf "~{sample_id}.mapping_mode.aux_info.tar.gz" "aux_info"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.mapping_mode.quant.sf" \
			-o "~{sample_id}.mapping_mode.cmd_info.json" \
			-o "~{sample_id}.mapping_mode.aux_info.tar.gz"
	>>>

	output {
		String quant_file = "~{raw_data_path}/~{sample_id}.mapping_mode.quant.sf"
		String command_info_json = "~{raw_data_path}/~{sample_id}.mapping_mode.cmd_info.json"
		String aux_info_tar_gz = "~{raw_data_path}/~{sample_id}.mapping_mode.aux_info.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/salmon:1.10.1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
