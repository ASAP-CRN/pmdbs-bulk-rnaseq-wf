version 1.0

# Map and quantify reads (mapping-based mode)

workflow pseudo_mapping_quantification {
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

	call mapping_quantification {
		input:
			sample_id = sample_id,
			salmon_genome_dir_tar_gz = salmon_genome_dir_tar_gz,
			trimmed_fastq_R1s = trimmed_fastq_R1s,
			trimmed_fastq_R2s = trimmed_fastq_R2s,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		# Salmon mapping and quantification
		File quant_tar_gz = mapping_quantification.quant_tar_gz #!FileCoercion
	}

	meta {
		description: "Maps and quantifies trimmed reads directly against the transcriptome reference using Salmon in mapping-based mode."
	}

	parameter_meta {
		sample_id: {help: "Generated ASAP sample ID; used to name output files."}
		salmon_genome_dir_tar_gz: {help: "The indexed concatenated transcriptome and genome files required for Salmon."}
		trimmed_fastq_R1s: {help: "Adapter-trimmed forward (R1) FASTQ files for the sample."}
    	trimmed_fastq_R2s: {help: "Adapter-trimmed reverse (R2) FASTQ files for the sample."}
		raw_data_path: {help: "Raw data bucket path for pseudo-mapping and quantification outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/upstream/pseudo_mapping_quantification`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
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
	Int disk_size = ceil((size(salmon_genome_dir_tar_gz, "GB") + size(flatten([trimmed_fastq_R1s, trimmed_fastq_R2s]), "GB")) * 2 + 50)

	command <<<
		set -euo pipefail

		tar -xzvf ~{salmon_genome_dir_tar_gz}

		salmon quant \
			--index salmon_index \
			--libType A \
			--mates1 ~{sep=' ' trimmed_fastq_R1s} \
			--mates2 ~{sep=' ' trimmed_fastq_R2s} \
			--output ~{sample_id}_salmon_quant \
			--validateMappings \
			--threads ~{threads} \
			--rangeFactorizationBins 4 \
			--gcBias

		# Outputs must remain in folder and unmodified for downstream analysis
		# Outputs include: quant.sf, cmd_info.json, and aux_info folder
		tar -czvf "~{sample_id}.mapping_mode.salmon_quant.tar.gz" "~{sample_id}_salmon_quant"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.mapping_mode.salmon_quant.tar.gz"
	>>>

	output {
		String quant_tar_gz = "~{raw_data_path}/~{sample_id}.mapping_mode.salmon_quant.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/salmon:1.10.3"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}

	meta {
		description: "Maps trimmed paired-end reads to the Salmon index and quantifies transcript abundances using Salmon in mapping-based mode."
	}

	parameter_meta {
		sample_id: {help: "Generated ASAP sample ID; used to name output files."}
		salmon_genome_dir_tar_gz: {help: "The indexed concatenated transcriptome and genome files required for Salmon."}
		trimmed_fastq_R1s: {help: "Adapter-trimmed forward (R1) FASTQ files for the sample."}
    	trimmed_fastq_R2s: {help: "Adapter-trimmed reverse (R2) FASTQ files for the sample."}
		raw_data_path: {help: "Raw data bucket path for pseudo-mapping and quantification outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/upstream/pseudo_mapping_quantification/<pseudo_mapping_quantification_workflow_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
