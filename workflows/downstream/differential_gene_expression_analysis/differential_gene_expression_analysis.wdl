version 1.0

# Differential gene expression analysis based on the negative binomial distribution

workflow differential_gene_expression_analysis {
	input {
		String team_id
		Array[Array[String]] project_sample_ids
		
		File metadata_csv
		File gene_map_csv
		File gene_ids_and_names_json

		String salmon_mode
		Array[File] salmon_quant_tar_gz

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	call differential_gene_expression {
		input:
			team_id = team_id,
			project_sample_ids = project_sample_ids,
			metadata_csv = metadata_csv,
			gene_map_csv = gene_map_csv,
			gene_ids_and_names_json = gene_ids_and_names_json,
			salmon_mode = salmon_mode,
			salmon_quant_tar_gz = salmon_quant_tar_gz,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		# PyDESeq2 DGE Analysis
		File dds_object_pkl = differential_gene_expression.dds_object_pkl #!FileCoercion
		Array[File] significant_genes_csv = differential_gene_expression.significant_genes_csv #!FileCoercion
		Array[File] volcano_plot_png = differential_gene_expression.volcano_plot_png #!FileCoercion
	}

	meta {
		description: "Performs differential gene expression analysis on Salmon quantification outputs using PyDESeq2 by comparing non-controls and controls."
	}

	parameter_meta {
		team_id: {help: "Name of the CRN Team; used to name output files."}
		project_sample_ids: {help: "Associated team ID, sample ID, and dataset DOI URL; used to generate a sample list."}
	    metadata_csv: {help: "ASAP-generated CSV containing all sample information including batch, condition, etc. used for DESeq2 pairwise condition ('PD', 'Control'). For the `batch` column, there must be at least two distinct values."}
	    gene_map_csv: {help: "CSV containing mapped transcript IDs and gene IDs that must be in this order."}
	    gene_ids_and_names_json: {help: "JSON file containing mapped gene IDs and gene names created from the gene annotation GTF."}
	    salmon_mode: {help: "Salmon quantification mode; either 'alignment_mode' or 'mapping_mode'."}
		salmon_quant_tar_gz: {help: "Tar-gzipped Salmon quantification output directories, one per sample."}
		raw_data_path: {help: "Raw data bucket path for DGE outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/downstream/<downstream_version>/<salmon_mode>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task differential_gene_expression {
	input {
		String team_id
		Array[Array[String]] project_sample_ids
		
		File metadata_csv
		File gene_map_csv
		File gene_ids_and_names_json

		String salmon_mode
		Array[File] salmon_quant_tar_gz

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 4
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil((size([metadata_csv, gene_map_csv], "GB") + size(flatten([salmon_quant_tar_gz]), "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		while read -r quant_tar_gz || [[ -n "${quant_tar_gz}" ]]; do
			tar -xzvf "${quant_tar_gz}"
		done < ~{write_lines(salmon_quant_tar_gz)}

		python3 /opt/scripts/dge_analysis.py \
			--team-id ~{team_id} \
			--sample-ids ~{write_tsv(project_sample_ids)} \
			--metadata ~{metadata_csv} \
			--gene-map ~{gene_map_csv} \
			--gene-ids-and-names ~{gene_ids_and_names_json} \
			--salmon-mode ~{salmon_mode}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{team_id}.~{salmon_mode}.dds.pkl"

		# Upload per-contrast significant genes CSVs and volcano plots
		for f in *.pydeseq2_significant_genes.csv; do
			upload_outputs \
				-b ~{billing_project} \
				-d ~{raw_data_path} \
				-i ~{write_tsv(workflow_info)} \
				-o "$f"
			echo "~{raw_data_path}/$f" >> significant_genes_csv_paths.txt
		done

		for f in *.volcano_plot.png; do
			upload_outputs \
				-b ~{billing_project} \
				-d ~{raw_data_path} \
				-i ~{write_tsv(workflow_info)} \
				-o "$f"
			echo "~{raw_data_path}/$f" >> volcano_plot_png_paths.txt
		done
	>>>

	output {
		String dds_object_pkl = "~{raw_data_path}/~{team_id}.~{salmon_mode}.dds.pkl"
		Array[String] significant_genes_csv = read_lines("significant_genes_csv_paths.txt")
		Array[String] volcano_plot_png = read_lines("volcano_plot_png_paths.txt")
	}
	runtime {
		docker: "~{container_registry}/pydeseq2:0.5.2_1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}

	meta {
		description: "Performs differential gene expression analysis on Salmon quantification outputs using PyDESeq2 by comparing non-controls and controls."
	}

	parameter_meta {
		team_id: {help: "Name of the CRN Team; used to name output files."}
		project_sample_ids: {help: "Associated team ID, sample ID, and dataset DOI URL; used to generate a sample list."}
	    metadata_csv: {help: "ASAP-generated CSV containing all sample information including batch, condition, etc. used for DESeq2 pairwise condition ('PD', 'Control'). For the `batch` column, there must be at least two distinct values."}
	    gene_map_csv: {help: "CSV containing mapped transcript IDs and gene IDs that must be in this order."}
	    gene_ids_and_names_json: {help: "JSON file containing mapped gene IDs and gene names created from the gene annotation GTF."}
	    salmon_mode: {help: "Salmon quantification mode; either 'alignment_mode' or 'mapping_mode'."}
		salmon_quant_tar_gz: {help: "Tar-gzipped Salmon quantification output directories, one per sample."}
		raw_data_path: {help: "Raw data bucket path for DGE outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/downstream/<downstream_version>/<salmon_mode>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
