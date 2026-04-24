version 1.0

# Downstream analysis including generating a single HTML report for multiple bioinformatics analyses across many samples and DGE analysis across multiple groups

import "../../wf-common/wdl/tasks/multiqc.wdl" as Multiqc
import "differential_gene_expression_analysis/differential_gene_expression_analysis.wdl" as DifferentialGeneExpressionAnalysis

workflow downstream {
	input {
		String team_id
		String dataset_id
		Array[Array[String]] project_sample_ids

		Array[File] output_files
		String output_name

		File metadata_csv
		File gene_map_csv
		File gene_ids_and_names_json

		String salmon_mode
		Array[File] salmon_quant_tar_gz

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	String sub_workflow_name = "downstream"
	String sub_workflow_version = "2.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{salmon_mode}/~{run_timestamp}"

	call Multiqc.multiqc {
		input:
			team_id = team_id,
			output_files = output_files,
			output_name = output_name,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call DifferentialGeneExpressionAnalysis.differential_gene_expression_analysis {
		input:
			team_id = team_id,
			dataset_id = dataset_id,
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
		# MultiQC report
		File multiqc_report_html = multiqc.multiqc_report_html #!FileCoercion
		File multiqc_data_zip = multiqc.multiqc_data_zip #!FileCoercion

		# PyDESeq2 DGE Analysis
		File dds_object_pkl = differential_gene_expression_analysis.dds_object_pkl #!FileCoercion
		Array[File] significant_genes_csv = differential_gene_expression_analysis.significant_genes_csv
		Array[File] volcano_plot_png = differential_gene_expression_analysis.volcano_plot_png
	}

	meta {
		description: "Aggregates upstream QC and alignment metrics into a MultiQC report and performs differential gene expression analysis using PyDESeq2."
	}

	parameter_meta {
		team_id: {help: "Name of the CRN Team; used to name output files."}
		dataset_id: {help: "Name of the ASAP-generated unique identifier for dataset; used to name output files."}
		project_sample_ids: {help: "Associated team ID, sample ID, and dataset DOI URL; used to generate a sample list."}
		output_files: {help: "Upstream output files to pass to MultiQC for report generation."}
	    output_name: {help: "Base name for the MultiQC report output file."}
	    metadata_csv: {help: "ASAP-generated CSV containing all sample information including batch, condition, etc. used for DESeq2 pairwise condition ('PD', 'Control'). For the `batch` column, there must be at least two distinct values."}
	    gene_map_csv: {help: "CSV containing mapped transcript IDs and gene IDs that must be in this order."}
	    gene_ids_and_names_json: {help: "JSON file containing mapped gene IDs and gene names created from the gene annotation GTF."}
	    salmon_mode: {help: "Salmon quantification mode; either 'alignment_mode' or 'mapping_mode'."}
		salmon_quant_tar_gz: {help: "Tar-gzipped Salmon quantification output directories, one per sample."}
		workflow_name: {help: "Workflow name; stored in the file-level manifest and final manifest with all saved files."}
		workflow_version: {help: "Workflow version; stored in the file-level manifest and final manifest with all saved files."}
		workflow_release: {help: "GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		run_timestamp: {help: "UTC timestamp; stored in the file-level manifest and final manifest with all saved files."}
		raw_data_path_prefix: {help: "Raw data bucket path prefix; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/downstream`)."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
