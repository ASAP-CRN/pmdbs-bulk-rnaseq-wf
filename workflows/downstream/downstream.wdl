version 1.0

# Downstream analysis including generating a single HTML report for multiple bioinformatics analyses across many samples and DGE analysis across multiple groups

import "../wf-common/wdl/tasks/multiqc.wdl" as Multiqc
import "differential_gene_expression_analysis/differential_gene_expression_analysis.wdl" as DifferentialGeneExpressionAnalysis

workflow downstream {
	input {
		String project_id

		Array[File] output_files
		String output_name

		File metadata_csv
		File gene_map_csv
		File blacklist_genes_bed
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
	String sub_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{salmon_mode}/~{run_timestamp}"

	call Multiqc.multiqc {
		input:
			project_id = project_id,
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
			project_id = project_id,
			metadata_csv = metadata_csv,
			gene_map_csv = gene_map_csv,
			blacklist_genes_bed = blacklist_genes_bed,
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
		File significant_genes_csv = differential_gene_expression_analysis.significant_genes_csv #!FileCoercion
		File volcano_plot_png = differential_gene_expression_analysis.volcano_plot_png #!FileCoercion
	}
}
