version 1.0

# Downstream analysis including DGE analysis across multiple groups

import "../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "differential_gene_expression_analysis/differential_gene_expression_analysis.wdl" as DifferentialGeneExpressionAnalysis
import "../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow downstream {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] upstream_output_file_paths = []

		File metadata_csv
		File gene_map_csv
		File blacklist_genes_bed
		File cell_type_markers_list #TODO

		String salmon_mode
		Array[File] salmon_quant_tar_gz

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		Array[String] staging_data_buckets
		String billing_project
		String container_registry
		String zones
	}

	String sub_workflow_name = "downstream"
	String sub_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{run_timestamp}"

	call WriteCohortSampleList.write_cohort_sample_list {
		input:
			cohort_id = cohort_id,
			project_sample_ids = project_sample_ids,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call DifferentialGeneExpressionAnalysis.differential_gene_expression_analysis {
		input:
			cohort_id = cohort_id,
			cohort_sample_list = write_cohort_sample_list.cohort_sample_list, #!FileCoercion
			metadata_csv = metadata_csv,
			gene_map_csv = gene_map_csv,
			blacklist_genes_bed = blacklist_genes_bed,
			salmon_mode = salmon_mode,
			salmon_quant_tar_gz = salmon_quant_tar_gz,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call UploadFinalOutputs.upload_final_outputs as upload_upstream_files {
		input:
			output_file_paths = upstream_output_file_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = "upstream",
			billing_project = billing_project,
			zones = zones
	}

	Array[String] downstream_final_output_paths = flatten([
		[
			differential_gene_expression_analysis.significant_genes_csv,
			differential_gene_expression_analysis.pca_plot_png,
			differential_gene_expression_analysis.volcano_plot_png
		]
	]) #!StringCoercion

	call UploadFinalOutputs.upload_final_outputs as upload_downstream_files {
		input:
			output_file_paths = downstream_final_output_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = sub_workflow_name,
			billing_project = billing_project,
			zones = zones
	}

	output {
		File cohort_sample_list = write_cohort_sample_list.cohort_sample_list #!FileCoercion

		# PyDESeq2 DGE Analysis
		File significant_genes_csv = differential_gene_expression_analysis.significant_genes_csv #!FileCoercion
		File pca_plot_png = differential_gene_expression_analysis.pca_plot_png #!FileCoercion
		File volcano_plot_png = differential_gene_expression_analysis.volcano_plot_png #!FileCoercion

		Array[File] upstream_manifest_tsvs = upload_upstream_files.manifests #!FileCoercion
		Array[File] downstream_manifest_tsvs = upload_downstream_files.manifests #!FileCoercion
	}
}
