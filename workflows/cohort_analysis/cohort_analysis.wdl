version 1.0

# Identify overlapping significantly differentially expressed genes and create plots

import "../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] upstream_output_file_paths = []
		Array[String] downstream_output_file_paths = []

		Array[File] significant_genes_csv
		Array[File] dds_object_pkl

		String salmon_mode

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

	String sub_workflow_name = "cohort_analysis"
	String sub_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{salmon_mode}/~{run_timestamp}"
	String downstream_staging_data_path = "downstream/~{salmon_mode}"
	String cohort_analysis_staging_data_path = "~{sub_workflow_name}/~{salmon_mode}"

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

	call degs_and_plot {
		input:
			cohort_id = cohort_id,
			cohort_sample_list = write_cohort_sample_list.cohort_sample_list, #!FileCoercion
			significant_genes_csv = significant_genes_csv,
			dds_object_pkl = dds_object_pkl,
			salmon_mode = salmon_mode,
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

	call UploadFinalOutputs.upload_final_outputs as upload_downstream_files {
		input:
			output_file_paths = downstream_output_file_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = downstream_staging_data_path,
			billing_project = billing_project,
			zones = zones
	}

	Array[String] cohort_analysis_final_output_paths = flatten([
		[
			write_cohort_sample_list.cohort_sample_list
		],
		[
			select_first([degs_and_plot.overlapping_significant_genes_csv]),
			degs_and_plot.pca_plot_png
		]
	]) #!StringCoercion

	call UploadFinalOutputs.upload_final_outputs as upload_cohort_analysis_files {
		input:
			output_file_paths = cohort_analysis_final_output_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = cohort_analysis_staging_data_path,
			billing_project = billing_project,
			zones = zones
	}

	output {
		File cohort_sample_list = write_cohort_sample_list.cohort_sample_list #!FileCoercion

		# Overlapping differentially expressed genes
		File? overlapping_significant_genes_csv = degs_and_plot.overlapping_significant_genes_csv #!FileCoercion
		# PCA plot
		File pca_plot_png = degs_and_plot.pca_plot_png #!FileCoercion

		Array[File] upstream_manifest_tsvs = upload_upstream_files.manifests #!FileCoercion
		Array[File] downstream_manifest_tsvs = upload_downstream_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
	}
}

task degs_and_plot {
	input {
		String cohort_id
		File cohort_sample_list

		Array[File] significant_genes_csv
		Array[File] dds_object_pkl

		String salmon_mode

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil((size(cohort_sample_list, "GB") + size(flatten([significant_genes_csv, dds_object_pkl]), "GB")) * 2 + 50)

	command <<<
		set -euo pipefail

		project_id_column_no=$(tr '\t' '\n' < <(head -1 ~{cohort_sample_list}) | grep -n "sample_id" | cut -d : -f 1)
		project_ids=$(cut -f "$project_id_column_no" ~{cohort_sample_list} | sed -r '/^\s*$/d' | tail -n +2 | uniq | tr '\n' ' ') # Do not sort to preserve order

		python3 /opt/scripts/cohort_analysis.py \
			--cohort-id ~{cohort_id} \
			--project-id "${project_ids}" \
			--degs ~{sep=' ' significant_genes_csv} \
			--dds-object ~{sep=' ' dds_object_pkl} \
			--salmon-mode ~{salmon_mode}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.~{salmon_mode}.pca_plot.png"

		project_id_length=$(echo "$project_ids" | wc -w)
		if [[ "$project_id_length" -gt 1 ]]; then
			upload_outputs \
				-b ~{billing_project} \
				-d ~{raw_data_path} \
				-i ~{write_tsv(workflow_info)} \
				-o "~{cohort_id}.~{salmon_mode}.overlapping_significant_genes.csv"
		fi
	>>>

	output {
		String? overlapping_significant_genes_csv = "~{raw_data_path}/~{cohort_id}.~{salmon_mode}.overlapping_significant_genes.csv" #!UnnecessaryQuantifier
		String pca_plot_png = "~{raw_data_path}/~{cohort_id}.~{salmon_mode}.pca_plot.png"
	}
	runtime {
		docker: "~{container_registry}/pydeseq2:0.4.11"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
