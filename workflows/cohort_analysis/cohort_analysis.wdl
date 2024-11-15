version 1.0

# Identify overlapping significantly differentially expressed genes and create plots

import "../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "../../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[String] team_ids
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
	String staging_data_path_prefix = "~{workflow_name}"
	String upstream_staging_data_path = "~{staging_data_path_prefix}/upstream/~{salmon_mode}"
	String downstream_staging_data_path = "~{staging_data_path_prefix}/downstream/~{salmon_mode}"
	String cohort_analysis_staging_data_path = "~{staging_data_path_prefix}/~{sub_workflow_name}/~{salmon_mode}"

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
			team_ids = team_ids,
			n_teams = length(team_ids),
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
			staging_data_path = upstream_staging_data_path,
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
		select_all([
			degs_and_plot.overlapping_significant_genes_csv
		]),
		[
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

		# Overlapping differentially expressed genes only for cross_team_cohort_analysis
		File? overlapping_significant_genes_csv = degs_and_plot.overlapping_significant_genes_csv #!FileCoercion
		# PCA plots
		File pca_plot_png = degs_and_plot.pca_plot_png #!FileCoercion

		Array[File] upstream_manifest_tsvs = upload_upstream_files.manifests #!FileCoercion
		Array[File] downstream_manifest_tsvs = upload_downstream_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
	}
}

task degs_and_plot {
	input {
		String cohort_id
		Array[String] team_ids
		Int n_teams

		Array[File] significant_genes_csv
		Array[File] dds_object_pkl

		String salmon_mode

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones

		# Purposefully unset
		String? my_none
	}

	Int threads = 4
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(flatten([significant_genes_csv, dds_object_pkl]), "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/cohort_analysis.py \
			--cohort-id ~{cohort_id} \
			--project-ids ~{sep=' ' team_ids} \
			--n-teams ~{n_teams} \
			--degs ~{sep=' ' significant_genes_csv} \
			--dds-object ~{sep=' ' dds_object_pkl} \
			--salmon-mode ~{salmon_mode}

		if [[ ~{n_teams} -gt 1 ]]; then
			upload_outputs \
				-b ~{billing_project} \
				-d ~{raw_data_path} \
				-i ~{write_tsv(workflow_info)} \
				-o "~{cohort_id}.~{salmon_mode}.overlapping_significant_genes.csv"
		fi

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.~{salmon_mode}.pca_plot.png"
	>>>

	output {
		String? overlapping_significant_genes_csv = if (n_teams > 1) then "~{raw_data_path}/~{cohort_id}.~{salmon_mode}.overlapping_significant_genes.csv" else my_none
		String pca_plot_png = "~{raw_data_path}/~{cohort_id}.~{salmon_mode}.pca_plot.png"
	}
	runtime {
		docker: "~{container_registry}/pydeseq2:0.4.11"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}
