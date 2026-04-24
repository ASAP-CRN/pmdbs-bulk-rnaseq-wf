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

		Array[Array[File]] significant_genes_csv
		Array[File] dds_object_pkl

		String salmon_mode

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String crn_release_version
		String raw_data_path_prefix
		Array[String] staging_data_buckets
		String billing_project
		String container_registry
		String zones
	}

	String sub_workflow_name = "cohort_analysis"
	String sub_workflow_version = "2.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{salmon_mode}/~{run_timestamp}"
	String staging_data_path_prefix = "~{workflow_name}/release/~{crn_release_version}"
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
			significant_genes_csv = flatten(significant_genes_csv),
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
		flatten(
			select_all([
				degs_and_plot.overlapping_significant_genes_csv
			])
		),
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

		# Overlapping differentially expressed genes per contrast, only for cross_team_cohort_analysis
		Array[File]? overlapping_significant_genes_csv = degs_and_plot.overlapping_significant_genes_csv #!FileCoercion
		# PCA plots
		File pca_plot_png = degs_and_plot.pca_plot_png #!FileCoercion

		Array[File] upstream_manifest_tsvs = upload_upstream_files.manifests #!FileCoercion
		Array[File] downstream_manifest_tsvs = upload_downstream_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
	}

	meta {
		description: "Identifies overlapping differentially expressed genes across teams and generates a PCA plot from PyDESeq2 objects."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		team_ids: {help: "Array of CRN Teams included in cohort analysis."}
		project_sample_ids: {help: "Associated team ID, sample ID, and dataset DOI URL; used to generate a sample list."}
		upstream_output_file_paths: {help: "Selected upstream output files to upload to the staging bucket alongside selected cohort analysis output files."}
		downstream_output_file_paths: {help: "Selected downstream output files to upload to the staging bucket alongside selected cohort analysis output files."}
		significant_genes_csv: {help: "Per-team, per-contrast CSV files of significantly differentially expressed genes from PyDESeq2."}
    	dds_object_pkl: {help: "Per-team pickled PyDESeq2 dataset objects used for cross-team PCA analysis."}
	    salmon_mode: {help: "Salmon quantification mode; either 'alignment_mode' or 'mapping_mode'."}
		workflow_name: {help: "Workflow name; stored in the file-level manifest and final manifest with all saved files."}
		workflow_version: {help: "Workflow version; stored in the file-level manifest and final manifest with all saved files."}
		workflow_release: {help: "GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		run_timestamp: {help: "UTC timestamp; stored in the file-level manifest and final manifest with all saved files."}
		crn_release_version: {help: "CRN Cloud release version; used to organize outputs and for the CRN Cloud release."}
		raw_data_path_prefix: {help: "Raw data bucket path prefix; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis`)."}
		staging_data_buckets: {help: "Array of staging data buckets to upload intermediate files to (i.e., DEV or UAT buckets depending on internal QC status)."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
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
		Array[String]? my_none
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
			for f in *.overlapping_significant_genes.csv; do
				upload_outputs \
					-b ~{billing_project} \
					-d ~{raw_data_path} \
					-i ~{write_tsv(workflow_info)} \
					-o "$f"
				echo "~{raw_data_path}/$f" >> overlapping_significant_genes_csv_paths.txt
			done
		fi

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.~{salmon_mode}.pca_plot.png"
	>>>

	output {
		Array[String]? overlapping_significant_genes_csv = if (n_teams > 1) then read_lines("overlapping_significant_genes_csv_paths.txt") else my_none
		String pca_plot_png = "~{raw_data_path}/~{cohort_id}.~{salmon_mode}.pca_plot.png"
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
		description: "Identifies overlapping significantly differentially expressed genes across teams and generates a cross-team PCA plot from PyDESeq2 objects."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		team_ids: {help: "Array of CRN Teams included in cohort analysis."}
		n_teams: {help: "Number of CRN Teams in the cohort; overlapping DEG analysis only runs when greater than 1."}
		significant_genes_csv: {help: "Per-team, per-contrast CSV files of significantly differentially expressed genes from PyDESeq2."}
    	dds_object_pkl: {help: "Per-team pickled PyDESeq2 dataset objects used for cross-team PCA analysis."}
	    salmon_mode: {help: "Salmon quantification mode; either 'alignment_mode' or 'mapping_mode'."}
		raw_data_path: {help: "Raw data bucket path for DGE outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<salmon_mode>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
