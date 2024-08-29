version 1.0

# Harmonized human PMDBS bulk RNAseq workflow entrypoint

import "pmdbs_bulk_rnaseq_analysis_structs.wdl"
import "wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "preprocess/preprocess.wdl" as Preprocess
import "alignment_quantification/alignment_quantification.wdl" as AlignmentQuantification

workflow pmdbs_bulk_rnaseq_analysis {
	input {
		String cohort_id
		Array[Project] projects

		Boolean run_index_ref_genome = false
		ReferenceData? reference
		File? star_genome_dir_tar_gz

		# Cohort analysis
		Boolean run_cross_team_cohort_analysis = false
		String cohort_raw_data_bucket
		Array[String] cohort_staging_data_buckets

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	String workflow_execution_path = "workflow_execution"
	String workflow_name = "pmdbs_bulk_rnaseq"
	String workflow_version = "v1.0.0"
	String workflow_release = "https://github.com/ASAP-CRN/pmdbs-bulk-rnaseq-wf/releases/tag/pmdbs_bulk_rnaseq_analysis-~{workflow_version}"

	call GetWorkflowMetadata.get_workflow_metadata {
		input:
			zones = zones
	}

	scatter (project in projects) {
		# Add workflow name or else intermediate files will be uploaded to the same folder as harmonized_pmdbs
		String project_raw_data_path_prefix = "~{project.raw_data_bucket}/~{workflow_execution_path}/~{workflow_name}"

		call Preprocess.preprocess {
			input:
				project_id = project.project_id,
				samples = project.samples,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		Array[String] preprocessing_output_file_paths = flatten([
			preprocess.fastqc_reports_tar_gz,
			preprocess.trimmed_fastq_R1s,
			preprocess.trimmed_fastq_R2s,
			preprocess.failed_paired_fastq,
			preprocess.report_html,
			preprocess.trimmed_fastqc_reports_tar_gz
		]) #!StringCoercion

		call AlignmentQuantification.alignment_quantification {
			input:
				samples = project.samples,
				run_index_ref_genome = run_index_ref_genome,
				reference = select_first([reference]),
				star_genome_dir_tar_gz = select_first([star_genome_dir_tar_gz]),
				trimmed_fastq_R1s = preprocess.trimmed_fastq_R1s,
				trimmed_fastq_R2s = preprocess.trimmed_fastq_R2s,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}
	}

	output {
		# Sample-level outputs
		# Sample list
		Array[Array[Array[String]]] project_sample_ids = preprocess.project_sample_ids

		# Preprocess
		Array[Array[File]] fastqc_raw_reads_reports_tar_gz = preprocess.fastqc_reports_tar_gz
		Array[Array[Array[File]]] fastp_trimmed_fastq_R1s = preprocess.trimmed_fastq_R1s
		Array[Array[Array[File]]] fastp_trimmed_fastq_R2s = preprocess.trimmed_fastq_R2s
		Array[Array[Array[File]]] fastp_failed_paired_fastq = preprocess.failed_paired_fastq
		Array[Array[Array[File]]] fastp_report_html = preprocess.report_html
		Array[Array[File]] fastqc_trimmed_reads_reports_tar_gz = preprocess.trimmed_fastqc_reports_tar_gz

		# Alignment and quantification
		Array[Array[IndexData]] star_aligned_bams = alignment_quantification.aligned_bams
		Array[Array[File]] star_unmapped_mate1 = alignment_quantification.unmapped_mate1
		Array[Array[File]] star_unmapped_mate2 = alignment_quantification.unmapped_mate2
		Array[Array[File]] star_log = alignment_quantification.log
		Array[Array[File]] star_final_log = alignment_quantification.final_log
		Array[Array[File]] star_progress_log = alignment_quantification.progress_log
		Array[Array[File]] star_sj_out_tab = alignment_quantification.sj_out_tab
	}

	meta {
		description: "Harmonized human postmortem-derived brain sequencing (PMDBS) bulk RNA-seq workflow"
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files during cross-team cohort analysis."}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level cohort analysis."}
		run_index_ref_genome: {help: "Option to index reference genome with STAR. [false]"}
		reference: {help: "The primary assembly FASTA and gene annotation GTF from GENCODE."}
		star_genome_dir_tar_gz: {help: "The indexed reference genome files required for STAR."}
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (cellranger and generating the initial adata object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team cohort intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team cohort analysis outputs in."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in."}
	}
}
