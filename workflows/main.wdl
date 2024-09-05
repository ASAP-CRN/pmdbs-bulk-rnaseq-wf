version 1.0

# Harmonized human PMDBS bulk RNAseq workflow entrypoint

import "pmdbs_bulk_rnaseq_analysis_structs.wdl"
import "wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "preprocess/preprocess.wdl" as Preprocess
import "alignment_quantification/alignment_quantification.wdl" as AlignmentQuantification
import "pseudo_mapping_quantification/pseudo_mapping_quantification.wdl" as PseudoMappingQuantification

workflow pmdbs_bulk_rnaseq_analysis {
	input {
		String cohort_id
		Array[Project] projects

		ReferenceData? reference

		# STAR and salmon quantification
		Boolean run_alignment_quantification = true
		Boolean run_star_index_ref_genome = false
		File? star_genome_dir_tar_gz

		# Salmon mapping and quantification
		Boolean run_pseudo_mapping_quantification = false
		Boolean run_salmon_index_ref_genome = false
		File? salmon_genome_dir_tar_gz

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
			flatten(preprocess.trimmed_fastq_R1s),
			flatten(preprocess.trimmed_fastq_R2s),
			preprocess.failed_paired_fastqs_tar_gz,
			preprocess.reports_html_tar_gz,
			preprocess.trimmed_fastqc_reports_tar_gz
		]) #!StringCoercion

		if (run_alignment_quantification) {
			call AlignmentQuantification.alignment_quantification {
				input:
					trimmed_samples = preprocess.trimmed_samples,
					run_index_ref_genome = run_star_index_ref_genome,
					reference = select_first([reference]),
					star_genome_dir_tar_gz = select_first([star_genome_dir_tar_gz]),
					workflow_name = workflow_name,
					workflow_version = workflow_version,
					workflow_release = workflow_release,
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = project_raw_data_path_prefix,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					zones = zones
			}

			Array[String] alignment_quantification_output_file_paths = flatten([
				alignment_quantification.aligned_bam,
				alignment_quantification.aligned_bam_index,
				alignment_quantification.unmapped_mate1,
				alignment_quantification.unmapped_mate2,
				alignment_quantification.log,
				alignment_quantification.final_log,
				alignment_quantification.progress_log,
				alignment_quantification.sj_out_tab,
				alignment_quantification.quant_tar_gz
			]) #!StringCoercion
		}

		if (run_pseudo_mapping_quantification) {
			call PseudoMappingQuantification.pseudo_mapping_quantification {
				input:
					trimmed_samples = preprocess.trimmed_samples,
					run_index_ref_genome = run_salmon_index_ref_genome,
					reference = select_first([reference]),
					salmon_genome_dir_tar_gz = select_first([salmon_genome_dir_tar_gz]),
					workflow_name = workflow_name,
					workflow_version = workflow_version,
					workflow_release = workflow_release,
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = project_raw_data_path_prefix,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					zones = zones
			}

			Array[String] pseudo_mapping_quantification_output_file_paths = flatten([
				pseudo_mapping_quantification.quant_tar_gz
			]) #!StringCoercion
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
		Array[Array[File]] fastp_failed_paired_fastqs_tar_gz = preprocess.failed_paired_fastqs_tar_gz
		Array[Array[File]] fastp_report_html_tar_gz = preprocess.reports_html_tar_gz
		Array[Array[File]] fastqc_trimmed_reads_reports_tar_gz = preprocess.trimmed_fastqc_reports_tar_gz

		# Alignment and quantification
		Array[Array[File]?] star_aligned_bam = alignment_quantification.aligned_bam
		Array[Array[File]?] star_aligned_bam_index = alignment_quantification.aligned_bam_index
		Array[Array[File]?] star_unmapped_mate1 = alignment_quantification.unmapped_mate1
		Array[Array[File]?] star_unmapped_mate2 = alignment_quantification.unmapped_mate2
		Array[Array[File]?] star_log = alignment_quantification.log
		Array[Array[File]?] star_final_log = alignment_quantification.final_log
		Array[Array[File]?] star_progress_log = alignment_quantification.progress_log
		Array[Array[File]?] star_sj_out_tab = alignment_quantification.sj_out_tab
		Array[Array[File]?] salmon_alignment_mode_quant_tar_gz = alignment_quantification.quant_tar_gz

		# Direct quantification with mapping
		Array[Array[File]?] salmon_mapping_mode_quant_tar_gz = pseudo_mapping_quantification.quant_tar_gz
	}

	meta {
		description: "Harmonized human postmortem-derived brain sequencing (PMDBS) bulk RNA-seq workflow"
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files during cross-team cohort analysis."}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level cohort analysis."}
		reference: {help: "The primary assembly FASTA and gene annotation GTF from GENCODE."}
		run_alignment_quantification: {help: "Option to align raw reads with STAR and quantify aligned reads with Salmon. This and/or '' must be set to true. [true]"}
		run_star_index_ref_genome: {help: "Option to index reference genome with STAR. [false]"}
		star_genome_dir_tar_gz: {help: "The indexed reference genome files required for STAR."}
		run_pseudo_mapping_quantification: {help: "Option to map and directly quantify raw reads with Salmon. This and/or '' must be set to true. [false]"}
		run_salmon_index_ref_genome: {help: "Option to create decoy sequences (from genome), concatenating transcriptome and genome, and index concatenated genome with Salmon. [false]"}
		salmon_genome_dir_tar_gz: {help: "The indexed concatenated transcriptome and genome files required for Salmon."}
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (cellranger and generating the initial adata object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team cohort intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team cohort analysis outputs in."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in."}
	}
}
