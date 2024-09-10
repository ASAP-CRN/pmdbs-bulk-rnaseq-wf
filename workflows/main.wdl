version 1.0

# Harmonized human PMDBS bulk RNAseq workflow entrypoint

import "pmdbs_bulk_rnaseq_analysis_structs.wdl"
import "wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "index_ref_genome/index_ref_genome.wdl" as IndexRefGenome
import "upstream/upstream.wdl" as Upstream

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

		File metadata_csv
		File gene_map_csv
		File blacklist_genes_bed

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	String workflow_execution_path = "workflow_execution"
	String workflow_name = "pmdbs_bulk_rnaseq_analysis"
	String workflow_version = "v1.0.0"
	String workflow_release = "https://github.com/ASAP-CRN/pmdbs-bulk-rnaseq-wf/releases/tag/pmdbs_bulk_rnaseq_analysis-~{workflow_version}"

	call GetWorkflowMetadata.get_workflow_metadata {
		input:
			zones = zones
	}

	call IndexRefGenome.index_ref_genome {
		input:
			reference = reference,
			run_star_index_ref_genome = run_star_index_ref_genome,
			run_salmon_index_ref_genome = run_salmon_index_ref_genome,
			container_registry = container_registry,
			zones = zones
	}

	scatter (project in projects) {
		# Add workflow name or else intermediate files will be uploaded to the same folder as harmonized_pmdbs
		String project_raw_data_path_prefix = "~{project.raw_data_bucket}/~{workflow_execution_path}/~{workflow_name}"

		call Upstream.upstream {
			input:
				project_id = project.project_id,
				samples = project.samples,
				transcripts_fasta = reference.transcripts_fasta,
				run_alignment_quantification = run_alignment_quantification,
				star_genome_dir_tar_gz = select_first([index_ref_genome.star_genome_dir_tar_gz, star_genome_dir_tar_gz]),
				run_pseudo_mapping_quantification = run_pseudo_mapping_quantification,
				salmon_genome_dir_tar_gz = select_first([index_ref_genome.salmon_genome_dir_tar_gz, salmon_genome_dir_tar_gz]),
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		Array[String] upstream_output_file_paths = flatten([
			upstream.fastqc_reports_tar_gz,
			flatten(upstream.trimmed_fastq_R1s),
			flatten(upstream.trimmed_fastq_R2s),
			upstream.failed_paired_fastqs_tar_gz,
			upstream.reports_html_tar_gz,
			upstream.trimmed_fastqc_reports_tar_gz,
			select_all([
				upstream.aligned_bam,
				upstream.aligned_bam_index,
				upstream.unmapped_mate1,
				upstream.unmapped_mate2,
				upstream.log,
				upstream.final_log,
				upstream.progress_log,
				upstream.sj_out_tab,
				upstream.alignment_mode_quant_tar_gz,
				upstream.alignment_mode_multiqc_report_html,
				upstream.alignment_mode_multiqc_data_tar_gz,
				upstream.mapping_mode_quant_tar_gz,
				upstream.mapping_mode_multiqc_report_html,
				upstream.mapping_mode_multiqc_data_tar_gz
			])
		]) #!StringCoercion
	}

	output {
		# Sample-level outputs
		# Sample list
		Array[Array[Array[String]]] project_sample_ids = upstream.project_sample_ids

		# Preprocess
		Array[Array[File]] fastqc_raw_reads_reports_tar_gz = upstream.fastqc_reports_tar_gz
		Array[Array[Array[File]]] fastp_trimmed_fastq_R1s = upstream.trimmed_fastq_R1s
		Array[Array[Array[File]]] fastp_trimmed_fastq_R2s = upstream.trimmed_fastq_R2s
		Array[Array[File]] fastp_failed_paired_fastqs_tar_gz = upstream.failed_paired_fastqs_tar_gz
		Array[Array[File]] fastp_report_html_tar_gz = upstream.reports_html_tar_gz
		Array[Array[File]] fastqc_trimmed_reads_reports_tar_gz = upstream.trimmed_fastqc_reports_tar_gz

		# Alignment and quantification
		Array[Array[File?]] star_aligned_bam = upstream.aligned_bam
		Array[Array[File?]] star_aligned_bam_index = upstream.aligned_bam_index
		Array[Array[File?]] star_unmapped_mate1 = upstream.unmapped_mate1
		Array[Array[File?]] star_unmapped_mate2 = upstream.unmapped_mate2
		Array[Array[File?]] star_log = upstream.log
		Array[Array[File?]] star_final_log = upstream.final_log
		Array[Array[File?]] star_progress_log = upstream.progress_log
		Array[Array[File?]] star_sj_out_tab = upstream.sj_out_tab
		Array[Array[File?]] salmon_alignment_mode_quant_tar_gz = upstream.alignment_mode_quant_tar_gz
		Array[File?] salmon_alignment_mode_multiqc_report_html = upstream.alignment_mode_multiqc_report_html
		Array[File?] salmon_alignment_mode_multiqc_data_tar_gz = upstream.alignment_mode_multiqc_data_tar_gz

		# Direct quantification with mapping
		Array[Array[File?]] salmon_mapping_mode_quant_tar_gz = upstream.mapping_mode_quant_tar_gz
		Array[File?] salmon_mapping_mode_multiqc_report_html = upstream.mapping_mode_multiqc_report_html
		Array[File?] salmon_mapping_mode_multiqc_data_tar_gz = upstream.mapping_mode_multiqc_data_tar_gz
	}

	meta {
		description: "Harmonized human postmortem-derived brain sequencing (PMDBS) bulk RNA-seq workflow"
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files during cross-team cohort analysis."}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level cohort analysis."}
		reference: {help: "The primary assembly FASTA and gene annotation GTF from GENCODE."}
		run_alignment_quantification: {help: "Option to align raw reads with STAR and quantify aligned reads with Salmon. This and/or 'run_pseudo_mapping_quantification' must be set to true. [true]"}
		run_star_index_ref_genome: {help: "Option to index reference genome with STAR. If set to false, star_genome_dir_tar_gz must be provided. [false]"}
		star_genome_dir_tar_gz: {help: "The indexed reference genome files required for STAR."}
		run_pseudo_mapping_quantification: {help: "Option to map and directly quantify raw reads with Salmon. This and/or 'run_alignment_quantification' must be set to true. [false]"}
		run_salmon_index_ref_genome: {help: "Option to create decoy sequences (from genome), concatenating transcriptome and genome, and index concatenated genome with Salmon. If set to false, salmon_genome_dir_tar_gz must be provided [false]"}
		salmon_genome_dir_tar_gz: {help: "The indexed concatenated transcriptome and genome files required for Salmon."}
		metadata_csv: {help: "CSV containing all sample information including batch, condition, etc."}
		gene_map_csv: {help: "CSV containing mapped transcript IDs and gene IDs that must be in this order."}
		blacklist_genes_bed: {help: "BED file containing the ENCODE Blacklist genes."}
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (cellranger and generating the initial adata object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team cohort intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team cohort analysis outputs in."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in."}
	}
}
