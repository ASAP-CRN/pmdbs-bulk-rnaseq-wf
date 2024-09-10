version 1.0

# Harmonized human PMDBS bulk RNAseq workflow entrypoint

import "pmdbs_bulk_rnaseq_analysis_structs.wdl"
import "wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "index_ref_genome/index_ref_genome.wdl" as IndexRefGenome
import "upstream/upstream.wdl" as Upstream
import "downstream/downstream.wdl" as Downstream

workflow pmdbs_bulk_rnaseq_analysis {
	input {
		String cohort_id
		Array[Project] projects

		ReferenceData reference

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

		Array[String] upstream_output_file_paths = select_all(
			flatten([
				upstream.fastqc_reports_tar_gz,
				flatten(upstream.trimmed_fastq_R1s),
				flatten(upstream.trimmed_fastq_R2s),
				upstream.failed_paired_fastqs_tar_gz,
				upstream.reports_html_tar_gz,
				upstream.trimmed_fastqc_reports_tar_gz,
				upstream.aligned_bam,
				upstream.aligned_bam_index,
				upstream.unmapped_mate1,
				upstream.unmapped_mate2,
				upstream.log,
				upstream.final_log,
				upstream.progress_log,
				upstream.sj_out_tab,
				upstream.alignment_mode_quant_tar_gz,
				upstream.mapping_mode_quant_tar_gz,
				[
					upstream.alignment_mode_multiqc_report_html,
					upstream.alignment_mode_multiqc_data_tar_gz,
					upstream.mapping_mode_multiqc_report_html,
					upstream.mapping_mode_multiqc_data_tar_gz
				]
			])
		)
		#!StringCoercion

		if (run_alignment_quantification) {
			if (project.run_project_cohort_analysis) {
				call Downstream.downstream as alignment_mode_project_downstream {
					input:
						cohort_id = project.project_id,
						project_sample_ids = upstream.project_sample_ids,
						upstream_output_file_paths = upstream_output_file_paths,
						metadata_csv = metadata_csv,
						gene_map_csv = gene_map_csv,
						blacklist_genes_bed = blacklist_genes_bed,
						salmon_mode = "alignment_mode",
						salmon_quant_tar_gz = select_all(upstream.alignment_mode_quant_tar_gz),
						workflow_name = workflow_name,
						workflow_version = workflow_version,
						workflow_release = workflow_release,
						run_timestamp = get_workflow_metadata.timestamp,
						raw_data_path_prefix = project_raw_data_path_prefix,
						staging_data_buckets = project.staging_data_buckets,
						billing_project = get_workflow_metadata.billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}
		}

		if (run_pseudo_mapping_quantification) {
			if (project.run_project_cohort_analysis) {
				call Downstream.downstream as mapping_mode_project_downstream {
					input:
						cohort_id = project.project_id,
						project_sample_ids = upstream.project_sample_ids,
						upstream_output_file_paths = upstream_output_file_paths,
						metadata_csv = metadata_csv,
						gene_map_csv = gene_map_csv,
						blacklist_genes_bed = blacklist_genes_bed,
						salmon_mode = "mapping_mode",
						salmon_quant_tar_gz = select_all(upstream.mapping_mode_quant_tar_gz),
						workflow_name = workflow_name,
						workflow_version = workflow_version,
						workflow_release = workflow_release,
						run_timestamp = get_workflow_metadata.timestamp,
						raw_data_path_prefix = project_raw_data_path_prefix,
						staging_data_buckets = project.staging_data_buckets,
						billing_project = get_workflow_metadata.billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}
		}
	}

	if (run_cross_team_cohort_analysis) {
		String cohort_raw_data_path_prefix = "~{cohort_raw_data_bucket}/~{workflow_execution_path}/~{workflow_name}"

		if (run_alignment_quantification) {
			call Downstream.downstream as alignment_mode_cross_team_downstream {
				input:
					cohort_id = cohort_id,
					project_sample_ids = flatten(upstream.project_sample_ids),
					upstream_output_file_paths = upstream_output_file_paths,
					metadata_csv = metadata_csv,
					gene_map_csv = gene_map_csv,
					blacklist_genes_bed = blacklist_genes_bed,
					salmon_mode = "alignment_mode",
					salmon_quant_tar_gz = select_all(flatten(upstream.alignment_mode_quant_tar_gz)),
					workflow_name = workflow_name,
					workflow_version = workflow_version,
					workflow_release = workflow_release,
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = cohort_raw_data_path_prefix,
					staging_data_buckets = cohort_staging_data_buckets,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		if (run_pseudo_mapping_quantification) {
			call Downstream.downstream  as mapping_mode_cross_team_downstream {
				input:
					cohort_id = cohort_id,
					project_sample_ids = flatten(upstream.project_sample_ids),
					upstream_output_file_paths = upstream_output_file_paths,
					metadata_csv = metadata_csv,
					gene_map_csv = gene_map_csv,
					blacklist_genes_bed = blacklist_genes_bed,
					salmon_mode = "mapping_mode",
					salmon_quant_tar_gz = select_all(flatten(upstream.mapping_mode_quant_tar_gz)),
					workflow_name = workflow_name,
					workflow_version = workflow_version,
					workflow_release = workflow_release,
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = cohort_raw_data_path_prefix,
					staging_data_buckets = cohort_staging_data_buckets,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}
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

		# Project downstream analysis outputs
		## List of samples included in the cohort. Both modes produce the same sample list
		Array[File?] project_cohort_sample_list = alignment_mode_project_downstream.cohort_sample_list

		# DGE analysis with Salmon mapping-mode counts
		Array[File?] project_pydeseq2_alignment_mode_significant_genes_csv = alignment_mode_project_downstream.significant_genes_csv
		Array[File?] project_pydeseq2_alignment_mode_pca_plot_png = alignment_mode_project_downstream.pca_plot_png
		Array[File?] project_pydeseq2_alignment_mode_volcano_plot_png = alignment_mode_project_downstream.volcano_plot_png

		Array[Array[File]?] alignment_mode_upstream_manifests = alignment_mode_project_downstream.upstream_manifest_tsvs
		Array[Array[File]?] alignment_mode_project_manifests = alignment_mode_project_downstream.downstream_manifest_tsvs

		# DGE analysis with Salmon mapping-mode counts
		Array[File?] project_pydeseq2_mapping_mode_significant_genes_csv = mapping_mode_project_downstream.significant_genes_csv
		Array[File?] project_pydeseq2_mapping_mode_pca_plot_png = mapping_mode_project_downstream.pca_plot_png
		Array[File?] project_pydeseq2_mapping_mode_volcano_plot_png = mapping_mode_project_downstream.volcano_plot_png

		# TODO this will replace the alignment_mode_project_downstream MANIFEST.tsv's
		Array[Array[File]?] mapping_mode_upstream_manifests = mapping_mode_project_downstream.upstream_manifest_tsvs
		Array[Array[File]?] mapping_mode_project_manifests = mapping_mode_project_downstream.downstream_manifest_tsvs

		# Cross-team downstream outputs
		## List of samples included in the cohort. Both modes produce the same sample list
		File? cohort_sample_list = alignment_mode_cross_team_downstream.cohort_sample_list

		# DGE analysis with Salmon mapping-mode counts
		File? cohort_pydeseq2_alignment_mode_significant_genes_csv = alignment_mode_cross_team_downstream.significant_genes_csv
		File? cohort_pydeseq2_alignment_mode_pca_plot_png = alignment_mode_cross_team_downstream.pca_plot_png
		File? cohort_pydeseq2_alignment_mode_volcano_plot_png = alignment_mode_cross_team_downstream.volcano_plot_png

		Array[File]? alignment_mode_cohort_manifests = alignment_mode_cross_team_downstream.downstream_manifest_tsvs

		# DGE analysis with Salmon mapping-mode counts
		File? cohort_pydeseq2_mapping_mode_significant_genes_csv = mapping_mode_cross_team_downstream.significant_genes_csv
		File? cohort_pydeseq2_mapping_mode_pca_plot_png = mapping_mode_cross_team_downstream.pca_plot_png
		File? cohort_pydeseq2_mapping_mode_volcano_plot_png = mapping_mode_cross_team_downstream.volcano_plot_png

		Array[File]?  mapping_mode_cohort_manifests = mapping_mode_cross_team_downstream.downstream_manifest_tsvs
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
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only upstream steps (QC, align/map, and quantify) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team downstream intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team downstream outputs in."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in."}
	}
}
