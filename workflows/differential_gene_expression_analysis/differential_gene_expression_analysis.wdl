version 1.0

# Differential gene expression analysis based on the negative binomial distribution

import "../pmdbs_bulk_rnaseq_analysis_structs.wdl"

workflow differential_gene_expression_analysis {
	input {
		String cohort_id
		Array[String] sample_ids

		File metadata_csv
		File gene_map_csv
		File blacklist_genes_bed

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

	# Task and subworkflow versions
	String sub_workflow_name = "differential_gene_expression_analysis"
	String differential_gene_expression_task_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String pydeseq2_raw_data_path = "~{workflow_raw_data_path_prefix}/differential_gene_expression/~{differential_gene_expression_task_version}"

	# TODO add step to modify metadata in task

	call differential_gene_expression {
		input:
			cohort_id = cohort_id,
			sample_ids = sample_ids,
			metadata_csv = metadata_csv,
			gene_map_csv = gene_map_csv,
			blacklist_genes_bed = blacklist_genes_bed,
			salmon_mode = salmon_mode,
			salmon_quant_tar_gz = salmon_quant_tar_gz,
			raw_data_path = pydeseq2_raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		# PyDESeq2 DGE Analysis
		File significant_genes_csv = differential_gene_expression.significant_genes_csv #!FileCoercion
		File pca_plot_png = differential_gene_expression.pca_plot_png #!FileCoercion
		File volcano_plot_png = differential_gene_expression.volcano_plot_png #!FileCoercion
	}
}

task differential_gene_expression {
	input {
		String cohort_id
		Array[String] sample_ids

		File metadata_csv
		File gene_map_csv
		File blacklist_genes_bed

		String salmon_mode
		Array[File] salmon_quant_tar_gz

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil((size([metadata_csv, gene_map_csv, blacklist_genes_bed], "GB") + size(flatten([salmon_quant_tar_gz]), "GB")) * 2 + 50)

	command <<<
		set -euo pipefail

		while read -r quant_tar_gz || [[ -n "${quant_tar_gz}" ]]; do
			tar -xzvf "${quant_tar_gz}"
		done < ~{write_lines(salmon_quant_tar_gz)}

		python3 pydeseq2.py \
			--cohort-id ~{cohort_id} \
			--samples ~{sep = '' sample_ids} \
			--metadata ~{metadata_csv} \
			--gene-map ~{gene_map_csv} \
			--blacklist-genes ~{blacklist_genes_bed} \
			--salmon-mode ~{salmon_mode}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.~{salmon_mode}.pydeseq2_significant_genes.csv" \
			-o "~{cohort_id}.~{salmon_mode}.pca_plot.png" \
			-o "~{cohort_id}.~{salmon_mode}.volcano_plot.png"
	>>>

	output {
		String significant_genes_csv = "~{raw_data_path}/~{cohort_id}.~{salmon_mode}.pydeseq2_significant_genes.csv"
		String pca_plot_png = "~{raw_data_path}/~{cohort_id}.~{salmon_mode}.pca_plot.png"
		String volcano_plot_png = "~{raw_data_path}/~{cohort_id}.~{salmon_mode}.volcano_plot.png"
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
