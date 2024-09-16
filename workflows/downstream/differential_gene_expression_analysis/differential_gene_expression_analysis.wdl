version 1.0

# Differential gene expression analysis based on the negative binomial distribution

workflow differential_gene_expression_analysis {
	input {
		String project_id
		
		File metadata_csv # TODO - ASAP_sample_id will be metadata once DTi processes them (SAMPLE.csv)
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

	call differential_gene_expression {
		input:
			project_id = project_id,
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

	output {
		# PyDESeq2 DGE Analysis
		File dds_object_pkl = differential_gene_expression.dds_object_pkl #!FileCoercion
		File significant_genes_csv = differential_gene_expression.significant_genes_csv #!FileCoercion
		File volcano_plot_png = differential_gene_expression.volcano_plot_png #!FileCoercion
	}
}

task differential_gene_expression {
	input {
		String project_id
		
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
			--cohort-id ~{project_id} \
			--metadata ~{metadata_csv} \
			--gene-map ~{gene_map_csv} \
			--blacklist-genes ~{blacklist_genes_bed} \
			--salmon-mode ~{salmon_mode}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{project_id}.~{salmon_mode}.dds.pkl" \
			-o "~{project_id}.~{salmon_mode}.pydeseq2_significant_genes.csv" \
			-o "~{project_id}.~{salmon_mode}.volcano_plot.png"
	>>>

	output {
		String dds_object_pkl = "~{raw_data_path}/~{project_id}.~{salmon_mode}.dds.pkl"
		String significant_genes_csv = "~{raw_data_path}/~{project_id}.~{salmon_mode}.pydeseq2_significant_genes.csv"
		String volcano_plot_png = "~{raw_data_path}/~{project_id}.~{salmon_mode}.volcano_plot.png"
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
