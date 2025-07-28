version 1.0

# Differential gene expression analysis based on the negative binomial distribution

workflow differential_gene_expression_analysis {
	input {
		String team_id
		Array[Array[String]] project_sample_ids
		
		File metadata_csv
		File condition_csv
		File gene_map_csv
		File gene_ids_and_names_json

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
			team_id = team_id,
			project_sample_ids = project_sample_ids,
			metadata_csv = metadata_csv,
			condition_csv = condition_csv,
			gene_map_csv = gene_map_csv,
			gene_ids_and_names_json = gene_ids_and_names_json,
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
		String team_id
		Array[Array[String]] project_sample_ids
		
		File metadata_csv
		File condition_csv
		File gene_map_csv
		File gene_ids_and_names_json

		String salmon_mode
		Array[File] salmon_quant_tar_gz

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 4
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil((size([metadata_csv, condition_csv, gene_map_csv], "GB") + size(flatten([salmon_quant_tar_gz]), "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		while read -r quant_tar_gz || [[ -n "${quant_tar_gz}" ]]; do
			tar -xzvf "${quant_tar_gz}"
		done < ~{write_lines(salmon_quant_tar_gz)}

		python3 /opt/scripts/dge_analysis.py \
			--team-id ~{team_id} \
			--sample-ids ~{write_tsv(project_sample_ids)} \
			--condition-dict ~{condition_csv} \
			--metadata ~{metadata_csv} \
			--gene-map ~{gene_map_csv} \
			--gene-ids-and-names ~{gene_ids_and_names_json} \
			--salmon-mode ~{salmon_mode}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{team_id}.~{salmon_mode}.dds.pkl" \
			-o "~{team_id}.~{salmon_mode}.pydeseq2_significant_genes.csv" \
			-o "~{team_id}.~{salmon_mode}.volcano_plot.png"
	>>>

	output {
		String dds_object_pkl = "~{raw_data_path}/~{team_id}.~{salmon_mode}.dds.pkl"
		String significant_genes_csv = "~{raw_data_path}/~{team_id}.~{salmon_mode}.pydeseq2_significant_genes.csv"
		String volcano_plot_png = "~{raw_data_path}/~{team_id}.~{salmon_mode}.volcano_plot.png"
	}
	runtime {
		docker: "~{container_registry}/pydeseq2:0.5.2"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		maxRetries: 2
		bootDiskSizeGb: 30
		zones: zones
	}
}
