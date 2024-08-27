version 1.0

# Align reads

import "../../pmdbs_bulk_rnaseq_analysis_structs.wdl"

workflow alignment {
	input {
		String sample_id

		Boolean run_index_ref_genome
		ReferenceData reference
		File? star_genome_dir_tar_gz

		Array[File] trimmed_fastq_R1s
		Array[File] trimmed_fastq_R2s

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	if (run_index_ref_genome) {
		call index_ref_genome {
		input:
			primary_assembly_fasta = reference.primary_assembly_fasta,
			gene_annotation_gtf = reference.gene_annotation_gtf,
			container_registry = container_registry,
			zones = zones
		}
	}

	File genome_dir_tar_gz = select_first([
		index_ref_genome.star_genome_dir_tar_gz,
		star_genome_dir_tar_gz
	])

	call align {
		input:
			sample_id = sample_id,
			star_genome_dir_tar_gz = genome_dir_tar_gz,
			trimmed_fastq_R1s = trimmed_fastq_R1s,
			trimmed_fastq_R2s = trimmed_fastq_R2s,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		File aligned_bam = align.aligned_bam #!FileCoercion
		File aligned_bam_index = align.aligned_bam_index #!FileCoercion
		File unmapped_mate1 = align.unmapped_mate1 #!FileCoercion
		File unmapped_mate2 = align.unmapped_mate2 #!FileCoercion
		File log = align.log #!FileCoercion
		File final_log = align.final_log #!FileCoercion
		File progress_log = align.progress_log #!FileCoercion
		File sj_out_tab = align.sj_out_tab #!FileCoercion
	}
}

task index_ref_genome {
	input {
		File primary_assembly_fasta
		File gene_annotation_gtf
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([primary_assembly_fasta, gene_annotation_gtf], "GB") *2 + 100)

	command <<<
		set -euo pipefail

		mkdir -p star_genome_dir

		STAR \
			--runThreadN ~{threads - 1} \
			--runMode genomeGenerate \
			--genomeDir star_genome_dir \
			--genomeFastaFiles ~{primary_assembly_fasta} \
			--sjdbGTFfile ~{gene_annotation_gtf}

		tar -czvf star_genome_dir.tar.gz star_genome_dir
	>>>

	output {
		File star_genome_dir_tar_gz = "star_genome_dir.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/star_samtools:2.7.11b_1.20"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}

task align {
	input {
		String sample_id

		File star_genome_dir_tar_gz

		Array[File] trimmed_fastq_R1s
		Array[File] trimmed_fastq_R2s

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(flatten([trimmed_fastq_R1s, trimmed_fastq_R2s]), "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		tar -xzvf ~{star_genome_dir_tar_gz}

		STAR \
			--runThreadN ~{threads - 1} \
			--genomeDir star_genome_dir \
			--readFilesIn ~{sep=',' trimmed_fastq_R1s} ~{sep=',' trimmed_fastq_R2s} \
			--readFilesCommand zcat \
			--outFileNamePrefix ~{sample_id}_ \
			--outReadsUnmapped Fastx \
			--outSAMtype BAM SortedByCoordinate \
			--outFilterType BySJout \
			--alignIntronMax 1000000 \
			--alignMatesGapMax 1000000 \
			--twopassMode Basic

		samtools index \
			-@ ~{threads} \
			~{sample_id}_Aligned.sortedByCoord.out.bam

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}_Aligned.sortedByCoord.out.bam" \
			-o "~{sample_id}_Aligned.sortedByCoord.out.bam.bai" \
			-o "~{sample_id}_Unmapped.out.mate1" \
			-o "~{sample_id}_Unmapped.out.mate2" \
			-o "~{sample_id}_Log.out" \
			-o "~{sample_id}_Log.final.out" \
			-o "~{sample_id}_Log.progress.out" \
			-o "~{sample_id}_SJ.out.tab"
	>>>

	output {
		String aligned_bam = "~{raw_data_path}/~{sample_id}_Aligned.sortedByCoord.out.bam"
		String aligned_bam_index = "~{raw_data_path}/~{sample_id}_Aligned.sortedByCoord.out.bam.bai"
		String unmapped_mate1 = "~{raw_data_path}/~{sample_id}_Unmapped.out.mate1"
		String unmapped_mate2 = "~{raw_data_path}/~{sample_id}_Unmapped.out.mate2"
		String log = "~{raw_data_path}/~{sample_id}_Log.out"
		String final_log = "~{raw_data_path}/~{sample_id}_Log.final.out"
		String progress_log = "~{raw_data_path}/~{sample_id}_Log.progress.out"
		String sj_out_tab = "~{raw_data_path}/~{sample_id}_SJ.out.tab"
	}

	runtime {
		docker: "~{container_registry}/star_samtools:2.7.11b_1.20"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
