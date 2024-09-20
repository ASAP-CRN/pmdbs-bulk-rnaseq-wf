version 1.0

# Align reads and quantify (alignment-based mode)

workflow alignment_quantification {
	input {
		String sample_id

		File transcripts_fasta
		File star_genome_dir_tar_gz

		Array[File] trimmed_fastq_R1s
		Array[File] trimmed_fastq_R2s

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	call alignment {
		input:
			sample_id = sample_id,
			star_genome_dir_tar_gz = star_genome_dir_tar_gz,
			trimmed_fastq_R1s = trimmed_fastq_R1s,
			trimmed_fastq_R2s = trimmed_fastq_R2s,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}
		
	call quantification {
		input:
			sample_id = sample_id,
			transcripts_fasta = transcripts_fasta,
			aligned_bam = alignment.aligned_bam, #!FileCoercion
			aligned_bam_index = alignment.aligned_bam_index, #!FileCoercion
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		# STAR alignment
		File aligned_bam = alignment.aligned_bam #!FileCoercion
		File aligned_bam_index = alignment.aligned_bam_index #!FileCoercion
		File unmapped_mate1 = alignment.unmapped_mate1 #!FileCoercion
		File unmapped_mate2 = alignment.unmapped_mate2 #!FileCoercion
		File log = alignment.log #!FileCoercion
		File final_log = alignment.final_log #!FileCoercion
		File progress_log = alignment.progress_log #!FileCoercion
		File sj_out_tab = alignment.sj_out_tab #!FileCoercion

		# Salmon quantification
		File quant_tar_gz = quantification.quant_tar_gz #!FileCoercion
	}
}

task alignment {
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

	Int threads = 56
	Int mem_gb = ceil(threads * 3)
	Int disk_size = ceil((size(star_genome_dir_tar_gz, "GB") + size(flatten([trimmed_fastq_R1s, trimmed_fastq_R2s]), "GB")) * 5 + 200)

	command <<<
		set -euo pipefail

		tar -xzvf ~{star_genome_dir_tar_gz}

		STAR \
			--runThreadN ~{threads - 1} \
			--genomeDir star_genome_dir \
			--readFilesIn ~{sep=',' trimmed_fastq_R1s} ~{sep=',' trimmed_fastq_R2s} \
			--readFilesCommand zcat \
			--outFileNamePrefix ~{sample_id}. \
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
			-o "~{sample_id}.Aligned.sortedByCoord.out.bam" \
			-o "~{sample_id}.Aligned.sortedByCoord.out.bam.bai" \
			-o "~{sample_id}.Unmapped.out.mate1" \
			-o "~{sample_id}.Unmapped.out.mate2" \
			-o "~{sample_id}.Log.out" \
			-o "~{sample_id}.Log.final.out" \
			-o "~{sample_id}.Log.progress.out" \
			-o "~{sample_id}.SJ.out.tab"
	>>>

	output {
		String aligned_bam = "~{raw_data_path}/~{sample_id}.Aligned.sortedByCoord.out.bam"
		String aligned_bam_index = "~{raw_data_path}/~{sample_id}.Aligned.sortedByCoord.out.bam.bai"
		String unmapped_mate1 = "~{raw_data_path}/~{sample_id}.Unmapped.out.mate1"
		String unmapped_mate2 = "~{raw_data_path}/~{sample_id}.Unmapped.out.mate2"
		String log = "~{raw_data_path}/~{sample_id}.Log.out"
		String final_log = "~{raw_data_path}/~{sample_id}.Log.final.out"
		String progress_log = "~{raw_data_path}/~{sample_id}.Log.progress.out"
		String sj_out_tab = "~{raw_data_path}/~{sample_id}.SJ.out.tab"
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

task quantification {
	input {
		String sample_id

		File transcripts_fasta

		File aligned_bam
		File aligned_bam_index

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([transcripts_fasta, aligned_bam, aligned_bam_index], "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		salmon quant \
			--transcripts ~{transcripts_fasta} \
			--libType A \
			--alignments ~{aligned_bam} \
			--output ~{sample_id}_salmon_quant \
			--validateMappings \
			--threads ~{threads} \
			--rangeFactorizationBins 4 \
			--useErrorModel \
			--gcBias

		# Outputs must remain in folder and unmodified for downstream analysis
		# Outputs include: quant.sf, cmd_info.json, and aux_info folder
		tar -czvf "~{sample_id}.alignment_mode.salmon_quant.tar.gz" "~{sample_id}_salmon_quant"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.alignment_mode.salmon_quant.tar.gz"
	>>>

	output {
		String quant_tar_gz = "~{raw_data_path}/~{sample_id}.alignment_mode.salmon_quant.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/salmon:1.10.3"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
