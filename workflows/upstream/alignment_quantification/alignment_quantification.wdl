version 1.0

# Align reads and quantify (alignment-based mode)

workflow alignment_quantification {
	input {
		String dataset_sample_id
		String sample_id

		File all_transcripts_fasta
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
			dataset_sample_id = dataset_sample_id,
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
			dataset_sample_id = dataset_sample_id,
			sample_id = sample_id,
			all_transcripts_fasta = all_transcripts_fasta,
			aligned_to_transcriptome_bam = alignment.aligned_to_transcriptome_bam, #!FileCoercion
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
		File aligned_to_transcriptome_bam = alignment.aligned_to_transcriptome_bam #!FileCoercion
		File unmapped_mate1 = alignment.unmapped_mate1 #!FileCoercion
		File unmapped_mate2 = alignment.unmapped_mate2 #!FileCoercion
		File log = alignment.log #!FileCoercion
		File final_log = alignment.final_log #!FileCoercion
		File progress_log = alignment.progress_log #!FileCoercion
		File sj_out_tab = alignment.sj_out_tab #!FileCoercion

		# Salmon quantification
		File quant_tar_gz = quantification.quant_tar_gz #!FileCoercion
	}

	meta {
		description: "Aligns trimmed reads to the reference genome with STAR and quantifies aligned reads with Salmon in alignment-based mode."
	}

	parameter_meta {
		dataset_sample_id: {help: "Generated ASAP dataset ID and sample ID; used to name output files."}
		sample_id: {help: "Generated ASAP sample ID; used to name output files."}
		all_transcripts_fasta: {help: "Manually generated all transcripts on the reference chromosomes with the `primary_assembly_fasta` and `gene_annotation_gtf`."}
		star_genome_dir_tar_gz: {help: "The indexed reference genome files required for STAR."}
		trimmed_fastq_R1s: {help: "Adapter-trimmed forward (R1) FASTQ files for the sample."}
    	trimmed_fastq_R2s: {help: "Adapter-trimmed reverse (R2) FASTQ files for the sample."}
		raw_data_path: {help: "Raw data bucket path for alignment and quantification outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/upstream/alignment_quantification`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task alignment {
	input {
		String dataset_sample_id

		File star_genome_dir_tar_gz

		Array[File] trimmed_fastq_R1s
		Array[File] trimmed_fastq_R2s

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 48
	Int mem_gb = ceil(threads * 2)
	Int sort_bam_mem_bytes = (mem_gb - 20) * 1024 * 1024 * 1024
	Int disk_size = ceil((size(star_genome_dir_tar_gz, "GB") + size(flatten([trimmed_fastq_R1s, trimmed_fastq_R2s]), "GB")) * 5 + 500)

	command <<<
		set -euo pipefail

		tar -xzvf ~{star_genome_dir_tar_gz}

		/usr/bin/time \
		STAR \
			--runThreadN ~{threads - 1} \
			--genomeDir star_genome_dir \
			--readFilesIn ~{sep=',' trimmed_fastq_R1s} ~{sep=',' trimmed_fastq_R2s} \
			--readFilesCommand zcat \
			--outFileNamePrefix ~{dataset_sample_id}. \
			--outReadsUnmapped Fastx \
			--outSAMtype BAM SortedByCoordinate \
			--outFilterType BySJout \
			--alignIntronMax 1000000 \
			--alignMatesGapMax 1000000 \
			--twopassMode Basic \
			--quantMode TranscriptomeSAM \
			--limitBAMsortRAM ~{sort_bam_mem_bytes}

		echo "Validating aligned and sorted BAM"
		samtools quickcheck "~{dataset_sample_id}.Aligned.sortedByCoord.out.bam"

		echo "Indexing aligned and sorted BAM"
		samtools index \
			-@ ~{threads} \
			~{dataset_sample_id}.Aligned.sortedByCoord.out.bam

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{dataset_sample_id}.Aligned.sortedByCoord.out.bam" \
			-o "~{dataset_sample_id}.Aligned.sortedByCoord.out.bam.bai" \
			-o "~{dataset_sample_id}.Aligned.toTranscriptome.out.bam" \
			-o "~{dataset_sample_id}.Unmapped.out.mate1" \
			-o "~{dataset_sample_id}.Unmapped.out.mate2" \
			-o "~{dataset_sample_id}.Log.out" \
			-o "~{dataset_sample_id}.Log.final.out" \
			-o "~{dataset_sample_id}.Log.progress.out" \
			-o "~{dataset_sample_id}.SJ.out.tab"
	>>>

	output {
		String aligned_bam = "~{raw_data_path}/~{dataset_sample_id}.Aligned.sortedByCoord.out.bam"
		String aligned_bam_index = "~{raw_data_path}/~{dataset_sample_id}.Aligned.sortedByCoord.out.bam.bai"
		String aligned_to_transcriptome_bam = "~{raw_data_path}/~{dataset_sample_id}.Aligned.toTranscriptome.out.bam"
		String unmapped_mate1 = "~{raw_data_path}/~{dataset_sample_id}.Unmapped.out.mate1"
		String unmapped_mate2 = "~{raw_data_path}/~{dataset_sample_id}.Unmapped.out.mate2"
		String log = "~{raw_data_path}/~{dataset_sample_id}.Log.out"
		String final_log = "~{raw_data_path}/~{dataset_sample_id}.Log.final.out"
		String progress_log = "~{raw_data_path}/~{dataset_sample_id}.Log.progress.out"
		String sj_out_tab = "~{raw_data_path}/~{dataset_sample_id}.SJ.out.tab"
	}

	runtime {
		docker: "~{container_registry}/star_samtools:2.7.11b_1.20"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} SSD"
		preemptible: 3
		zones: zones
	}

	meta {
		description: "Aligns trimmed paired-end reads to the reference genome using STAR two-pass mode."
	}

	parameter_meta {
		dataset_sample_id: {help: "Generated ASAP dataset ID and sample ID; used to name output files."}
		star_genome_dir_tar_gz: {help: "The indexed reference genome files required for STAR."}
		trimmed_fastq_R1s: {help: "Adapter-trimmed forward (R1) FASTQ files for the sample."}
    	trimmed_fastq_R2s: {help: "Adapter-trimmed reverse (R2) FASTQ files for the sample."}
		raw_data_path: {help: "Raw data bucket path for alignment outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/upstream/alignment_quantification/<alignment_quantification_workflow_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task quantification {
	input {
		String dataset_sample_id
		String sample_id

		File all_transcripts_fasta

		File aligned_to_transcriptome_bam

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([all_transcripts_fasta, aligned_to_transcriptome_bam], "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		salmon quant \
			--targets ~{all_transcripts_fasta} \
			--libType A \
			--alignments ~{aligned_to_transcriptome_bam} \
			--output ~{sample_id}_salmon_quant \
			--threads ~{threads} \
			--gcBias

		# Outputs must remain in folder and unmodified for downstream analysis
		# Outputs include: quant.sf, cmd_info.json, and aux_info folder
		tar -czvf "~{dataset_sample_id}.alignment_mode.salmon_quant.tar.gz" "~{sample_id}_salmon_quant"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{dataset_sample_id}.alignment_mode.salmon_quant.tar.gz"
	>>>

	output {
		String quant_tar_gz = "~{raw_data_path}/~{dataset_sample_id}.alignment_mode.salmon_quant.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/salmon:1.10.3"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}

	meta {
		description: "Quantifies transcript abundances from a STAR transcriptome-aligned BAM using Salmon in alignment-based mode."
	}

	parameter_meta {
		dataset_sample_id: {help: "Generated ASAP dataset ID and sample ID; used to name output files."}
		sample_id: {help: "Generated ASAP sample ID; used to name output files."}
		all_transcripts_fasta: {help: "Manually generated all transcripts on the reference chromosomes with the `primary_assembly_fasta` and `gene_annotation_gtf`."}
		aligned_to_transcriptome_bam: {help: "BAM file aligned to the transcriptome, output from STAR with --quantMode TranscriptomeSAM."}
		raw_data_path: {help: "Raw data bucket path for quantification outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/upstream/alignment_quantification/<alignment_quantification_workflow_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
