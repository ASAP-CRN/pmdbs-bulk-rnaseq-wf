version 1.0

# Align reads and quantify (alignment-based mode)

import "../pmdbs_bulk_rnaseq_analysis_structs.wdl"

workflow alignment_quantification {
	input {
		Array[Sample] samples

		Boolean run_index_ref_genome
		ReferenceData reference
		File star_genome_dir_tar_gz

		Array[File] trimmed_fastq_R1s
		Array[File] trimmed_fastq_R2s

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
	String sub_workflow_name = "alignment_quantification"
	String alignment_task_version = "1.0.0"
	String quantification_task_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String star_raw_data_path = "~{workflow_raw_data_path_prefix}/alignment/~{alignment_task_version}"
	String salmon_raw_data_path = "~{workflow_raw_data_path_prefix}/quantification/~{quantification_task_version}"

	scatter (sample_object in samples) {
		String alignment_output = "~{star_raw_data_path}/~{sample_object.sample_id}.Aligned.sortedByCoord.out.bam"
		String quantification_output = "~{salmon_raw_data_path}/~{sample_object.sample_id}.alignment_mode.quant.sf"
	}

	# For each sample, outputs an array of true/false: [alignment_complete, quantification_complete]
	call check_output_files_exist {
		input:
			alignment_output_files = alignment_output,
			quantification_output_files = quantification_output,
			billing_project = billing_project,
			zones = zones
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

	scatter (sample_index in range(length(samples))) {
		Sample sample = samples[sample_index]

		String alignment_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][0]
		String quantification_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][1]

		String star_aligned_bam = "~{star_raw_data_path}/~{sample.sample_id}.Aligned.sortedByCoord.out.bam"
		String star_aligned_bam_index = "~{star_raw_data_path}/~{sample.sample_id}.Aligned.sortedByCoord.out.bam.bai"
		String star_unmapped_mate1 = "~{star_raw_data_path}/~{sample.sample_id}.Unmapped.out.mate1"
		String star_unmapped_mate2 = "~{star_raw_data_path}/~{sample.sample_id}.Unmapped.out.mate2"
		String star_log = "~{star_raw_data_path}/~{sample.sample_id}.Log.out"
		String star_final_log = "~{star_raw_data_path}/~{sample.sample_id}.Log.final.out"
		String star_progress_log = "~{star_raw_data_path}/~{sample.sample_id}.Log.progress.out"
		String star_sj_out_tab = "~{star_raw_data_path}/~{sample.sample_id}.SJ.out.tab"

		if (alignment_complete == "false") {
			call alignment {
				input:
					sample_id = sample.sample_id,
					star_genome_dir_tar_gz = genome_dir_tar_gz,
					trimmed_fastq_R1s = trimmed_fastq_R1s,
					trimmed_fastq_R2s = trimmed_fastq_R2s,
					raw_data_path = star_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		IndexData aligned_bam_output = {
			"data": select_first([alignment.aligned_bam, star_aligned_bam]), #!FileCoercion
			"data_index": select_first([alignment.aligned_bam_index, star_aligned_bam_index]) #!FileCoercion
		}
		File unmapped_mate1_output = select_first([alignment.unmapped_mate1, star_unmapped_mate1]) #!FileCoercion
		File unmapped_mate2_output = select_first([alignment.unmapped_mate2, star_unmapped_mate2]) #!FileCoercion
		File log_output = select_first([alignment.log, star_log]) #!FileCoercion
		File final_log_output = select_first([alignment.final_log, star_final_log]) #!FileCoercion
		File progress_log_output = select_first([alignment.progress_log, star_progress_log]) #!FileCoercion
		File sj_out_tab_output = select_first([alignment.sj_out_tab, star_sj_out_tab]) #!FileCoercion

		String salmon_quant_file = "~{salmon_raw_data_path}/~{sample.sample_id}.alignment_mode.quant.sf"
		String salmon_command_info_json = "~{salmon_raw_data_path}/~{sample.sample_id}.alignment_mode.cmd_info.json"
		String salmon_aux_info_tar_gz = "~{salmon_raw_data_path}/~{sample.sample_id}.alignment_mode.aux_info.tar.gz"

		if (quantification_complete == "false") {
			call quantification {
				input:
					sample_id = sample.sample_id,
					transcripts_fasta = reference.transcripts_fasta,
					aligned_bam = aligned_bam_output.data,
					aligned_bam_index = aligned_bam_output.data_index,
					raw_data_path = salmon_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File quant_file_output = select_first([quantification.quant_file, salmon_quant_file]) #!FileCoercion
		File command_info_json_output = select_first([quantification.command_info_json, salmon_command_info_json]) #!FileCoercion
		File aux_info_tar_gz_output = select_first([quantification.aux_info_tar_gz, salmon_aux_info_tar_gz]) #!FileCoercion
	}

	output {
		# STAR alignment
		Array[IndexData] aligned_bams = aligned_bam_output #!FileCoercion
		Array[File] unmapped_mate1 = unmapped_mate1_output #!FileCoercion
		Array[File] unmapped_mate2 = unmapped_mate2_output #!FileCoercion
		Array[File] log = log_output #!FileCoercion
		Array[File] final_log = final_log_output #!FileCoercion
		Array[File] progress_log = progress_log_output #!FileCoercion
		Array[File] sj_out_tab = sj_out_tab_output #!FileCoercion

		# Salmon quantification
		Array[File] quant_file = quant_file_output #!FileCoercion
		Array[File] command_info_json = command_info_json_output #!FileCoercion
		Array[File] aux_info_tar_gz = aux_info_tar_gz_output #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] alignment_output_files
		Array[String] quantification_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			alignment_file=$(echo "${output_files}" | cut -f 1)
			quantification_file=$(echo "${output_files}" | cut -f 2)

			if gsutil -u ~{billing_project} ls "${alignment_file}"; then
				if gsutil -u ~{billing_project} ls "${quantification_file}"; then
					# If we find all outputs, don't rerun anything
					echo -e "true\ttrue" >> sample_preprocessing_complete.tsv
				else
					# If we find alignment outputs, but not quantification outputs, then run (or rerun) quantification
					echo -e "true\tfalse" >> sample_preprocessing_complete.tsv
				fi
			else
				# If we don't find alignment outputs, we must also need to run (or rerun) alignment_quantification
				echo -e "false\tfalse" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(alignment_output_files)} ~{write_lines(quantification_output_files)})
	>>>

	output {
		Array[Array[String]] sample_preprocessing_complete = read_tsv("sample_preprocessing_complete.tsv")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
		zones: zones
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
			--output salmon_quant \
			--validateMappings \
			--threads ~{threads}

		cd salmon_quant
		mv quant.sf ~{sample_id}.alignment_mode.quant.sf
		mv cmd_info.json ~{sample_id}.alignment_mode.cmd_info.json
		tar -czvf "~{sample_id}.alignment_mode.aux_info.tar.gz" "aux_info"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.alignment_mode.quant.sf" \
			-o "~{sample_id}.alignment_mode.cmd_info.json" \
			-o "~{sample_id}.alignment_mode.aux_info.tar.gz"
	>>>

	output {
		String quant_file = "~{raw_data_path}/~{sample_id}.alignment_mode.quant.sf"
		String command_info_json = "~{raw_data_path}/~{sample_id}.alignment_mode.cmd_info.json"
		String aux_info_tar_gz = "~{raw_data_path}/~{sample_id}.alignment_mode.aux_info.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/salmon:1.10.1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
