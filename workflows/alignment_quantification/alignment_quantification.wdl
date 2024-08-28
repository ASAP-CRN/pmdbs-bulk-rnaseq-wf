version 1.0

# Align reads

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
		String alignment_output = "~{star_raw_data_path}/~{sample_object.sample_id}_Aligned.sortedByCoord.out.bam"
		String quantification_output = "~{salmon_raw_data_path}/" #TODO
	}

	# For each sample, outputs an array of true/false: [alignment_complete, quantification_complete]
	call check_output_files_exist {
		input:
			alignment_output_files = alignment_output,
			quantification_output_files = quantification_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (sample_index in range(length(samples))) {
		Sample sample = samples[sample_index]

		String alignment_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][0]
		String quantification_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][1]

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

		String star_aligned_bam = "~{star_raw_data_path}/~{sample.sample_id}_Aligned.sortedByCoord.out.bam"
		String star_aligned_bam_index = "~{star_raw_data_path}/~{sample.sample_id}_Aligned.sortedByCoord.out.bam.bai"
		String star_unmapped_mate1 = "~{star_raw_data_path}/~{sample.sample_id}_Unmapped.out.mate1"
		String star_unmapped_mate2 = "~{star_raw_data_path}/~{sample.sample_id}_Unmapped.out.mate2"
		String star_log = "~{star_raw_data_path}/~{sample.sample_id}_Log.out"
		String star_final_log = "~{star_raw_data_path}/~{sample.sample_id}_Log.final.out"
		String star_progress_log = "~{star_raw_data_path}/~{sample.sample_id}_Log.progress.out"
		String star_sj_out_tab = "~{star_raw_data_path}/~{sample.sample_id}_SJ.out.tab"

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
