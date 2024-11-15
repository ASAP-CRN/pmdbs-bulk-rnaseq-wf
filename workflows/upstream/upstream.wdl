version 1.0

# Perform QC, trim, align, and quantify reads

import "../pmdbs_bulk_rnaseq_analysis_structs.wdl"
import "../../wf-common/wdl/tasks/fastqc.wdl" as Fastqc
import "alignment_quantification/alignment_quantification.wdl" as AlignmentQuantification
import "pseudo_mapping_quantification/pseudo_mapping_quantification.wdl" as PseudoMappingQuantification

workflow upstream {
	input {
		String team_id
		Array[Sample] samples

		File all_transcripts_fasta

		# STAR and salmon quantification option
		Boolean run_alignment_quantification
		File? star_genome_dir_tar_gz

		# Salmon mapping and quantification option
		Boolean run_pseudo_mapping_quantification
		File? salmon_genome_dir_tar_gz

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
	String sub_workflow_name = "upstream"
	String fastqc_task_version = "1.0.0"
	String trim_and_qc_task_version = "1.0.1"
	String alignment_quantification_workflow_version = "1.0.0"
	String pseudo_mapping_quantification_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String fastqc_raw_reads_raw_data_path = "~{workflow_raw_data_path_prefix}/fastqc_raw_reads/~{fastqc_task_version}"
	String fastp_raw_data_path = "~{workflow_raw_data_path_prefix}/trim_and_qc/~{trim_and_qc_task_version}"
	String fastqc_trimmed_reads_raw_data_path = "~{workflow_raw_data_path_prefix}/fastqc_trimmed_reads/~{fastqc_task_version}"
	String star_and_salmon_alignment_mode_raw_data_path = "~{workflow_raw_data_path_prefix}/alignment_quantification/~{alignment_quantification_workflow_version}"
	String salmon_mapping_mode_raw_data_path = "~{workflow_raw_data_path_prefix}/mapping_quantification/~{pseudo_mapping_quantification_workflow_version}"

	scatter (sample_object in samples) {
		String fastqc_raw_reads_output = "~{fastqc_raw_reads_raw_data_path}/~{sample_object.sample_id}.fastqc_reports.tar.gz"
		String fastqc_trimmed_reads_output = "~{fastqc_trimmed_reads_raw_data_path}/~{sample_object.sample_id}.trimmed_fastqc_reports.tar.gz"
		String alignment_quantification_output = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample_object.sample_id}.alignment_mode.salmon_quant.tar.gz"
		String pseudo_mapping_quantification_output = "~{salmon_mapping_mode_raw_data_path}/~{sample_object.sample_id}.mapping_mode.salmon_quant.tar.gz"
	}

	# For each sample, outputs an array of true/false: [fastqc_raw_reads_complete, fastqc_trimmed_reads_complete, alignment_quantification_complete, pseudo_mapping_quantification_complete]
	call check_output_files_exist {
		input:
			fastqc_raw_reads_output_files = fastqc_raw_reads_output,
			fastqc_trimmed_reads_output_files = fastqc_trimmed_reads_output,
			alignment_quantification_output_files = alignment_quantification_output,
			pseudo_mapping_quantification_output_files = pseudo_mapping_quantification_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (sample_index in range(length(samples))) {
		Sample sample = samples[sample_index]

		Array[String] project_sample_id = [team_id, sample.sample_id]

		String fastqc_raw_reads_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][0]
		String fastqc_trimmed_reads_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][1]
		String alignment_quantification_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][2]
		String pseudo_mapping_quantification_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][3]

		String fastqc_raw_reads_reports_tar_gz = "~{fastqc_raw_reads_raw_data_path}/~{sample.sample_id}.fastqc_reports.tar.gz"

		if (fastqc_raw_reads_complete == "false") {
			call Fastqc.fastqc as fastqc_raw_reads {
				input:
					sample_id = sample.sample_id,
					fastq_R1s = sample.fastq_R1s,
					fastq_R2s = sample.fastq_R2s,
					raw_data_path = fastqc_raw_reads_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File fastqc_reports_tar_gz_output = select_first([fastqc_raw_reads.fastqc_reports_tar_gz, fastqc_raw_reads_reports_tar_gz]) #!FileCoercion

		scatter (fastq_index in range(length(sample.fastq_R1s))) {
			String fastq_R1_basename = sub(basename(sample.fastq_R1s[fastq_index]), "(\\.fastq|\\.fq)(\\.gz|)", "")
			String fastq_R2_basename = sub(basename(sample.fastq_R2s[fastq_index]), "(\\.fastq|\\.fq)(\\.gz|)", "")

			String fastp_trimmed_fastq_R1 = "~{fastp_raw_data_path}/~{fastq_R1_basename}.trimmed.fastq.gz"
			String fastp_trimmed_fastq_R2 = "~{fastp_raw_data_path}/~{fastq_R2_basename}.trimmed.fastq.gz"
		}

		String fastp_failed_paired_fastqs = "~{fastp_raw_data_path}/~{sample.sample_id}.fastp_failed_paired_fastqs.tar.gz"
		String fastp_reports_html = "~{fastp_raw_data_path}/~{sample.sample_id}.fastp_reports.tar.gz"
		String fastp_jsons = "~{fastp_raw_data_path}/~{sample.sample_id}.fastp_json.tar.gz"

		if (fastqc_trimmed_reads_complete == "false") {
			call trim_and_qc {
				input:
					sample_id = sample.sample_id,
					fastq_R1s = sample.fastq_R1s,
					fastq_R2s = sample.fastq_R2s,
					raw_data_path = fastp_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		Array[File] trimmed_fastq_R1s_output = if (fastqc_trimmed_reads_complete == "false") then select_first([trim_and_qc.trimmed_fastq_R1s]) else fastp_trimmed_fastq_R1 #!FileCoercion
		Array[File] trimmed_fastq_R2s_output = if (fastqc_trimmed_reads_complete == "false") then select_first([trim_and_qc.trimmed_fastq_R2s]) else fastp_trimmed_fastq_R2 #!FileCoercion
		File qc_failed_paired_fastqs_tar_gz_output = select_first([trim_and_qc.qc_failed_paired_fastqs_tar_gz, fastp_failed_paired_fastqs]) #!FileCoercion
		File qc_reports_html_tar_gz_output = select_first([trim_and_qc.qc_reports_html_tar_gz, fastp_reports_html]) #!FileCoercion
		File qc_json_tar_gz_output = select_first([trim_and_qc.qc_json_tar_gz, fastp_jsons]) #!FileCoercion

		String fastqc_trimmed_reads_reports_tar_gz = "~{fastqc_trimmed_reads_raw_data_path}/~{sample.sample_id}.trimmed_fastqc_reports.tar.gz"

		if (fastqc_trimmed_reads_complete == "false") {
			call Fastqc.fastqc as fastqc_trimmed_reads {
				input:
					sample_id = sample.sample_id,
					fastq_R1s = trimmed_fastq_R1s_output,
					fastq_R2s = trimmed_fastq_R2s_output,
					raw_data_path = fastqc_trimmed_reads_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File trimmed_fastqc_reports_tar_gz_output = select_first([fastqc_trimmed_reads.trimmed_fastqc_reports_tar_gz, fastqc_trimmed_reads_reports_tar_gz]) #!FileCoercion

		if (run_alignment_quantification) {
			String star_aligned_bam = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.Aligned.sortedByCoord.out.bam"
			String star_aligned_bam_index = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.Aligned.sortedByCoord.out.bam.bai"
			String star_aligned_to_transcriptome_bam = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.Aligned.toTranscriptome.out.bam"
			String star_unmapped_mate1 = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.Unmapped.out.mate1"
			String star_unmapped_mate2 = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.Unmapped.out.mate2"
			String star_log = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.Log.out"
			String star_final_log = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.Log.final.out"
			String star_progress_log = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.Log.progress.out"
			String star_sj_out_tab = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.SJ.out.tab"
			String salmon_alignment_mode_quant_tar_gz = "~{star_and_salmon_alignment_mode_raw_data_path}/~{sample.sample_id}.alignment_mode.salmon_quant.tar.gz"

			if (alignment_quantification_complete == "false") {
				call AlignmentQuantification.alignment_quantification {
					input:
						sample_id = sample.sample_id,
						all_transcripts_fasta = all_transcripts_fasta,
						star_genome_dir_tar_gz = select_first([star_genome_dir_tar_gz]),
						trimmed_fastq_R1s = trimmed_fastq_R1s_output,
						trimmed_fastq_R2s = trimmed_fastq_R2s_output,
						raw_data_path = star_and_salmon_alignment_mode_raw_data_path,
						workflow_info = workflow_info,
						billing_project = billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}

			File aligned_bam_output = select_first([alignment_quantification.aligned_bam, star_aligned_bam]) #!FileCoercion
			File aligned_bam_index_output = select_first([alignment_quantification.aligned_bam_index, star_aligned_bam_index]) #!FileCoercion
			File aligned_to_transcriptome_bam_output = select_first([alignment_quantification.aligned_to_transcriptome_bam, star_aligned_to_transcriptome_bam]) #!FileCoercion
			File unmapped_mate1_output = select_first([alignment_quantification.unmapped_mate1, star_unmapped_mate1]) #!FileCoercion
			File unmapped_mate2_output = select_first([alignment_quantification.unmapped_mate2, star_unmapped_mate2]) #!FileCoercion
			File log_output = select_first([alignment_quantification.log, star_log]) #!FileCoercion
			File final_log_output = select_first([alignment_quantification.final_log, star_final_log]) #!FileCoercion
			File progress_log_output = select_first([alignment_quantification.progress_log, star_progress_log]) #!FileCoercion
			File sj_out_tab_output = select_first([alignment_quantification.sj_out_tab, star_sj_out_tab]) #!FileCoercion
			File salmon_alignment_mode_quant_tar_gz_output = select_first([alignment_quantification.quant_tar_gz, salmon_alignment_mode_quant_tar_gz]) #!FileCoercion
		}

		if (run_pseudo_mapping_quantification) {
			String salmon_mapping_mode_quant_tar_gz = "~{salmon_mapping_mode_raw_data_path}/~{sample.sample_id}.mapping_mode.salmon_quant.tar.gz"

			if (pseudo_mapping_quantification_complete == "false") {
				call PseudoMappingQuantification.pseudo_mapping_quantification {
					input:
						sample_id = sample.sample_id,
						salmon_genome_dir_tar_gz = select_first([salmon_genome_dir_tar_gz]),
						trimmed_fastq_R1s = trimmed_fastq_R1s_output,
						trimmed_fastq_R2s = trimmed_fastq_R2s_output,
						raw_data_path = salmon_mapping_mode_raw_data_path,
						workflow_info = workflow_info,
						billing_project = billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}

			File salmon_mapping_mode_quant_tar_gz_output = select_first([pseudo_mapping_quantification.quant_tar_gz, salmon_mapping_mode_quant_tar_gz]) #!FileCoercion
		}
	}

	output {
		# Sample list
		Array[Array[String]] project_sample_ids = project_sample_id

		# FastQC reports on raw reads
		Array[File] fastqc_reports_tar_gz = fastqc_reports_tar_gz_output #!FileCoercion

		# Trimmed reads
		Array[Array[File]] trimmed_fastq_R1s = trimmed_fastq_R1s_output #!FileCoercion
		Array[Array[File]] trimmed_fastq_R2s = trimmed_fastq_R2s_output #!FileCoercion
		Array[File] qc_failed_paired_fastqs_tar_gz = qc_failed_paired_fastqs_tar_gz_output #!FileCoercion
		Array[File] qc_reports_html_tar_gz = qc_reports_html_tar_gz_output #!FileCoercion
		Array[File] qc_json_tar_gz = qc_json_tar_gz_output #!FileCoercion

		# FastQC reports on trimmed reads
		Array[File] trimmed_fastqc_reports_tar_gz = trimmed_fastqc_reports_tar_gz_output #!FileCoercion

		# STAR alignment
		Array[File?] aligned_bam = aligned_bam_output #!FileCoercion
		Array[File?] aligned_bam_index = aligned_bam_index_output #!FileCoercion
		Array[File?] aligned_to_transcriptome_bam = aligned_to_transcriptome_bam_output #!FileCoercion
		Array[File?] unmapped_mate1 = unmapped_mate1_output #!FileCoercion
		Array[File?] unmapped_mate2 = unmapped_mate2_output #!FileCoercion
		Array[File?] log = log_output #!FileCoercion
		Array[File?] final_log = final_log_output #!FileCoercion
		Array[File?] progress_log = progress_log_output #!FileCoercion
		Array[File?] sj_out_tab = sj_out_tab_output #!FileCoercion
		# Salmon alignment-mode quantification
		Array[File?] alignment_mode_quant_tar_gz = salmon_alignment_mode_quant_tar_gz_output #!FileCoercion

		# Salmon mapping and quantification
		Array[File?] mapping_mode_quant_tar_gz = salmon_mapping_mode_quant_tar_gz_output #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] fastqc_raw_reads_output_files
		Array[String] fastqc_trimmed_reads_output_files
		Array[String] alignment_quantification_output_files
		Array[String] pseudo_mapping_quantification_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			fastqc_raw_reads_reports_file=$(echo "${output_files}" | cut -f 1)
			fastqc_trimmed_reads_reports_file=$(echo "${output_files}" | cut -f 2)
			alignment_quantification_file=$(echo "${output_files}" | cut -f 3)
			pseudo_mapping_quantification_file=$(echo "${output_files}" | cut -f 4)

			if gsutil -u ~{billing_project} ls "${fastqc_raw_reads_reports_file}"; then
				if gsutil -u ~{billing_project} ls "${fastqc_trimmed_reads_reports_file}"; then
					if gsutil -u ~{billing_project} ls "${alignment_quantification_file}"; then
						if gsutil -u ~{billing_project} ls "${pseudo_mapping_quantification_file}"; then
							# If we find all outputs, don't rerun anything
							echo -e "true\ttrue\ttrue\ttrue" >> sample_preprocessing_complete.tsv
						else
							# If we find all outputs except pseudo_mapping_quantification outputs, then run (or rerun) pseudo_mapping_quantification
							echo -e "true\ttrue\ttrue\tfalse" >> sample_preprocessing_complete.tsv
						fi
					else
						if gsutil -u ~{billing_project} ls "${pseudo_mapping_quantification_file}"; then
							# If we find all outputs except alignment_quantification outputs, then run (or rerun) alignment_quantification
							echo -e "true\ttrue\tfalse\ttrue" >> sample_preprocessing_complete.tsv
						else
							# If we can't find any alignment/mapping and quantification files, then run (or rerun) alignment_quantification and/or pseudo_mapping_quantification
							echo -e "true\ttrue\tfalse\tfalse" >> sample_preprocessing_complete.tsv
						fi
					fi
				else
					# If we can't find trimmed reads, then run (or rerun) trimming, aligning/mapping, and quantification
					echo -e "true\tfalse\tfalse\tfalse" >> sample_preprocessing_complete.tsv
				fi
			else
				# If we can't find fastqc reports on raw reads, then run (or rerun) everything
				echo -e "false\tfalse\tfalse\tfalse" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(fastqc_raw_reads_output_files)} ~{write_lines(fastqc_trimmed_reads_output_files)} ~{write_lines(alignment_quantification_output_files)} ~{write_lines(pseudo_mapping_quantification_output_files)})
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

task trim_and_qc {
	input {
		String sample_id

		Array[File] fastq_R1s
		Array[File] fastq_R2s

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 16
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(flatten([fastq_R1s, fastq_R2s]), "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		while read -r fastq_R1 fastq_R2 || [[ -n "${fastq_R1}" ]] || [[ -n "${fastq_R2}" ]]; do
			fastq_R1_basename=$(basename "${fastq_R1}" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
			fastq_R2_basename=$(basename "${fastq_R2}" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
			fastq_basename=$(echo "${fastq_R1_basename}" | sed -E 's/(R1|r1)[_-]*$//; s/[_-]+$//')

			fastp \
				--in1 "${fastq_R1}" \
				--out1 "${fastq_R1_basename}.trimmed.fastq.gz" \
				--in2 "${fastq_R2}" \
				--out2 "${fastq_R2_basename}.trimmed.fastq.gz" \
				--failed_out "${fastq_basename}.failed.fastq.gz" \
				--detect_adapter_for_pe \
				--correction \
				--overrepresentation_analysis \
				--html "${fastq_basename}.fastp.html" \
				--report_title "${fastq_basename}" \
				--thread ~{threads - 1}

			mv fastp.json "${fastq_basename}.fastp.json"

			# Trimmed fastqs will live alongside sample-level upstream outputs in workflow execution raw bucket
			upload_outputs \
				-b ~{billing_project} \
				-d ~{raw_data_path} \
				-i ~{write_tsv(workflow_info)} \
				-o "${fastq_R1_basename}.trimmed.fastq.gz" \
				-o "${fastq_R2_basename}.trimmed.fastq.gz"

			# Get list of trimmed fastqs in order to save into Array[String]
			echo "~{raw_data_path}/${fastq_R1_basename}.trimmed.fastq.gz" > trimmed_fastq_R1s_path.txt
			echo "~{raw_data_path}/${fastq_R2_basename}.trimmed.fastq.gz" > trimmed_fastq_R2s_path.txt
		done < <(paste ~{write_lines(fastq_R1s)} ~{write_lines(fastq_R2s)})

		mkdir -p ~{sample_id}_fastp_failed_paired_fastqs
		find . -maxdepth 1 -type f -name "*failed.fastq.gz" -exec mv {} ~{sample_id}_fastp_failed_paired_fastqs/ \;
		tar -czvf "~{sample_id}.fastp_failed_paired_fastqs.tar.gz" "~{sample_id}_fastp_failed_paired_fastqs"

		mkdir -p ~{sample_id}_fastp_reports
		find . -maxdepth 1 -type f -name "*fastp.html" -exec mv {} ~{sample_id}_fastp_reports/ \;
		tar -czvf "~{sample_id}.fastp_reports.tar.gz" "~{sample_id}_fastp_reports"

		mkdir -p ~{sample_id}_fastp_json
		find . -maxdepth 1 -type f -name "*fastp.json" -exec mv {} ~{sample_id}_fastp_json/ \;
		tar -czvf "~{sample_id}.fastp_json.tar.gz" "~{sample_id}_fastp_json"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.fastp_failed_paired_fastqs.tar.gz" \
			-o "~{sample_id}.fastp_reports.tar.gz" \
			-o "~{sample_id}.fastp_json.tar.gz"
	>>>

	output {
		Array[String] trimmed_fastq_R1s = read_lines("trimmed_fastq_R1s_path.txt")
		Array[String] trimmed_fastq_R2s = read_lines("trimmed_fastq_R2s_path.txt")
		String qc_failed_paired_fastqs_tar_gz = "~{raw_data_path}/~{sample_id}.fastp_failed_paired_fastqs.tar.gz"
		String qc_reports_html_tar_gz = "~{raw_data_path}/~{sample_id}.fastp_reports.tar.gz"
		String qc_json_tar_gz = "~{raw_data_path}/~{sample_id}.fastp_json.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/fastp:0.23.4"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
