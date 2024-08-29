version 1.0

# Perform QC and trim reads

import "../pmdbs_bulk_rnaseq_analysis_structs.wdl"
import "../wf-common/wdl/tasks/fastqc.wdl" as Fastqc

workflow preprocess {
	input {
		String project_id
		Array[Sample] samples

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
	String sub_workflow_name = "preprocess"
	String sub_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String fastqc_raw_reads_raw_data_path = "~{workflow_raw_data_path_prefix}/fastqc_raw_reads/~{sub_workflow_version}"
	String fastp_raw_data_path = "~{workflow_raw_data_path_prefix}/trim_and_qc/~{sub_workflow_version}"
	String fastqc_trimmed_reads_raw_data_path = "~{workflow_raw_data_path_prefix}/fastqc_trimmed_reads/~{sub_workflow_version}"

	scatter (sample_object in samples) {
		String fastqc_raw_reads_output = "~{fastqc_raw_reads_raw_data_path}/~{sample_object.sample_id}_fastqc_reports.tar.gz"
		String fastqc_trimmed_reads_output = "~{fastqc_trimmed_reads_raw_data_path}/~{sample_object.sample_id}_trimmed_fastqc_reports.tar.gz"
	}

	# For each sample, outputs an array of true/false: [fastqc_raw_reads_complete, fastqc_trimmed_reads_complete]
	call check_output_files_exist {
		input:
			fastqc_raw_reads_output_files = fastqc_raw_reads_output,
			fastqc_trimmed_reads_output_files = fastqc_trimmed_reads_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (sample_index in range(length(samples))) {
		Sample sample = samples[sample_index]

		Array[String] project_sample_id = [project_id, sample.sample_id]

		String fastqc_raw_reads_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][0]
		String fastqc_trimmed_reads_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][1]

		String fastqc_raw_reads_reports_tar_gz = "~{fastqc_raw_reads_raw_data_path}/~{sample.sample_id}_fastqc_reports.tar.gz"

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
			String fastq_basename = sub(basename(fastq_R1_basename), "(R1|r1)[_-]*$", "")

			String fastp_trimmed_fastq_R1 = "~{fastp_raw_data_path}/~{fastq_R1_basename}_trimmed.fastq.gz"
			String fastp_trimmed_fastq_R2 = "~{fastp_raw_data_path}/~{fastq_R2_basename}_trimmed.fastq.gz"
			String fastp_failed_paired_fastq = "~{fastp_raw_data_path}/~{fastq_basename}_failed.fastq.gz"
			String fastp_report_html = "~{fastp_raw_data_path}/~{fastq_basename}_fastp.html"

			if (fastqc_trimmed_reads_complete == "false") {
				call trim_and_qc {
					input:
						fastq_R1 = sample.fastq_R1s[fastq_index],
						fastq_R2 = sample.fastq_R2s[fastq_index],
						fastq_R1_basename = fastq_R1_basename,
						fastq_R2_basename = fastq_R2_basename,
						fastq_basename = fastq_basename,
						raw_data_path = fastp_raw_data_path,
						workflow_info = workflow_info,
						billing_project = billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}

			File selected_trimmed_fastq_R1 = select_first([trim_and_qc.trimmed_fastq_R1, fastp_trimmed_fastq_R1]) #!FileCoercion
		    File selected_trimmed_fastq_R2 = select_first([trim_and_qc.trimmed_fastq_R2, fastp_trimmed_fastq_R2]) #!FileCoercion
		    File selected_failed_paired_fastq = select_first([trim_and_qc.failed_paired_fastq, fastp_failed_paired_fastq]) #!FileCoercion
		    File selected_report_html = select_first([trim_and_qc.report_html, fastp_report_html]) #!FileCoercion
		}

		Array[File] trimmed_fastq_R1s_output = selected_trimmed_fastq_R1
		Array[File] trimmed_fastq_R2s_output = selected_trimmed_fastq_R2
		Array[File] failed_paired_fastq_output = selected_failed_paired_fastq
		Array[File] report_html_output = selected_report_html

		String fastqc_trimmed_reads_reports_tar_gz = "~{fastqc_raw_reads_raw_data_path}/~{sample.sample_id}_trimmed_fastqc_reports.tar.gz"

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
		
		TrimmedSample trimmed_fastq_samples_output = {
			"sample_id": sample.sample_id,
			"batch": sample.batch,
			"fastq_R1s": trimmed_fastq_R1s_output,
			"fastq_R2s": trimmed_fastq_R2s_output,
			"fastq_I1s": sample.fastq_I1s,
			"fastq_I2s": sample.fastq_I2s
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
		Array[Array[File]] failed_paired_fastq = failed_paired_fastq_output #!FileCoercion
		Array[Array[File]] report_html = report_html_output #!FileCoercion

		# FastQC reports on trimmed reads
		Array[File] trimmed_fastqc_reports_tar_gz = trimmed_fastqc_reports_tar_gz_output #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] fastqc_raw_reads_output_files
		Array[String] fastqc_trimmed_reads_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			fastqc_raw_reads_reports_file=$(echo "${output_files}" | cut -f 1)
			fastqc_trimmed_reads_reports_file=$(echo "${output_files}" | cut -f 2)

			if gsutil -u ~{billing_project} ls "${fastqc_raw_reads_reports_file}"; then
				if gsutil -u ~{billing_project} ls "${fastqc_trimmed_reads_reports_file}"; then
					# If we find all outputs, don't rerun anything
					echo -e "true\ttrue" >> sample_preprocessing_complete.tsv
				else
					# If we find fastqc outputs for raw reads, but not trimmed reads, then run (or rerun) trim_and_qc and fastqc_trimmed_reads
					echo -e "true\tfalse" >> sample_preprocessing_complete.tsv
				fi
			else
				# If we don't find fastqc output for raw reads, we must also need to run (or rerun) preprocessing
				echo -e "false\tfalse" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(fastqc_raw_reads_output_files)} ~{write_lines(fastqc_trimmed_reads_output_files)})
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
		File fastq_R1
		File fastq_R2

		String fastq_R1_basename
		String fastq_R2_basename
		String fastq_basename

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 8
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([fastq_R1, fastq_R2], "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		fastp \
			--in1 ~{fastq_R1} \
			--out1 ~{fastq_R1_basename}_trimmed.fastq.gz \
			--in2 ~{fastq_R2} \
			--out2 ~{fastq_R2_basename}_trimmed.fastq.gz \
			--failed_out ~{fastq_basename}_failed.fastq.gz \
			--detect_adapter_for_pe \
			--correction \
			--overrepresentation_analysis \
			--html ~{fastq_basename}_fastp.html \
			--report_title ~{fastq_basename} \
			--thread ~{threads - 1}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{fastq_R1_basename}_trimmed.fastq.gz" \
			-o "~{fastq_R2_basename}_trimmed.fastq.gz" \
			-o "~{fastq_basename}_failed.fastq.gz" \
			-o "~{fastq_basename}_fastp.html"
	>>>

	output {
		String trimmed_fastq_R1 = "~{raw_data_path}/~{fastq_R1_basename}_trimmed.fastq.gz"
		String trimmed_fastq_R2 = "~{raw_data_path}/~{fastq_R2_basename}_trimmed.fastq.gz"
		String failed_paired_fastq = "~{raw_data_path}/~{fastq_basename}_failed.fastq.gz"
		String report_html = "~{raw_data_path}/~{fastq_basename}_fastp.html"
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
