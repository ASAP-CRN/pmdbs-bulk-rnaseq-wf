version 1.0

struct Sample {
	String sample_id
	String? batch

	Array[File]+ fastq_R1s
	Array[File]+ fastq_R2s
	Array[File] fastq_I1s
	Array[File] fastq_I2s
}

struct Project {
	String team_id
	String dataset_id
	String dataset_doi_url
	Array[Sample] samples

	File project_sample_metadata_csv

	Boolean run_project_cohort_analysis

	String raw_data_bucket
	Array[String] staging_data_buckets
}

struct ReferenceData {
	File primary_assembly_fasta
	File gene_annotation_gtf
	File transcripts_fasta
	File all_transcripts_fasta
}
