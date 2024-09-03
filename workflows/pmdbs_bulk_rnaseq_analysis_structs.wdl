version 1.0

import "wf-common/wdl/structs.wdl"

struct ReferenceData {
	File primary_assembly_fasta
	File gene_annotation_gtf
	File transcripts_fasta
}

struct OutputSample {
	Array[String] sample_id
	Array[String]? batch

	Array[File]+ fastq_R1s
	Array[File]+ fastq_R2s
	Array[File] fastq_I1s
	Array[File] fastq_I2s
}
