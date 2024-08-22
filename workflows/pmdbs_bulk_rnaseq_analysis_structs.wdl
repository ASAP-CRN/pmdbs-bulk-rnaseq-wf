version 1.0

import "wf-common/wdl/structs.wdl"

struct ReferenceData {
	File primary_assembly_fasta
	File gene_annotation_gtf
	File transcripts_fasta
}
