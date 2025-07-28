version 1.0

# Index reference genome for STAR and/or Salmon (pseudo)mapping-mode

import "../structs.wdl"

workflow index_ref_genome {
	input {
		ReferenceData reference

		Boolean run_star_index_ref_genome
		Boolean run_salmon_index_ref_genome

		String container_registry
		String zones
	}

	# Prepare reference files for alignment/mapping and quantification
	if (run_star_index_ref_genome) {
		call star_index_ref_genome {
		input:
			primary_assembly_fasta = reference.primary_assembly_fasta,
			gene_annotation_gtf = reference.gene_annotation_gtf,
			container_registry = container_registry,
			zones = zones
		}
	}

	if (run_salmon_index_ref_genome) {
		call generate_decoy {
		input:
			primary_assembly_fasta = reference.primary_assembly_fasta,
			transcripts_fasta = reference.transcripts_fasta,
			container_registry = container_registry,
			zones = zones
		}

		call salmon_index_ref_genome {
		input:
			gentrome_fasta = generate_decoy.gentrome_fasta,
			decoys_txt = generate_decoy.decoys_txt,
			container_registry = container_registry,
			zones = zones
		}
	}

	output {
		File? star_genome_dir_tar_gz = star_index_ref_genome.star_genome_dir_tar_gz
		File? salmon_genome_dir_tar_gz = salmon_index_ref_genome.salmon_genome_dir_tar_gz
	}
}

task star_index_ref_genome {
	input {
		File primary_assembly_fasta
		File gene_annotation_gtf
		String container_registry
		String zones
	}

	String primary_assembly_fasta_basename = basename(primary_assembly_fasta, ".gz")
	String gene_annotation_gtf_basename = basename(gene_annotation_gtf, ".gz")

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([primary_assembly_fasta, gene_annotation_gtf], "GB") * 2 + 100)

	command <<<
		set -euo pipefail

		gunzip ~{primary_assembly_fasta} ~{gene_annotation_gtf}

		ref_path=$(dirname ~{primary_assembly_fasta})

		STAR \
			--runThreadN ~{threads - 1} \
			--runMode genomeGenerate \
			--genomeDir star_genome_dir \
			--genomeFastaFiles "$ref_path"/~{primary_assembly_fasta_basename} \
			--sjdbGTFfile "$ref_path"/~{gene_annotation_gtf_basename}

		tar -czvf star_genome_dir.tar.gz star_genome_dir
	>>>

	output {
		File star_genome_dir_tar_gz = "star_genome_dir.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/star_samtools:2.7.11b_1.22.1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		maxRetries: 2
		zones: zones
	}
}

task generate_decoy {
	input {
		File primary_assembly_fasta
		File transcripts_fasta
		String container_registry
		String zones
	}

	String ref_basename = basename(transcripts_fasta, ".transcripts.fa.gz")

	Int threads = 2
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([primary_assembly_fasta, transcripts_fasta], "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		grep "^>" <(gunzip -c ~{primary_assembly_fasta}) | cut -d " " -f 1 > ~{ref_basename}.decoys.txt
		sed -i -e 's/>//g' ~{ref_basename}.decoys.txt
		cat ~{transcripts_fasta} ~{primary_assembly_fasta} > ~{ref_basename}.gentrome.fa.gz
	>>>

	output {
		File decoys_txt = "~{ref_basename}.decoys.txt"
		File gentrome_fasta = "~{ref_basename}.gentrome.fa.gz"
	}

	runtime {
		docker: "~{container_registry}/salmon:1.10.3"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		maxRetries: 2
		zones: zones
	}
}

task salmon_index_ref_genome {
	input {
		File gentrome_fasta
		File decoys_txt
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([gentrome_fasta, decoys_txt], "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		salmon index \
			--transcripts ~{gentrome_fasta} \
			--index salmon_index \
			--kmerLen 31 \
			--gencode \
			--decoys ~{decoys_txt} \
			--threads ~{threads}

		tar -czvf salmon_genome_dir.tar.gz salmon_index
	>>>

	output {
		File salmon_genome_dir_tar_gz = "salmon_genome_dir.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/salmon:1.10.3"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		maxRetries: 2
		zones: zones
	}
}
