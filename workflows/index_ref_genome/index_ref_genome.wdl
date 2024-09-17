version 1.0

# Index reference genome for STAR and/or Salmon (pseudo)mapping-mode

import "../pmdbs_bulk_rnaseq_analysis_structs.wdl"

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
			--sjdbGTFfile ~{gene_annotation_gtf} \
			--readFilesCommand zcat

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

task generate_decoy {
	input {
		File primary_assembly_fasta
		File transcripts_fasta
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([primary_assembly_fasta, transcripts_fasta], "GB") *2 + 100)

	command <<<
		set -euo pipefail

		grep "^>" <(gunzip -c ~{primary_assembly_fasta}) | cut -d " " -f 1 > decoys.txt
		sed -i.bak -e 's/>//g' decoys.txt
		zcat ~{transcripts_fasta} ~{primary_assembly_fasta} | gzip > gentrome.fa.gz
	>>>

	output {
		File decoys_txt = "decoys.txt"
		File gentrome_fasta = "gentrome.fa.gz"
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

task salmon_index_ref_genome {
	input {
		File gentrome_fasta
		File decoys_txt
		String container_registry
		String zones
	}

	Int threads = 24
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([gentrome_fasta, decoys_txt], "GB") *2 + 100)

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
		docker: "~{container_registry}/salmon:1.10.1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
