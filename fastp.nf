#!/usr/bin/env nextflow

params.outdir = "/proj/breedmap/trimmed_reads"
reads_ch = channel.fromFilePairs('/proj/breedmap/hgen_bull_fastq/BTA126_S2_L001_R{1,2}_001.fastq.gz')

process fastp {
	publishDir "${params.outdir}/trimmed_BTA126", mode: 'copy'
	
	input:
	tuple sample_id, file(x) from reads_ch

	output:
	file("*") 

	script:
	"""
	module load fastp
	fastp -i ${x[0]} -I ${x[1]} -o ${sample_id}1.fq.gz -O ${sample_id}2.fq.gz --detect_adapter_for_pe -h ${sample_id}.html -j ${sample_id}.json
	
	"""


}
