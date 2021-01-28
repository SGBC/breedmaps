#!/usr/bin/env/ nextflow
params.ref = "/proj/breedmap/scripts/indexes/ARS-UCD1.2_Btau5.0.1Y.fa.gz"
reads_ch = channel.fromFilePairs('/proj/breedmap/trimmed_reads/trimmed_BTA125/BTA125_S1_L003_R{1,2}.fq.gz')
params.outdir = "/proj/breedmap/results"
params.REF = "/proj/breedmap/REF/ARS-UCD1.2_Btau5.0.1Y.fa.gz"
/*
process sorting_fastq {
	publishDir "${outdir}/sorted_fastq", mode: 'copy'
	input:
	file(reads) from sort_reads_ch
	output:
	file("*")
	script:
	"""
	zcat ${reads}|paste----|sort -k1,1 -S 3G|tr '\t' '\n'|gzip > sort${reads}

	"""
}

*/


process indexing_aligning {
	publishDir "${params.outdir}/indexed_aligned", mode: 'copy'
	input :
	
	path ref from params.ref
	tuple pair_id, file(reads) from reads_ch
	output: 
	file("*")
	set val(pair_id),file("${pair_id}_aligned_reads.sam") into aligned_ch
	script:
	"""
	module load bwa
	bwa mem /proj/breedmap/scripts/indexes/ARS-UCD1.2_Btau5.0.1Y.fa.gz ${reads} > ${pair_id}_aligned_reads.sam
	
	"""
} 

process samtobam_sort {
	publishDir "${params.outdir}/sam_sorted_bam", mode: 'copy'
	input:
	set pair_id, file(reads) from aligned_ch
	output:
	set pair_id, file("${pair_id}_aligned_reads_sorted.bam"), \
	file("${pair_id}_aligned_reads_sorted.bam.bai") into bam_ch	
	script:
	"""
	module load samtools
	samtools view -h -b -S ${reads} > ${pair_id}_aligned_reads.bam
	samtools sort ${pair_id}_aligned_reads.bam -o ${pair_id}_aligned_reads_sorted.bam
	samtools index ${pair_id}_aligned_reads_sorted.bam
	"""
}




process markedups {
	publishDir "${params.outdir}/markedup", mode:'copy' 
	input:
	set val(pair_id), file(bamfile), file(baifile) from bam_ch
	output:
	set pair_id, file("${pair_id}_sorted_duplicates_rm.bam") into markedup_bam_ch
	set pair_id, file("${pair_id}_sorted_duplicates_rm.bam") into markedup_bam_ch2
	set pair_id, file("${pair_id}_markedup_metrics.txt") into marked_metrics_ch
	file("${pair_id}_sorted_duplicates_rm.bam.bai") into marked_index_ch
	script:
	"""
	java -jar /opt/sw/picard/1.137/picard.jar MarkDuplicates I=${bamfile} O=${pair_id}_sorted_duplicates_rm.bam M=${pair_id}_markedup_metrics.txt 
	module load samtools
	samtools index ${pair_id}_sorted_duplicates_rm.bam
	
	"""
}

/*
process BQSR_table {
	publishDir "${params.outdir}/BQSR", mode: 'copy'
	input:
	set pair_id, path(bqsr) from markedup_bam_ch
	path ref from params.REF
	path known_vcf from params.known_sites
	output:
	file("${pair_id}.recal.table") into recal_table_ch
	file("${pair_id}.recal.table") into recal_table_ch2
	script:
	"""
	module load gatk
	java -Xmx80G jar /opt/sw/gatk/3.4.42/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 \
	-R ${ref} -I ${bqsr} -knownSites:vcf ${known_vcf} --bqsrBAQGapOpenPenalty 45 \
	-o ${pair_id}.recal.table
	"""

}		


process PrintReads {
	publishDir "${params.outdir}/BQSR", mode:'copy'
	input:
	set pair_id, path(reads) from marked_bam_ch2
	path ref from params.REF
	file(recal_table) from recal_table_ch
	output:
	set pair_id, file("${pair_id}_dedup_recal.bam") into printreads_ch
	script:
	"""
	module load gatk 
	java -Xmx80G -jar /opt/sw/gatk/3.4.42/GenomeAnalysisTK.jar -T PrintReads -nct 8 \
	-R ${ref} -I ${reads} -BQSR ${recal_table} -o ${pair_id}_dedup_recal.bam
	"""

}
process Analyze_covariates {
	publishDir "${params.outdir}/BQSR", mode: 'copy'
	input:
	file(recal_table) from recal_table_ch2
	path ref from params.REF
	output:
	file("*")
	script:
	"""
	module load gatk
	java -Xmx80G -jar /opt/sw/3.4.42/GenomeAnalysisTK.jar AnalyzeCovariates -R ${ref}\
	-before ${recal_table} -after after_recal.table -plots ${pair_id}_recal_plots.pdf
	"""
}
*/
process delly_variants {
	publishDir "${params.outdir}/bcf_output", mode:'copy'
	input:
	set pair_id, path(marked_bam) from markedup_bam_ch
	path ref from params.REF
	file(baifile) from marked_index_ch
	output:
	set pair_id, file("${pair_id}_delly.bcf") into delly_bcf_ch
	script:
	"""
	delly call -o ${pair_id}_delly.bcf -g ${ref} ${marked_bam}
	
	"""
}



process bcf_to_vcf {
	publishDir "${params.outdir}/variant_calling", mode:'copy'
	input: 
	set pair_id, file(reads_bcf) from delly_bcf_ch
	output:
	set pair_id, file("${pair_id}_SV.vcf") into vcf_ch	
	script:
	"""
	module load bcftools
	bcftools view ${reads_bcf} > ${pair_id}_SV.vcf	
	"""
}
/*
process build graph {
	publishDir "${outdir}/vg_graph", mode: 'copy'
	input:
	output:
	script:
	"""
	code
	"""
}
*/
