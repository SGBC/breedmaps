#!/usr/bin/env nextflow

params.ref="/home/phanindra/TOY/transcriptome.fa" //Reference genome
params.outdir="/home/phanindra/TOY/result_alignment" //output directory
params.out = "${params.outdir}/out"


	log.info """\
         ============================================    
         || GRAPH GENOME - N F   P I P E L I N E   ||
         ============================================
         
         Reference Genome: ${params.ref}
         Read_file(s)    : ${params.reads}
         outdir          : ${params.outdir}
         """
         .stripIndent()

// Process to index the reference genome

process index {
	publishDir "${params.out}/indexed_ref", mode:'copy'
	//container defining.
    input:
    path transcriptome from params.ref
     
    output:
    //path 'index' into index_ch
    file("*") into index_ch

    script:       
    """
    bwa index ${transcriptome}
    """

}

params.accession = 'SRP043510' //string resembling the accession number
reads = Channel.fromSRA(params.accession) //creating channel

// Another option to create a channel 

reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

// Process to perform quality control and outputs the files necessary for analysis

process fastqc {

	publishDir "${params.out}/fastqc", mode:'copy'
	//container defining.	
	input:
	tuple sample_id, file(sample_files) from reads_ch
	output:
	file("fastqc_${sample_id}_logs") into fastqc_ch

	script:
	"""
	mkdir fastqc_${sample_id}_logs
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}	  
	"""

}

reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

// Process for mapping the paired reads to the reference and outputs the sam file. 

process align {
	publishDir "${params.out}/aligned_reads", mode:'copy'
	//container defining.
	input:
	set pair_id, file(reads) from reads_ch //the channel creates the pair ids
	path(ref) from params.ref
	output:
	set val(pair_id), file("${pair_id}_aligned_reads.sam") \
	into aligned_reads_ch
	script:
	"""
	bwa mem ${ref} ${reads} > ${pair_id}_aligned_reads.sam #mapping
	"""

}

//Process for converting the sam to bam, sorting the bam file and finally indexing the bam file.

process sam_to_bam_sorted {
	publishDir "${params.out}/aligned_reads", mode:'copy'
	//container defining.
	
	input:
	set pair_id, file(aligned_reads) from aligned_reads_ch //the channel creates the pair ids
	
	output:
	set val(pair_id),  file("${pair_id}_aligned_reads_sorted.bam"), file("${pair_id}_aligned_reads_sorted.bam.bai")  into bam_for_sorting_bam_ch
	
	script:
	"""
	samtools view -S -b ${aligned_reads} > ${pair_id}_aligned_reads.bam
    samtools sort ${pair_id}_aligned_reads.bam -o ${pair_id}_aligned_reads_sorted.bam
    samtools index ${pair_id}_aligned_reads_sorted.bam
	"""

}

// Process for creating a dictionary for reference genome which will be used in Markduplicates process.

process genome_dict {
	
	publishDir "${params.out}/index_dict", mode:'copy'
	input:
    path transcriptome from params.ref
    
    output:
    file("*.dict") into picard_dict_ch
    file("*.fai") into samtools_index_ch, gen_indx_var_filtering_ch
    
    script:
    """
    echo "java -jar /opt/Core/picard/picard.jar \
            CreateSequenceDictionary R= ${transcriptome} \
            O=${transcriptome}.dict" > ${transcriptome}.dict
            
    echo "samtools faidx ${transcriptome}" > ${transcriptome}.fai
    """
 
}

// Process for marking the duplicate reads in the bam file.

process mark_duplicates {
	publishDir "${params.out}/markedup", mode:'copy'
	container 'cbcrg/callings-nf@sha256:b65a7d721b9dd2da07d6bdd7f868b04039860f14fa514add975c59e68614c310'
	input:
    set val(pair_id), path(bamfile), path(baifile) from bam_for_sorting_bam_ch
	output:
	set val(pair_id), path("${pair_id}_sorted_duplicates_rm.bam") into bam_for_variant_calling_ch
    set val(pair_id), path("${pair_id}_sorted_duplicates_metrics.txt") into dedup_qc_ch
	script:
	"""
	echo "java -jar /opt/Core/picard/picard.jar \
        MarkDuplicates I=${bamfile} O=${pair_id}_sorted_duplicates_rm.bam \
        M=${pair_id}_sorted_duplicates_metrics.txt \
        REMOVE_SEQUENCING_DUPLICATES=true" > ${pair_id}_sorted_duplicates_rm.bam
        
    echo "${pair_id}" > ${pair_id}_sorted_duplicates_metrics.txt
	"""
}

// Process for creating a bed file for the reference genome

process bed_file{
	publishDir "${params.out}/bed_file", mode:'copy'
	//container defining
	input:
	output:
	script:
	"""
	"""

}

// Process for calling variants through delly

process delly_variants{
	publishDir "${params.out}/variant_calling", mode:'copy'
	//container defining
	input:
	set value(pair_id), path(MD_bamfile)from bam_for_variant_calling_ch
	path transcriptome from params.ref
	output:
	set value(pair_id), path("${pair_id}_delly.bcf") into delly_ch
	script:
	"""
	delly call -x ref.excl(BED-file) -o {pair_id}_delly.bcf -g ${transcriptome} ${MD_bamfile} # calling variants
	"""
}

process bcf_to_vcf{
	publishDir "${params.out}/variant_calling", mode:'copy'
	//container defining
	input:
	set value(pair_id), path(bcf_file) from delly_ch
	output:
	set value(pair_id), path("${pair_id}_delly.vcf") into delly_vcf_ch
	script:
	"""
	bcftools view ${pair_id}_delly.bcf > ${pair_id}_delly.vcf
	"""
}

process vg_graph{

	publishDir "${params.out}/graph_genome", mode:'copy'
	//container defining
	input:
	path transcriptome from params.ref
	output:
	file("*") into vg_ch
	script:
	"""
	vg construct -r ${transcriptome} -m 32 > ref.vg # constructing a graph
	"""


}
