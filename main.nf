#!/usr/bin/env/ nextflow
params.ref = "/home/phanindra/TOY/transcriptome.fa" // change the path accordingly
reads_ch = channel.fromFilePairs('/home/phanindra/TOY/*_{1,2}.fq') // creating the channel 
params.outdir = "/home/phanindra/results" //A directory for storing the results

process indexing_aligning {
        publishDir "${params.outdir}/indexed_aligned", mode: 'copy' // saves all the files available in the output into the specified directory
        input :
        path ref from params.ref
        tuple pair_id, file(reads) from reads_ch
        output:
        file("*")
        set val(pair_id),file("${pair_id}_aligned_reads.sam") into aligned_ch //creating a new channel for the next process
        script:
        """
        module load bwa 
        bwa index ${ref} // index the ref genome
        bwa mem ${ref} ${reads} > ${pair_id}_aligned_reads.sam // mapping the reads against the ref genome

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
        samtools view -h -b -S ${reads} > ${pair_id}_aligned_reads.bam // converting reads to bam format
        samtools sort ${pair_id}_aligned_reads.bam -o ${pair_id}_aligned_reads_sorted.bam // sorting the reads
        samtools index ${pair_id}_aligned_reads_sorted.bam // indexing the bam reads
        """
}

process markedups {
        publishDir "${params.outdir}/markedup", mode:'copy'
        input:
        set val(pair_id), file(bamfile), file(baifile) from bam_ch
        output:
        set pair_id, file("${pair_id}_sorted_duplicates_rm.bam") into markedup_bam_ch
        set pair_id, file("${pair_id}_markedup_metrics.txt") into marked_metrics_ch
        file("${pair_id}_sorted_duplicates_rm.bam.bai") into marked_index_ch
        script:
        """
        java -jar /opt/sw/picard/1.137/picard.jar MarkDuplicates \
        I=${bamfile} O=${pair_id}_sorted_duplicates_rm.bam M=${pair_id}_markedup_metrics.txt // marking the duplicates
        module load samtools
        samtools index ${pair_id}_sorted_duplicates_rm.bam // indexing the updated bam

        """
}

process delly_variants {
        publishDir "${params.outdir}/variant_calling", mode:'copy'
        input:
        set pair_id, path(marked_bam) from markedup_bam_ch
        path ref from params.ref
        file(baifile) from marked_index_ch
        output:
        set pair_id, file("${pair_id}_delly.bcf") into delly_bcf_ch
        script:
        """
        delly call -o ${pair_id}_delly.bcf -g ${ref} ${marked_bam} // calling the variants
        """

}
