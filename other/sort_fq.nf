#!/usr/bin/env nextflow
params.outdir = "/proj/breedmap/sorted_data"
reads_ch = Channel.fromFilePairs('/proj/breedmap/data/*{1,2}.fastq.gz', flat: true)

process sort_fq {
        publishDir "${params.outdir}", mode:'copy'

        input:
        set sample_id, file(sam1), file(sam2) from reads_ch

        output:
        file("*")

        script:
        """
        zcat ${sam1} |paste - - - - |sort -k1,1 -S 3G| tr '\t' '\n'| gzip > ${sample_id}sorted_1.fastq.gz
        zcat ${sam2} |paste - - - - |sort -k1,1 -S 3G| tr '\t' '\n'| gzip > ${sample_id}sorted_2.fastq.gz
        """

}
