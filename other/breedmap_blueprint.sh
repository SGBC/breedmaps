#!/bin/bash

module load bwa
bwa index transcriptome.fa #indexing
bwa mem transcriptome.fa gut_1.fq gut_2.fq > out_put.sam #mapping
module load samtools
samtools view -h -b -S out_put.sam > out_put.bam # sam to bam
samtools view -b -F 4 out_put.bam > out_put.mapped.bam # filtering mapped reads
samtools sort -m 1000000000 out_put.mapped.bam > out_put.mapped.sorted.bam # sorting the aligned reads

#MarkDuplicates(if needed)
java -jar picard.jar MarkDuplicates \
      I=input.bam \
      O=marked_duplicates.bam \
	  REMOVE_DUPLICATES=true \
      M=marked_dup_metrics.txt

#Variant calling:
conda init bash
#activating conda environment
module load delly
delly call -x ref.excl(BED-file) -o delly.bcf -g transcriptome.fa input.bam # calling variants
#conda deactivate
module load bcftools
bcftools view delly.bcf > delly.vcf # converting bcf to vcf

#building graph 
module load vg/1.6.0
vg construct -r transcriptome.fa -m 32 > ref.vg # constructing a graph
vg view transcriptome.ref.vg # displaying in GFA format 
vg view -j transcriptome.ref.vg # displaying in JSON format
vg view -d transcriptome.ref.vg # displaying in DOT format
#To work with the JSON output use the tool jq. To get all sequences in the graph
vg view -j tiny.ref.vg | jq '.node[].sequence'
# graphviz to layout the graph representation in DOT format.
vg view -d tiny.ref.vg | dot -Tpdf -o tiny.ref.pdf
# Build a new graph that has some variants built into it. First, take a look at at tiny/tiny.vcf.gz, which contains variants in (gzipped) VCF format.
vg construct -r transcriptome.fa -v delly.vcf -m 32 > ref_variants_graph.vg
#Visualizing the outcome.
#add this path to the visualization.
# original ref path
vg view -dp tiny.ref.vg | dot -Tpdf -o tiny.pdf 
