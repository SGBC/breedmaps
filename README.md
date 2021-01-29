# breedmaps
## Variant calling and creating a graph genome
Pipeline for variant calling and construction of graph genome.
Scripts: bash, Nextflow.

Tools 

Mapping :        1. BWA-MEM
                 2. VGtools
                 
Variant calling: 1. Delly 
                 2. GATK
                 3. VGtools

Graph Genome :   1. VGtools

Samtools, bcftools and BEDtools were used for additional tasks.

Before running main.nf perform trimming and analysis with fastp using the fastp.nf.
You can parallely run all the available files for trimming based on the available memory.

## Analysis of the structural variation
The analysis of the structural variations(VCF files) can be found in the folder SVs_identification. The analysis was done by Jenny Jakobsson.
