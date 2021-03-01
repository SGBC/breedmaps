# breedmaps
## Variant calling and creating a graph genome
- Pipeline for variant calling and construction of graph genome.
- Scripts: bash, Nextflow.
Key role of this pipeline is to call large structural variants. Various tools were used based on different SV strategies.
- BWA is used for indexing and mapping the reads against the reference genome (ARS-UCD1.2).
- All the parameters are set to default while using the following tools.
- "Nextflow run main.nf" to run the pipeline and "Nextflow run main.nf -resume" to resume the process.

## Analysis of the structural variation
The analysis of the structural variations(VCF files) can be found in the folder SVs_identification. The analysis was done by Jenny Jakobsson.
