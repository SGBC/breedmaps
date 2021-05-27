#!/bin/bash

Rscript ~/breedmaps/SVs_identification/scripts/gene_annotation/filtering_SVs.R

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir	Base directory		default="~/breedmaps/SVs_identification/"
#  -a, --annDir	Combined annotation	default="data/annotation/"
#  -s, --scriptDir	Script directory	default ="scripts/gene_annotation/"
#  -v, --vcfDir	VCF directory		default="data/vcf/"
#  -r, --resultsDir	Result directory	default="results/gene_annotation/"
#  -f, --functions	Function file name	default="functions.R"


Rscript ~/breedmaps/SVs_identification/scripts/gene_annotation/gene_ann.R

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir	Base directory			default="~/breedmaps/SVs_identification/"
#  -a, --annDir	Combined annotation		default="data/annotation/"
#  -s, --scriptDir	Script directory		default ="scripts/gene_annotation/"
#  -d, --dataDir	Dir for the filtered variants	default="results/gene_annotation/filtered_variants/"
#  -r, --resultsDir	Result directory		default="results/gene_annotation/annotated_variants"
#  -f, --functions	Function file name		default="functions.R"
#  -n, --annFile	ENSEMBL annotation file	default ="Bos_taurus.ARS-UCD1.2.103.gtf"

Rscript ~/breedmaps/SVs_identification/scripts/gene_annotation/sv_datasets_joining.R -e "~/breedmaps/SVs_identification/data/annotation/bos_taurus_structural_variations.gvf"

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir     Base directory			        default="~/breedmaps/SVs_identification/"
#  -a, --annDir         Combined annotation		        default="data/annotation/"
#  -s, --scriptDir      Script directory		        default ="scripts/gene_annotation/"
#  -d, --dataDir	    Dir for the filtered variants	default="results/gene_annotation/filtered_variants/"
#  -r, --resultsDir	    Result directory		        default="results/gene_annotation/annotated_variants"
#  -f, --functions	    Function file name		        default="functions.R"
#  -e, --datasetPath    Dataset path in gff/gvf format  default ="data/eva/remapped/"

for file in $(ls ~/breedmaps/SVs_identification/data/eva/remapped/*.gff)
do
    Rscript ~/breedmaps/SVs_identification/scripts/gene_annotation/sv_datasets_joining.R -e $file

done



