#!/bin/bash

Rscript ~/breedmaps/SVs_identification/scripts/filtering_SVs.R

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir	Base directory		default="~/breedmaps/SVs_identification/"
#  -a, --annDir	Combined annotation	default="data/annotation/"
#  -s, --scriptDir	Script directory	default ="scripts/"
#  -v, --vcfDir	VCF directory		default="data/vcf/"
#  -r, --resultsDir	Result directory	default="results/"
#  -f, --functions	Function file name	default="functions.R"


Rscript ~/breedmaps/SVs_identification/scripts/ensembl_ann.R

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir	Base directory			default="~/breedmaps/SVs_identification/"
#  -a, --annDir	Combined annotation		default="data/annotation/"
#  -s, --scriptDir	Script directory		default ="scripts/"
#  -d, --dataDir	Dir for the filtered variants	default="results/filtered_variants/"
#  -r, --resultsDir	Result directory		default="results/annotated_variants"
#  -f, --functions	Function file name		default="functions.R"
#  -n, --annFile	ENSEMBL annotation file	default ="Bos_taurus.ARS-UCD1.2.103.gtf"

Rscript ~/breedmaps/SVs_identification/scripts/regulatory_regions.R
# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir     Base directory			        default="~/breedmaps/SVs_identification/"
#  -a, --annDir         Combined annotation		        default="data/annotation/"
#  -s, --scriptDir      Script directory		        default ="scripts/"
#  -d, --dataDir	    Dir for the filtered variants	default="results/filtered_variants/"
#  -r, --resultsDir	    Result directory		        default="results/datasets"
#  -f, --functions	    Function file name		        default="functions.R"
#  -e, --datasetPath    Dataset path in gff/gvf format  default ="data/eva/remapped/"

Rscript ~/breedmaps/SVs_identification/scripts/sv_datasets.R -e "~/breedmaps/SVs_identification/data/annotation/bos_taurus_structural_variations.gvf"

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir     Base directory			        default="~/breedmaps/SVs_identification/"
#  -a, --annDir         Combined annotation		        default="data/annotation/"
#  -s, --scriptDir      Script directory		        default ="scripts/"
#  -d, --dataDir	    Dir for the filtered variants	default="results/datasets/"
#  -r, --resultsDir	    Result directory		        default="results/annotated_variants"
#  -f, --functions	    Function file name		        default="functions.R"
#  -e, --datasetPath    Dataset path in gff/gvf format  default ="data/eva/remapped/"


for file in $(ls ~/breedmaps/SVs_identification/data/eva/remapped/*.gff)
do
    Rscript ~/breedmaps/SVs_identification/scripts/sv_datasets.R -e $file

done


Rscript ~/breedmaps/SVs_identification/scripts/vep_filtering.R

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir     Base directory			        default="~/breedmaps/SVs_identification/"
#  -d, --dataDir	    Dir for the filtered variants	default="results/filtered_variants/"


cd ~/breedmaps/SVs_identification/results/vcf_analysis/
python3 ~/breedmaps/SVs_identification/scripts/vcf_analysis/sum_BTA_vcf.py
python3 ~/breedmaps/SVs_identification/scripts/vcf_analysis/sum_RDC_vcf.py
python3 ~/breedmaps/SVs_identification/scripts/vcf_analysis/sum_BTA_RDC.py

gzip ~/breedmaps/SVs_identification/results/filtered_variants/*.tsv
gzip ~/breedmaps/SVs_identification/results/annotated_variants/*.tsv
gzip ~/breedmaps/SVs_identification/results/regulatory_variants/*.tsv
gzip ~/breedmaps/SVs_identification/results/datasets/*.tsv


