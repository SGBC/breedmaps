#!/bin/bash

# Filters out LowQual and IMPRECISE SVs 
Rscript ~/breedmaps/SVs_identification/scripts/filtering_SVs.R

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir	Base directory		default="~/breedmaps/SVs_identification/"
#  -a, --annDir	Combined annotation	default="data/annotation/"
#  -s, --scriptDir	Script directory	default ="scripts/"
#  -v, --vcfDir	VCF directory		default="data/vcf/"
#  -r, --resultsDir	Result directory	default="results/"
#  -f, --functions	Function file name	default="functions.R"

# Overlaps SVs from filtering_SVs.R with annotation downloaded from Ensembl
# Output is stored in specified folder (or default results/annotated_variants/)
Rscript ~/breedmaps/SVs_identification/scripts/ensembl_ann.R

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir	Base directory			default="~/breedmaps/SVs_identification/"
#  -a, --annDir	Combined annotation		default="data/annotation/"
#  -s, --scriptDir	Script directory		default ="scripts/"
#  -d, --dataDir	Dir for the filtered variants	default="results/filtered_variants/"
#  -r, --resultsDir	Result directory		default="results/annotated_variants"
#  -f, --functions	Function file name		default="functions.R"
#  -n, --annFile	ENSEMBL annotation file	default ="Bos_taurus.ARS-UCD1.2.103.gtf"

# Overlaps SVs from filtering_SVs.R with regulatory regions from file taken from Fang et al.
# Output is stored in specified folder (or default results/regulatory_variants/)
Rscript ~/breedmaps/SVs_identification/scripts/regulatory_regions.R
# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir     Base directory			        default="~/breedmaps/SVs_identification/"
#  -a, --annDir         Combined annotation		        default="data/annotation/"
#  -s, --scriptDir      Script directory		        default ="scripts/"
#  -d, --dataDir	      Dir for the filtered variants	default="results/filtered_variants/"
#  -r, --resultsDir	    Result directory		        default="results/regulatory_variants/"
#  -f, --functions	    Function file name		        default="functions.R"
#  -e, --datasetPath    Dataset path in gff/gvf format  default ="data/eva/remapped/"

# Overlaps SVs from filtering_SVs.R with SVs registered in Ensembl (downloaded file)
# Output is stored in specified folder (or default results/datasets/)
Rscript ~/breedmaps/SVs_identification/scripts/sv_datasets.R -e "data/annotation/" -k "bos_taurus_structural_variations.gvf"

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir     Base directory			        default="~/breedmaps/SVs_identification/"
#  -a, --annDir         Combined annotation		        default="data/annotation/"
#  -s, --scriptDir      Script directory		        default ="scripts/"
#  -d, --dataDir	    Dir for the filtered variants	default="results/datasets/"
#  -r, --resultsDir	    Result directory		        default="results/annotated_variants"
#  -f, --functions	    Function file name		        default="functions.R"
#  -e, --datasetPath    Dataset directory path in gff/gvf format  default ="data/eva/remapped/"
#  -k --datasetName     Dataset file                default ="remapped_nstd119_Menzi_et_al_2016.2017-04-24.Bos_taurus_UMD_3.1.1.Submitted.gff"

# Overlaps all datasets downloaded from EVA with same script used above
# Each loaded SV file generates one output per dataset
for file in $(ls -h ~/breedmaps/SVs_identification/data/eva/remapped/)
do
    Rscript ~/breedmaps/SVs_identification/scripts/sv_datasets.R -k $file

done

# Counting output for each overlap per filtered SV
# This is a computational heavy script if there are a lot of SVs involved
# Input from this script is the output from previous all script(filtering_SVs.R, ensembl_ann.R, regulatory_regions.R, sv_datasets for 
# ensembl SVs and all 6 datsets)
for file in $(ls -h ~/breedmaps/SVs_identification/data/vcf/)
do
    Rscript ~/breedmaps/SVs_identification/scripts/counting_SV_occurence.R -n $file

done

#####################################################
# Extra analysis using manual input steps
#####################################################

Rscript ~/breedmaps/SVs_identification/scripts/vep_filtering.R

# Add flag to specify if the files in a different directory than the default
#  -w, --workingDir     Base directory			        default="~/breedmaps/SVs_identification/"
#  -d, --dataDir	    Dir for the filtered variants	default="results/filtered_variants/"


cd ~/breedmaps/SVs_identification/results/summery_count/
# The input for the python script are the vcf and the filtered SVs from filtering_SVs.R
# The file structure is relative to the home folder but can be changed in the scripts individually
python3 ~/breedmaps/SVs_identification/scripts/summery_count/sum_BTA_single_samples.py
python3 ~/breedmaps/SVs_identification/scripts/summery_count/sum_BTA_combined_samples.py
python3 ~/breedmaps/SVs_identification/scripts/summery_count/sum_RDC.py
python3 ~/breedmaps/SVs_identification/scripts/summery_count/sum_BTA_comb_RDC.py

gzip ~/breedmaps/SVs_identification/results/filtered_variants/*.tsv
gzip ~/breedmaps/SVs_identification/results/annotated_variants/*.tsv
gzip ~/breedmaps/SVs_identification/results/regulatory_variants/*.tsv
gzip ~/breedmaps/SVs_identification/results/datasets/*.tsv


