#!/bin/bash
# SURVIVOR was downloaded from github: https://github.com/fritzsedlazeck/SURVIVOR
#ls -d ~/breedmaps/SVs_identification/data/vcf/* > sample_files

#/Users/jj/SURVIVOR/Debug/SURVIVOR merge sample_files 1000 2 1 1 0 15 merged.vcf

# File with VCF names and paths
# max distance between breakpoints (0-1 percent of length, 1- number of bp) 
# Minimum number of supporting caller
# Take the type into account (1==yes, else no)
# Take the strands of SVs into account (1==yes, else no)
# Disabled.
# Minimum size of SVs to be taken into account.
# Output VCF filename

#/Users/jj/SURVIVOR/Debug/SURVIVOR stats merged.vcf -1 -1 -1 merged.stats
# vcf file
# Min SV size (disable: -1)
# Max SV size (disable: -1)
# Min number read support (disable: -1)
# output summary file
cd ~/breedmaps/SVs_identification/data/vcf/
for file in $(ls ~/breedmaps/SVs_identification/data/vcf/*.vcf)
do
    /Users/jj/SURVIVOR/Debug/SURVIVOR stats $file -1 -1 -1 "$file".stats
done

mv ~/breedmaps/SVs_identification/data/vcf/BTA1*.stats* ~/breedmaps/SVs_identification/results/SURVIVOR/