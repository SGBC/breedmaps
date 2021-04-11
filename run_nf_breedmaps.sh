#!/bin/bash

#$ -N phan_mapping
#$ -M raghumandala.phanindra.balaji@slu.se
#$ -m seab
#$ -cwd 
#$ -l h_rt=24:0:0,h_vmem=5G
#$ -j y
#$ -pe smp 24
#$ -e phan_mapping.log
#$ -o phan_mapping.log

echo "Starting the pipeline"
module load conda 
. ../../scripts/conda.sh
conda activate /home/phanindra/env/SVpipeline
../../scripts/nextflow run ../../scripts/main.nf
