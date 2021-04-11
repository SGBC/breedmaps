# breedmaps
## Variant calling
## About the pipeline
Various tools were used based on different SV strategies. Currently, only Delly is incorporated into the Nextflow pipeline for SV calling.
- Reads are trimmed using fastp with default parameters (Quality control step).
- FastQC is used to double check the quality of the reads.
- Reference genome(ARS-UCD1.2) is indexed using BWA.
- BWA-MEM is used for mapping the reads against the reference genome (ARS-UCD1.2).
- Produced SAM files are then converted to BAM files using the SAMTOOLS.
- Duplicated reads are then marked for the BAM file using PICARD software.
- Finally Delly software is used calling the SV's.
## How to run the pipeline
Before executing the pipeline there are three things that need to be fixed.
1. Path for the samples(R1,R2).
2. Path for the Reference genome(indexed reference).
3. Path for Storing the results.

- After updating the paths/location of the files accordingly.
### To run the pipeline on your local computer
`./nextflow run main.nf`
### To resume a run just add "-resume" flag 
`./nextflow run main.nf -resume`

### To run on the server
for example these scripts were ran on the planetsmasher server using QSUB job scheduler.
Details are available in the "run_nf_breedmap.sh" file, change the parameters to your need. 
And then

`qsub run_nf_breedmap.sh`

To check the status 

`qstat`

To kill a submission/run 

`qdel "JOB-ID"`

JOB-ID can be found using qstat 

All the results can be found in the specified OUTPUT folder. 



## Analysis of the structural variation
The analysis of the structural variations(VCF files) can be found in the folder SVs_identification. The analysis was done by Jenny Jakobsson.
