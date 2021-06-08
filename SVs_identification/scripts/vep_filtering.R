# Filtering web-based variant effect predictor (VEP) output
# Upload the VCF file to https://www.ensembl.org/Tools/VEP, before running script
# Download the result in txt format
# The output will be printed in the terminal, no output file

# Note! The SVs will not be filtered here (as in filtering_SVs.R)

# Check if the required packages are installed
if (!require("tidyverse")) {
  install.packages("tidyverse")
}
if (!require("optparse")) {
  install.packages("optparse")
}
library(optparse)
library(tidyverse)

################################################################################
################################################################################
################################################################################

# Create inputs from the command line and save the options
options <- list(
  make_option(c("-w", "--workingDir"), help = "Base directory", default =
                "~/breedmaps/SVs_identification/"),
  make_option(c("-d", "--dataDir"), help = "Data directory for VEP output", default =
                "results/variant_effect_predictor/"),
  make_option(c("-r", "--resultsDir"), help = "Result directory", default =
                "results/"),
  make_option(c("-f", "--functions"), help = "Function file name", default =
                "functions.R")
)

params <- parse_args(OptionParser(option_list = options))

################################################################################
################################################################################
################################################################################

# Load RDC VCF files
path = paste(params$workingDir, params$dataDir, sep = "")
files = list.files(path = path, pattern = "*.txt")
print("FILE             #SVs HIGH IMPACT")

for (i in 1:length(files)) {
  # Paste together the entire file path
  full_path = paste(path, files[i], sep = "")
  
  # Read vcf file
  df = read.table(
    file = full_path,
    quote = "",
    sep = "\t",
    header = F,
    stringsAsFactors = F
  )
  df_renamed = df %>% dplyr::rename(
    Uploaded_variation = V1,
    Location = V2,
    Allele = V3,
    Consequence = V4,
    IMPACT = V5,
    SYMBOL = V6,
    Gene = V7,
    Feature_type = V8,
    Feature	= V9,
    BIOTYPE = V10,
    EXON = V11,
    INTRON = V12,
    HGVSc = V13,
    HGVSp = V14,
    cDNA_position = V15,
    CDS_position = V16,
    Protein_position = V17,
    Amino_acids = V18,
    Codons = V19,
    Existing_variation = V20,
    DISTANCE = V21,
    STRAND = V22,
    FLAGS = V23,
    SYMBOL_SOURCE = V24,
    HGNC_ID = V25,
    MANE_SELECT = V26,
    MANE_PLUS_CLINICAL = V27,
    TSL = V28,
    APPRIS = V29,
    SWISSPROT = V30,
    TREMBL = V31,
    UNIPARC = V32,
    UNIPROT_ISOFORM = V33,
    SIFT = V34,
    CLIN_SIG = V35,
    SOMATIC = V36,
    PHENO = V37,
    PHENOTYPES = V38
  )
  df_cleaned = df_renamed %>% dplyr::select(Uploaded_variation, Location, Allele, Consequence, IMPACT)
  df_high_impact = df_cleaned %>% dplyr::filter(IMPACT == "HIGH") %>% unique()
  print(paste(files[i], "             ", dim(df_high_impact)[1], sep = ""))
  
}
