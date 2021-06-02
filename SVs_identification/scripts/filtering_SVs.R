# FILTERING SVs

# Check if the required packages are installed
if (!require("tidyverse")) {
  install.packages("tidyverse")
}
if (!require("optparse")) {
  install.packages("optparse")
}
library(optparse)
library(tidyverse)
library(ggplot2)

################################################################################
################################################################################
################################################################################

# Create inputs from the command line and save the options
options <- list(
  make_option(c("-w", "--workingDir"), help = "Base directory", default =
                "~/breedmaps/SVs_identification/"),
  make_option(c("-a", "--annDir"), help = "Combined annotation", default =
                "data/annotation/"),
  make_option(c("-s", "--scriptDir"), help = "Script directory", default =
                "scripts/"),
  make_option(c("-v", "--vcfDir"), help = "VCF directory", default = "data/vcf/"),
  make_option(c("-r", "--resultsDir"), help = "Result directory", default =
                "results/filtered_variants/"),
  make_option(c("-f", "--functions"), help = "Function file name", default =
                "functions.R")
)

params <- parse_args(OptionParser(option_list = options))

source(paste(params$workingDir, params$scriptDir, params$functions, sep =
               ""))

################################################################################
################################################################################
################################################################################

# Load RDC VCF files
vcf_path = paste(params$workingDir, params$vcfDir, sep = "")
vcf_files = list.files(path = vcf_path, pattern = "*.vcf")
vcf_names = list()

for (i in 1:length(vcf_files)) {
  # Paste together the entire file path
  vcf_file = paste(vcf_path, vcf_files[i], sep = "")
  
  # Save the file name to be able to use it when writing results files
  name = strsplit(vcf_files, "\\.")[[i]][1]
  vcf_names[[i]] = name
  
  # Read vcf file
  df = read_new_delly_vcf(vcf_file)
  
  # Filter out the IMPRECISE and LowQual SVs
  filt_var = df %>% dplyr::filter(FILTER == "PASS") %>% filter(PRECISE == TRUE)
  filt_var = filt_var %>% dplyr::mutate(SV_LENGTH = as.integer(END) - as.integer(POS))
  
  # Write results to files
  filt_results = paste(
    params$workingDir,
    params$resultsDir,
    "precise_",
    vcf_names[[i]],
    ".tsv",
    sep = ""
  )
  write.table(
    x = filt_var,
    file = filt_results,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
  png(
    filename = paste(
      params$workingDir,
      params$resultsDir,
      "hist_",
      vcf_names[[i]],
      ".png",
      sep = ""
    )
  )
  hist(filt_var$SV_LENGTH,main=paste(vcf_names[[i]]," SV_LENGTH",sep=""), xlab="ID", ylab="LENGTH")
  dev.off()

}



