# FILTERING SVs


if(!require("tidyverse")){
  install.packages("tidyverse")
}
if(!require("optparse")){
  install.packages("optparse")
}
library(optparse)
library(tidyverse)


# Create inputs from the command line and save the options
options <- list( 
  make_option(c("-b", "--workingDir"), help="Base directory", default="~/breedmaps/SVs_identification/"),
  make_option(c("-n", "--annDir"), help="Combined annotation", default="data/annotation/"),
  make_option(c("-s", "--scriptDir"), help="script directory", default ="scripts/gene_annotation/"),
  make_option(c("-v", "--vcfDir"), help="VCF directory", default ="data/vcf/"),
  make_option(c("-r", "--resultsDir"), help="Result directory", default="results/gene_annotation/"),
  make_option(c("-f", "--functions"), help="Function file name",default="functions.R") )

params <- parse_args(OptionParser(option_list=options))


source(paste(params$workingDir, params$scriptDir, params$functions, sep=""))

########################################
########################################
########################################
# Load RDC VCF files
vcf_path = paste(params$workingDir, params$vcfDir, sep = "")
vcf_files = list.files(path = vcf_path, pattern = "*.vcf")
vcf_names = list()
df_list = list()
all_bta_df = data.frame()

for (i in 1:length(vcf_files)) {
  vcf_file = paste(vcf_path, vcf_files[i], sep = "")
  
  # Save the file name to be able to use it when writing results files
  name = strsplit(vcf_files, "\\.")[[i]][1]
  vcf_names[[i]] = name
  
  #df = read_delly_vcf(vcf_file)
  df = read_new_delly_vcf(vcf_file)
  
  # Because of the ranges need to filter out the breakpoint
  filt_var = df %>% dplyr::filter(FILTER == "PASS") %>% filter(PRECISE == TRUE)
  filt_var_final = filt_var %>% dplyr::mutate(SV_LENGTH = as.numeric(END) -
                                                as.numeric(POS))
  df_list[[i]] = filt_var_final
  if (i == 1){
    filt_var_final$FILE = vcf_files[i]
    all_df = filt_var_final
  }
  else {
    filt_var_final$FILE = vcf_files[i]
    all_df = full_join(all_df, filt_var_final)
  }
}

########################################
########################################
########################################
for (i in 1:length(df_list)) {
  filt_results = paste(
    params$workingDir,
    params$resultsDir,
    "filtered_variants/",
    "precise_",
    vcf_names[[i]],
    ".tsv",
    sep = ""
  )
  write.table(
    x = df_list[[i]],
    file = filt_results,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
}

all_results = paste(
  params$workingDir,
  params$resultsDir,
  "filtered_variants/",
  "all_precise.tsv",
  sep = ""
)
write.table(
  x = all_df,
  file = all_results,
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)




