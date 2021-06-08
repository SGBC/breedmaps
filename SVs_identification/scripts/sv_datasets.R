# Overlapping with existing annotation

if (!require("tidyverse")) {
  install.packages("tidyverse")
}

if (!require("BiocManager")) {
  install.packages("BiocManager")
  BiocManager::install()
  library(BiocManager)
}
if (!require("GenomicRanges")) {
  BiocManager::install("GenomicRanges")
}
if (!require("rtracklayer")) {
  BiocManager::install("rtracklayer")
}

if (!require("optparse")) {
  install.packages("optparse")
}
library(optparse)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)


# Create inputs from the command line and save the options
options <- list(
  make_option(c("-w", "--workingDir"), help = "Base directory", default =
                "~/breedmaps/SVs_identification/"),
  make_option(c("-a", "--annDir"), help = "Annotation directory", default =
                "data/annotation/"),
  make_option(c("-e", "--datasetDir"), help = "Dataset directory", default =
                "data/eva/remapped/"),
  make_option(c("-k", "--datasetName"), help = "Dataset name", default =
                "remapped_nstd119_Menzi_et_al_2016.2017-04-24.Bos_taurus_UMD_3.1.1.Submitted.gff"),
  make_option(c("-s", "--scriptDir"), help = "script directory", default =
                "scripts/"),
  make_option(c("-d", "--dataDir"), help = "Data directory for the filtered variants", default =
                "results/filtered_variants/"),
  make_option(c("-r", "--resultsDir"), help = "Result directory", default =
                "results/datasets/"),
  make_option(c("-f", "--functions"), help = "Function file name", default =
                "functions.R")
)


params <- parse_args(OptionParser(option_list = options))
source(paste(params$workingDir, params$scriptDir, params$functions, sep =
               ""))

######################################################################
######################################################################
######################################################################
dataset_path = paste(params$workingDir, params$datasetDir, params$datasetName,sep="")

dataset = as.data.frame(rtracklayer::import(dataset_path))
dataset_cleaned = dataset %>% dplyr::rename(
  DATASET_CHROM = seqnames,
  DATASET_SOURCE = source,
  DATASET_SV_TYPE = type,
  DATASET_START = start,
  DATASET_END = end,
  DATASET_STRAND = strand,
  DATASET_SV_LENGTH = width
)
dataset_unique = dataset_cleaned %>% dplyr::select(DATASET_CHROM, DATASET_START, DATASET_END, DATASET_STRAND, DATASET_SV_TYPE) %>% unique()
# Creates ranges from the dataset using package GenomicRanges
dataset_range = makeGRangesFromDataFrame(
  dataset_unique,
  seqnames.field = "DATASET_CHROM",
  start.field = "DATASET_START",
  end.field = "DATASET_END",
  strand.field = "DATASET_STRAND"
)

######################################################################
######################################################################
######################################################################
# Load the filtered SVs outputted from filtering_SVs.R
data_path = paste(params$workingDir,
                  params$dataDir,
                  sep = "")
data_file_names = list.files(path = data_path, pattern = "precise_*")

overlap_list = list()
for (i in 1:(length(data_file_names))) {
  full_path = paste(params$workingDir,
                    params$dataDir,
                    data_file_names[[i]],
                    sep = "")
  svs = read.table(
    file = full_path,
    quote = "",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  ) %>% dplyr::rename(
    SV_CHROM = CHROM,
    SV_START = POS,
    SV_END = END,
    ID1 = ID
  )
  sv_range = makeGRangesFromDataFrame(
    svs,
    seqnames.field = "SV_CHROM",
    start.field = "SV_START",
    end.field = "SV_END",
  )
  dataset_cleaned$ID2 = rownames(dataset_cleaned)
  overlap = findoverlap_dataframe(sv_range, dataset_range, svs, dataset_cleaned)
  #sv_range = ID1, dataset = ID2
  overlap_list[[i]] = overlap %>% dplyr::rename(SV_ID = ID1) %>% dplyr::select(-ID2)
  filt_results = paste(
    params$workingDir,
    params$resultsDir,
    params$datasetName,
    "_",
    data_file_names[[i]],
    sep = ""
  )
  write.table(
    x = overlap_list[[i]],
    file = filt_results,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
}




