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
  make_option(c("-b", "--workingDir"), help = "Base directory", default =
                "~/breedmaps/SVs_identification/"),
  make_option(c("-a", "--annDir"), help = "Annotation directory", default =
                "data/annotation/"),
  make_option(c("-n", "--svsDir"), help = "SVs gff file from ENSEMBL", default =
                "Bos_taurus.ARS-UCD1.2.103.gtf"),
  make_option(c("-e", "--evaFile"), help = "EVA gff file", default =
                "data/eva/remapped/x-special/remapped_estd223_Boussaha_et_al_2015.2015-11-02.Bos_taurus_UMD_3.1.Submitted.gff"),
  make_option(c("-s", "--scriptDir"), help = "script directory", default =
                "scripts/gene_annotation/"),
  make_option(c("-d", "--dataDir"), help = "Data directory for the filtered variants", default =
                "results/gene_annotation/filtered_variants/"),
  make_option(c("-r", "--resultsDir"), help = "Result directory", default =
                "results/gene_annotation/ensembl_ann/"),
  make_option(c("-f", "--functions"), help = "Function file name", default =
                "functions.R")
)


params <- parse_args(OptionParser(option_list = options))
source(paste(params$workingDir, params$scriptDir, params$functions, sep =
               ""))

######################################################################
######################################################################
######################################################################

dataset_path = paste(params$workingDir,
                     params$evaFile,
                     sep = "")

dataset = as.data.frame(rtracklayer::import(dataset_path))
dataset_renamed = dataset %>% dplyr::rename(
  DATASET_CHROM = seqnames,
  DATASET_SOURCE = source,
  DATASET_SV_TYPE = type,
  DATASET_START = start,
  DATASET_END = end,
  DATSET_ID = ID,
  DATSET_STRAND = strand,
  DATASET_SV_LENGTH = width
)

dataset_range = makeGRangesFromDataFrame(
  dataset_renamed,
  seqnames.field = "DATASET_CHROM",
  start.field = "DATASET_START",
  end.field = "DATASET_END",
  strand.field = "DATASET_STRAND"
)

######################################################################
######################################################################
######################################################################

data_path = paste(params$workingDir,
                  params$resultsDir,
                  sep = "")
data_file_names = list.files(path = data_path, pattern = "precise_[A-Z]*")

overlap_list = list()
for (i in 1:(length(gene_file_names))) {
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
  )
  %>% dplyr::rename(SV_CHROM = CHROM,
                    SV_START = POS,
                    SV_END = END)
  svs = svs %>% dplyr::filter(SV_LENGTH < 50000000) %>% dplyr::rename(ID1 = ID)
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
    "dataset/",
    "dataset_",
    gene_file_names[[i]],
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