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
  make_option(c("-s", "--scriptDir"), help = "script directory", default =
                "scripts/gene_annotation/"),
  make_option(c("-d", "--dataDir"), help = "Data directory for the filtered variants", default =
                "results/gene_annotation/filtered_variants/"),
  make_option(c("-r", "--resultsDir"), help = "Result directory", default =
                "results/gene_annotation/ensembl_ann/"),
  make_option(c("-f", "--functions"), help = "Function file name", default =
                "functions.R"),
  make_option(c("-e", "--ensembl_svs"), help = "Ensembl annotation file name", default =
                "bos_taurus_structural_variations.gvf")
)


params <- parse_args(OptionParser(option_list = options))
source(paste(params$workingDir, params$scriptDir, params$functions, sep =
               ""))



######################################################################
######################################################################
######################################################################

path = paste(params$workingDir,
             params$annDir,
             params$ensembl_svs,
             sep = "")
ensembl_svs = read.table(
  file = path,
  quote = "",
  sep = "\t",
  header = F,
  stringsAsFactors = F
)
ensembl_svs = ensembl_svs %>% dplyr::rename(
  ENSEMBL_CHROM = V1,
  ENSEMBL_SOURCE = V2,
  ENSEMBL_SV_TYPE = V3,
  ENSEMBL_START = V4,
  ENSEMBL_END = V5,
  Empty1 = V6,
  ENSEMBL_STRAND = V7,
  Empty2 = V8,
  ENSEMBL_INFO = V9
)
ensembl_cleaned = ensembl_svs %>% dplyr::select(ENSEMBL_CHROM, ENSEMBL_START, ENSEMBL_END, ENSEMBL_STRAND, ENSEMBL_SOURCE, ENSEMBL_INFO)

ensembl_range = makeGRangesFromDataFrame(
  ensembl_cleaned,
  seqnames.field = "ENSEMBL_CHROM",
  start.field = "ENSEMBL_START",
  end.field = "ENSEMBL_END",
  strand.field = "ENSEMBL_STRAND"
)

path = paste(params$workingDir,
             params$resultsDir,
             sep = "")
gene_file_names = list.files(path = path, pattern = "precise_[A-Z]*")

overlap_list = list()
for (i in 1:(length(gene_file_names))) {
  full_path = paste(params$workingDir,
                    params$resultsDir,
                    gene_file_names[[i]],
                    sep = "")
  svs = read.table(
    file = full_path,
    quote = "",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  ) %>% dplyr::rename(SV_CHROM = CHROM,
                      SV_START = POS,
                      SV_END = END)
  sv_range = makeGRangesFromDataFrame(
    svs,
    seqnames.field = "SV_CHROM",
    start.field = "SV_START",
    end.field = "SV_END",
  )
  ensembl_cleaned$ID2 = rownames(ensembl_cleaned)
  overlap = findoverlap_dataframe(sv_range, ensembl_range, svs, ensembl_cleaned)
  #sv_range = ID1, ensembl = ID2
  overlap_list[[i]] = overlap %>% dplyr::rename(SV_ID = ID1) %>% dplyr::select(-ID2)
  filt_results = paste(
    params$workingDir,
    "results/",
    "gene_annotation/",
    "ensembl_ann/",
    "ensembl_",
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