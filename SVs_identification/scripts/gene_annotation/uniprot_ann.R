# Adding Uniprot annotation

if (!require("tidyverse")) {
  install.packages("tidyverse")
}

if (!require("BiocManager")) {
  install.packages("BiocManager")
  BiocManager::install()
  library(BiocManager)
}
if (!require("UniProt.ws")) {
  BiocManager::install("UniProt.ws")
}
if (!require("optparse")) {
  install.packages("optparse")
}
library(optparse)
library(tidyverse)
library(GenomicRanges)


# Create inputs from the command line and save the options
options <- list(
  make_option(c("-b", "--workingDir"), help = "Base directory", default =
                "~/breedmaps/SVs_identification/"),
  make_option(c("-a", "--annDir"), help = "Annotation directory", default =
                "data/annotation/"),
  make_option(c("-n", "--uniprotDir"), help = "Result directory for uniprot results", default =
                "results/gene_annotation/uniprot_ann/"),
  make_option(c("-s", "--scriptDir"), help = "script directory", default =
                "scripts/gene_annotation/"),
  make_option(c("-r", "--resultsDir"), help = "Directory for the result folder", default =
                "results/gene_annotation/annotated_variants/"),
  make_option(c("-f", "--functions"), help = "Function file name", default =
                "functions.R")
)

params <- parse_args(OptionParser(option_list = options))
source(paste(params$workingDir, params$scriptDir, params$functions, sep =
               ""))

#############################################
#############################################
#############################################

# TaxID 9913 = Bos taurus
up <- UniProt.ws(taxId = 9913)
up
cols = columns(up)
all_cols = c(
  cols[28],
  cols[31],
  cols[38],
  cols[43],
  cols[45],
  cols[49],
  cols[56],
  cols[63],
  cols[66],
  cols[69],
  cols[72],
  cols[86],
  cols[88],
  cols[102],
  cols[104:105],
  cols[117]
)
#############################################
#############################################
#############################################

gene_path = paste(params$workingDir, params$resultsDir, sep = "")
gene_file_names = list.files(path = gene_path, pattern = "filt_ann_gene_BTA126*")
joined_list = list()

for (i in 1:length(gene_file_names)) {
  full_path = paste(gene_path,
                    gene_file_names[[i]],
                    sep = "")
  gene_file = read.table(
    file = full_path,
    quote = "",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )
  ordered_by_length = gene_file[order(-gene_file$SV_LENGTH), ]
  # UniProt recommends max 100 keys per query or else it might fail so for files with more than 100 variants I chose the 100 longest SVs
  if (dim(ordered_by_length)[1] > 100) {
    n_lines = as.integer(dim(gene_file)[1] / 100)
    for (j in 1:n_lines) {
      if (j == 1) {
        gene_ids = ordered_by_length$gene_id[1:100]
        uniprot = UniProt.ws::select(up,
                                     keys = gene_ids,
                                     columns = all_cols,
                                     keytype = "ENSEMBL")
        uniprot$ENSEMBL = as.factor(uniprot$ENSEMBL)
        uniprot = uniprot %>% filter(REVIEWED == "reviewed")
        renamed = uniprot %>% dplyr::rename(gene_id = ENSEMBL, UNIPROT_ID = ID)
        final = renamed
      }
      else if (j == n_lines) {
        x = (j - 1) * 100 + 1
        y = dim(ordered_by_length)[1]
        gene_ids = ordered_by_length$gene_id[x:y]
        uniprot = UniProt.ws::select(up,
                                     keys = gene_ids,
                                     columns = all_cols,
                                     keytype = "ENSEMBL")
        uniprot$ENSEMBL = as.factor(uniprot$ENSEMBL)
        uniprot = uniprot %>% filter(REVIEWED == "reviewed")
        renamed = uniprot %>% dplyr::rename(gene_id = ENSEMBL, UNIPROT_ID = ID)
        final = merge(final, renamed)
      }
      else{
        x = (j - 1) * 100 + 1
        y = j * 100
        gene_ids = ordered_by_length$gene_id[x:y]
        uniprot = UniProt.ws::select(up,
                                     keys = gene_ids,
                                     columns = all_cols,
                                     keytype = "ENSEMBL")
        uniprot$ENSEMBL = as.factor(uniprot$ENSEMBL)
        uniprot = uniprot %>% filter(REVIEWED == "reviewed")
        renamed = uniprot %>% dplyr::rename(gene_id = ENSEMBL, UNIPROT_ID = ID)
        final = merge(final, renamed)
      }
      
    }
  }
  else {
    gene_ids = ordered_by_length$gene_id
    
    uniprot = UniProt.ws::select(up,
                                 keys = gene_ids,
                                 columns = all_cols,
                                 keytype = "ENSEMBL")
    uniprot$ENSEMBL = as.factor(uniprot$ENSEMBL)
    uniprot = uniprot %>% filter(REVIEWED == "reviewed")
    renamed = uniprot %>% dplyr::rename(gene_id = ENSEMBL, UNIPROT_ID = ID)
    final = renamed
  }
  joined = inner_join(gene_file, final, by = "gene_id")
  joined_list[[i]] = final
}
#############################################
#############################################
#############################################

for (i in 1:length(gene_file_names)) {
  filt_results = paste(params$workingDir,
                       params$uniprotDir,
                       "uniprot_",
                       gene_file_names[[i]],
                       sep = "")
  write.table(
    x = joined_list[[i]],
    file = filt_results,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
}
