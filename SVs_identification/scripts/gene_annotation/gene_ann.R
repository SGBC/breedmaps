# ADDING GENE ANNOTATION FROM ENSEMBL

if(!require("tidyverse")){
  install.packages("tidyverse")
}

if (!require("BiocManager")) {
  install.packages("BiocManager")
}
if(!require("GenomicRanges")){
  BiocManager::install("GenomicRanges")
}
if(!require("rtracklayer")){
  BiocManager::install("rtracklayer")
}
if(!require("optparse")){
  install.packages("optparse")
}
library(optparse)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)


# Create inputs from the command line and save the options
options <- list( 
  make_option(c("-b", "--workingDir"), help="Base directory", default="~/breedmaps/SVs_identification/"),
  make_option(c("-a", "--annDir"), help="Annotation directory", default="data/annotation/"),
  make_option(c("-n", "--gene_ann"), help="Annotation gtf file from ENSEMBL", default="Bos_taurus.ARS-UCD1.2.103.gtf"),
  make_option(c("-s", "--scriptDir"), help="script directory", default ="scripts/gene_annotation/"),
  make_option(c("-s", "--dataDir"), help="Data directory for the filtered variants", default ="results/gene_annotation/filtered_variants/"),
  make_option(c("-r", "--resultsDir"), help="Result directory", default="results/gene_annotation/"),
  make_option(c("-f", "--functions"), help="Function file name",default="functions.R") )

params <- parse_args(OptionParser(option_list=options))
source(paste(params$workingDir, params$scriptDir, params$functions, sep=""))

########################################
########################################
########################################
# Load VCF file
vcf_path = paste(params$workingDir, params$dataDir, sep = "")
vcf_files = list.files(path = vcf_path, pattern = "precise_[A-Z]*")
vcf_names = list()
df_list = list()
range_list = list()

for (i in 1:length(vcf_files)) {
  vcf_file = paste(vcf_path, vcf_files[i], sep = "")
  
  # Save the file name to be able to use it when writing results files
  name = strsplit(vcf_files, "\\.")[[i]][1]
  name = strsplit(name, "precise_")[[1]][2]
  vcf_names[[i]] = name
  
  df = read.table(
    file = vcf_file,
    quote = "",
    sep="\t",
    header = T,
    stringsAsFactors = F
  )
  # Because of the ranges need to filter out the breakpoint
  filt_var = df %>% dplyr::filter(ALT == "<DUP>" | ALT == "<DEL>" |
                                    ALT == "<INV>" | ALT == "<INS>") %>% dplyr::rename(ID1 = ID)
  # Filtering out all variants larger than 50 000 000 since it is too computer heavy to add annotation to these variants
  filt_var_final = filt_var %>% dplyr::filter(SV_LENGTH < 50000000)
  df_list[[i]] = filt_var_final
  vcf_range = makeGRangesFromDataFrame(
    filt_var_final,
    seqnames.field = "CHROM",
    start.field = "POS",
    end.field = "END",
  )
  range_list[[i]] = vcf_range
}

########################################
########################################
########################################
# Load gene info from ENSEMBL gtf file
ann_path = paste(params$workingDir, params$annDir, params$gene_ann, sep =
                   "")
ann_df = as.data.frame(rtracklayer::import(ann_path))
ann_df = ann_df %>% dplyr::rename(gene_chrom = seqnames, gene_start = start, gene_end = end, gene_strand = strand, gene_width = width, gene_type = type) %>% dplyr::select(-source)
########################################
########################################
########################################
# The overlap function requires that the ID for the two overlapping data sets have ID1 or ID2 for the ID
gene_df = ann_df %>% filter(gene_type == "gene") %>% dplyr::rename(ID2 = gene_id)
gene_range = makeGRangesFromDataFrame(
  gene_df,
  seqnames.field = "gene_chrom",
  start.field = "gene_start",
  end.field = "gene_end",
  strand.field = "gene_strand"
)

overlap_gene_list = list()
for (i in 1:length(range_list)) {
  overlap_gene = findoverlap_dataframe(range_list[[i]], gene_range, df_list[[i]], gene_df)
  overlap_gene_list[[i]] = overlap_gene %>% dplyr::rename(ID = ID1, gene_id = ID2)
}
########################################
########################################
########################################
# The overlap function requires that the ID for the two overlapping data sets have ID1 or ID2 for the ID
tran_df = ann_df %>% filter(gene_type == "transcript") %>% dplyr::rename(ID2 = gene_id)
tran_range = makeGRangesFromDataFrame(
  tran_df,
  seqnames.field = "gene_chrom",
  start.field = "gene_start",
  end.field = "gene_end",
  strand.field = "gene_strand"
)

overlap_tran_list = list()
for (i in 1:length(range_list)) {
  overlap_tran = findoverlap_dataframe(range_list[[i]], tran_range, df_list[[i]], tran_df)
  overlap_tran_list[[i]] = overlap_tran %>% dplyr::rename(ID = ID1, gene_id = ID2)
}
########################################
########################################
########################################
exon_df = ann_df %>% filter(gene_type == "exon") %>% dplyr::rename(ID2 = gene_id)
exon_range = makeGRangesFromDataFrame(
  exon_df,
  seqnames.field = "gene_chrom",
  start.field = "gene_start",
  end.field = "gene_end",
  strand.field = "gene_strand"
)

overlap_exon_list = list()
for (i in 1:length(range_list)) {
  overlap_exon = findoverlap_dataframe(range_list[[i]], exon_range, df_list[[i]], exon_df)
  overlap_exon_list[[i]] = overlap_exon %>% dplyr::rename(ID = ID1, gene_id = ID2)
}
########################################
########################################
########################################
rest_df = ann_df %>% filter(gene_type != "exon") %>% filter(gene_type != "gene") %>% filter(gene_type != "transcript") %>% dplyr::rename(ID2 = gene_id)
rest_range = makeGRangesFromDataFrame(
  rest_df,
  seqnames.field = "gene_chrom",
  start.field = "gene_start",
  end.field = "gene_end",
  strand.field = "gene_strand"
)

overlap_rest_list = list()
for (i in 1:length(range_list)) {
  overlap_rest = findoverlap_dataframe(range_list[[i]], rest_range, df_list[[i]], rest_df)
  overlap_rest_list[[i]] = overlap_rest %>% dplyr::rename(ID = ID1, gene_id = ID2)
}
########################################
########################################
########################################
# Save only the PRECISE variants
for (i in 1:length(range_list)) {
  filt_results = paste(
    params$workingDir,
    params$resultsDir,
    "annotated_variants/",
    "filt_ann_gene_",
    vcf_names[[i]],
    ".tsv",
    sep = ""
  )
  write.table(
    x = overlap_gene_list[[i]],
    file = filt_results,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
  
  filt_results = paste(
    params$workingDir,
    params$resultsDir,
    "annotated_variants/",
    "filt_ann_exon_",
    vcf_names[[i]],
    ".tsv",
    sep = ""
  )
  write.table(
    x = overlap_exon_list[[i]],
    file = filt_results,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
  
  filt_results = paste(
    params$workingDir,
    params$resultsDir,
    "annotated_variants/",
    "filt_ann_transcript_",
    vcf_names[[i]],
    ".tsv",
    sep = ""
  )
  write.table(
    x = overlap_tran_list[[i]],
    file = filt_results,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
  filt_results = paste(
    params$workingDir,
    params$resultsDir,
    "annotated_variants/",
    "filt_ann_rest_",
    vcf_names[[i]],
    ".tsv",
    sep = ""
  )
  write.table(
    x = overlap_rest_list[[i]],
    file = filt_results,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
}


