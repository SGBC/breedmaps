##################  Summarizing overlaps for each SVs  ############################
# This script summarize the number of overlaps for each external dataset, regulatory regions and Ensembl annotation
# It groups by SV_ID which is unique for each file/individual


if (!require("tidyverse")) {
  install.packages("tidyverse")
}
if (!require("optparse")) {
  install.packages("optparse")
}
library(optparse)
library(tidyverse)


options <- list(
  make_option(c("-b", "--workingDir"), help = "Base directory", default =
                "~/breedmaps/SVs_identification/"),
  make_option(c("-r", "--resDir"), help = "Results directory", default =
                "results/"),
  make_option(c("-f", "--filtDir"), help = "Filtered variants directory", default =
                "filtered_variants/"),
  make_option(c("-a", "--annDir"), help = "Annotated variants directory", default =
                "annotated_variants/"),
  make_option(c("-d", "--datasetDir"), help = "Datasets results directory", default =
                "datasets/"),
  make_option(c("-g", "--reglDir"), help = "Regulatory regions results directory", default =
                "regulatory_variants/"),
  make_option(c("-n", "--fileName"), help = "File name/individual", default =
                "combined_BTA125_S1_R_SV.vcf")
)


params <- parse_args(OptionParser(option_list = options))

##################  Load the annotated SVs  ############################
# Extract name for result file
name = strsplit(params$fileName, "\\.")[[1]][1]

sv_path = paste(
  params$workingDir,
  params$resDir,
  params$filtDir,
  "precise_",
  name,
  ".tsv",
  sep = ""
)
sv_df = read.table(
  file = sv_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(ID) %>% dplyr::rename(SV_ID = ID)


# Load the Ensembl annotated SVs, output from ensembl_ann.R
gene_path = paste(
  params$workingDir,
  params$resDir,
  params$annDir,
  "filt_ann_gene_",
  name,
  ".tsv",
  sep = ""
)
gene_df = read.table(
  file = gene_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(ID, gene_id) %>% dplyr::rename(SV_ID = ID, ID_gene = gene_id)

exon_path = paste(
  params$workingDir,
  params$resDir,
  params$annDir,
  "filt_ann_exon_",
  name,
  ".tsv",
  sep = ""
)

exon_df = read.table(
  file = exon_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(ID, gene_id) %>% dplyr::rename(SV_ID = ID, ID_exon = gene_id)

tran_path = paste(
  params$workingDir,
  params$resDir,
  params$annDir,
  "filt_ann_transcript_",
  name,
  ".tsv",
  sep = ""
)

tran_df = read.table(
  file = tran_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(ID, gene_id) %>% dplyr::rename(SV_ID = ID, ID_tran = gene_id)

rest_path = paste(
  params$workingDir,
  params$resDir,
  params$annDir,
  "filt_ann_rest_",
  name,
  ".tsv",
  sep = ""
)

rest_df = read.table(
  file = rest_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(ID, gene_id) %>% dplyr::rename(SV_ID = ID, ID_rest = gene_id)

# Load the overlapping regulatory regions, output from regulatory_regions.R
regl_path = paste(
  params$workingDir,
  params$resDir,
  params$reglDir,
  "regl_precise_",
  name,
  ".tsv",
  sep = ""
)

regl_df = read.table(
  file = regl_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) 
# The regulatory regions does not have ID so one is created to be able to count
regl_df$ID_regl = rownames(regl_df)
regl_df = regl_df %>% dplyr::select(SV_ID, ID_regl)

##################  Load the SVs overlapping with datasets  ############################
# Load the result file for each dataset(known SVs from Ensembl and studies in EVA), output from the sv_datasets.R
dataset_path = paste(params$workingDir,
                     params$resDir,
                     params$datasetDir,
                     sep = "")
d1_path = paste(
  dataset_path,
  "bos_taurus_structural_variations.gvf_precise_",
  name,
  ".tsv",
  sep = ""
)
d1_df = read.table(
  file = d1_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(SV_ID, ID) %>% dplyr::rename(ID_df1 = ID)

d2_path = paste(
  dataset_path,
  "remapped_estd223_Boussaha_et_al_2015.2015-11-02.Bos_taurus_UMD_3.1.Submitted.gff_precise_",
  name,
  ".tsv",
  sep = ""
)
d2_df = read.table(
  file = d2_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(SV_ID, ID) %>% dplyr::rename(ID_df2 = ID)

d3_path = paste(
  dataset_path,
  "remapped_estd234_Mesbah-Uddin_et_al_2017.2018-06-01.Bos_taurus_UMD_3.1.Submitted.gff_precise_",
  name,
  ".tsv",
  sep = ""
)
d3_df = read.table(
  file = d3_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(SV_ID, ID) %>% dplyr::rename(ID_df3 = ID)

d4_path = paste(
  dataset_path,
  "remapped_nstd56_Liu_et_al_2010.2016-07-08.Bos_taurus_UMD_3.1.1.Remapped.gff_precise_",
  name,
  ".tsv",
  sep = ""
)
d4_df = read.table(
  file = d4_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(SV_ID, ID) %>% dplyr::rename(ID_df4 = ID)

d5_path = paste(
  dataset_path,
  "remapped_nstd61_Hou_et_al_2011b.2016-07-11.Bos_taurus_UMD_3.1.1.Remapped.gff_precise_",
  name,
  ".tsv",
  sep = ""
)
d5_df = read.table(
  file = d5_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(SV_ID, ID) %>% dplyr::rename(ID_df5 = ID)

d6_path = paste(
  dataset_path,
  "remapped_nstd69_Bickhart_et_al_2012.2016-07-05.Bos_taurus_UMD_3.1.1.Remapped.gff_precise_",
  name,
  ".tsv",
  sep = ""
)
d6_df = read.table(
  file = d6_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(SV_ID, ID) %>% dplyr::rename(ID_df6 = ID)

d7_path = paste(
  dataset_path,
  "remapped_nstd135_Karimi_et_al_2016.2017-02-08.Bos_taurus_UMD_3.1.1.Remapped.gff_precise_",
  name,
  ".tsv",
  sep = ""
)
d7_df = read.table(
  file = d7_path,
  quote = "",
  sep = "\t",
  header = T,
  stringsAsFactors = F
) %>% dplyr::select(SV_ID, ID) %>% dplyr::rename(ID_df7 = ID)


##################  Count each file for the SVs  ############################

gene_count = as.data.frame(gene_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(Genes = n)
exon_count = as.data.frame(exon_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(Exons = n)
tran_count = as.data.frame(tran_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(Transcipts = n)
rest_count = as.data.frame(rest_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(Rest_of_ann = n)
regl_count = as.data.frame(regl_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(Regulatory_regions = n)

d1_count = as.data.frame(d1_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(Ensembl_SVs = n)
d2_count = as.data.frame(d2_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(estd223 = n)
d3_count = as.data.frame(d3_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(estd234 = n)
d4_count = as.data.frame(d4_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(nstd56 = n)
d5_count = as.data.frame(d5_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(nstd61 = n)
d6_count = as.data.frame(d6_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(nstd69 = n)
d7_count = as.data.frame(d7_df %>% group_by(SV_ID) %>% count()) %>% dplyr::rename(nstd135 = n)

##################  Join all of the counts for each SV  ############################
join_gene = full_join(sv_df ,gene_count, by="SV_ID")
join_exon = full_join(join_gene, exon_count, by="SV_ID")
join_tran = full_join(join_gene, tran_count, by="SV_ID")
join_rest = full_join(join_tran, rest_count, by="SV_ID")
join_regl = full_join(join_rest, regl_count, by="SV_ID")
join1 = full_join(join_regl, d1_count, by = "SV_ID")
join2 = full_join(join1, d2_count, by = "SV_ID")
join3 = full_join(join2, d3_count, by = "SV_ID")
join4 = full_join(join3, d4_count, by = "SV_ID")
join5 = full_join(join4, d5_count, by = "SV_ID")
join6 = full_join(join5, d6_count, by = "SV_ID")
join7 = full_join(join6, d7_count, by = "SV_ID")
join7[is.na(join7)] = 0

##################  Write the counts to file  ############################

result_path = paste(params$workingDir,
                    params$resDir,
                    "counts_per_file/",
                    "SV_overlap_count_",
                    name,
                    ".tsv",
                    sep = "")
write.table(
  x = join7,
  file = result_path,
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)

