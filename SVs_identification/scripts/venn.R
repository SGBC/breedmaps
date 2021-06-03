library(tidyverse)
library(VennDiagram)
library(RColorBrewer)

# Read single lane files
# Read combined files
path = "~/breedmaps/SVs_identification/results/filtered_variants/"
resultsDir = "~/breedmaps/SVs_identification/results/venn/"

source("~/breedmaps/SVs_identification/scripts/functions.R")

# The BTA individuals are named BTA125 -> BTA133
for (i in 25:33) {
  bta = list.files(path = path,
                   pattern = paste("precise_BTA1*", i, "_", sep = ""))
  bta_comb_path = list.files(path = path,
                   pattern = paste("precise_combined_BTA1*", i, "_", sep = ""))
  lane1_path = bta[1]
  lane2_path = bta[2]
  lane3_path = bta[3]
  lane4_path = bta[4]
  bta_comb = read.table(
    file = paste(path, bta_comb_path, sep = ""),
    quote = "",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )
  bta_comb_cleaned = bta_comb %>% dplyr::mutate(Location = paste(
    as.factor(CHROM),
    ":",
    as.factor(POS),
    "-",
    as.factor(END) ,
    sep = ""
  )) %>% dplyr::rename(ID_comb = ID) %>% dplyr::select(Location, ID_comb)
  
  lane1_df = read.table(
    file = paste(path, lane1_path, sep = ""),
    quote = "",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )
  lane1_cleaned = lane1_df %>% dplyr::mutate(Location = paste(
    as.factor(CHROM),
    ":",
    as.factor(POS),
    "-",
    as.factor(END) ,
    sep = ""
  )) %>% dplyr::select(Location, ID)
  
  lane2_df = read.table(
    file = paste(path, lane2_path, sep = ""),
    quote = "",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )
  lane2_cleaned = lane2_df %>% dplyr::mutate(Location = paste(
    as.factor(CHROM),
    ":",
    as.factor(POS),
    "-",
    as.factor(END) ,
    sep = ""
  )) %>% dplyr::select(Location, ID)
  
  lane3_df = read.table(
    file = paste(path, lane3_path, sep = ""),
    quote = "",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )
  lane3_cleaned = lane3_df %>% dplyr::mutate(Location = paste(
    as.factor(CHROM),
    ":",
    as.factor(POS),
    "-",
    as.factor(END) ,
    sep = ""
  )) %>% dplyr::select(Location, ID)
  
  lane4_df = read.table(
    file = paste(path, lane4_path, sep = ""),
    quote = "",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )
  lane4_cleaned = lane4_df %>% dplyr::mutate(Location = paste(
    as.factor(CHROM),
    ":",
    as.factor(POS),
    "-",
    as.factor(END) ,
    sep = ""
  )) %>% dplyr::select(Location, ID)
  lane1_join = full_join(lane1_cleaned,
                         lane2_cleaned,
                         by = "Location",
                         suffix = c("_L1", "_L2"))
  lane3_join = full_join(lane3_cleaned,
                         lane4_cleaned,
                         by = "Location",
                         suffix = c("_L3", "_L4"))
  both_joined = full_join(lane1_join, lane3_join, by = "Location")
  all_joined = full_join(both_joined, bta_comb_cleaned, by = "Location")
  
  comb = all_joined[!is.na(all_joined$ID_comb),"Location"]
  File1 = all_joined[!is.na(all_joined$ID_L1), "Location"]
  File2 = all_joined[!is.na(all_joined$ID_L2), "Location"]
  File3 = all_joined[!is.na(all_joined$ID_L3), "Location"]
  File4 = all_joined[!is.na(all_joined$ID_L4), "Location"]
  
  input = list(L1 = File1,
                L2 = File2,
                L3 = File3,
                L4 = File4,
               comb = comb)
  col = brewer.pal(5, "Set2")
  venn = venn.diagram(
    x = input,
    category.names = c("L1", "L2", "L3", "L4", "comb"),
    filename = paste(resultsDir, "BTA1", i, "_venn.jpg", sep = ""),
    fill = col,
    main = paste("BTA1", i, sep = "")
  )
  
}
