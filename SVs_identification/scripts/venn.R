library(dplyr)
library(VennDiagram)
library(RColorBrewer)

# Read single lane files
# Read combined files
path = "~/breedmaps/SVs_identification/results/filtered_variants/"
resultsDir = "~/breedmaps/SVs_identification/results/venn/"

# The BTA individuals are named BTA125 -> BTA133
for (i in 25:33) {
  bta = list.files(path = path,
                   pattern = paste("precise_BTA1*", i, "_", sep = ""))
  lane1_path = bta[1]
  lane2_path = bta[2]
  lane3_path = bta[3]
  lane4_path = bta[4]
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
  both_joined = full_join(lane1_join, lane3_join, by = c("Location"))
  
  File1 <- both_joined[!is.na(both_joined$ID_L1), "Location"]
  File2 <- both_joined[!is.na(both_joined$ID_L2), "Location"]
  File3 <- both_joined[!is.na(both_joined$ID_L3), "Location"]
  File4 <- both_joined[!is.na(both_joined$ID_L4), "Location"]
  
  input <- list(L1 = File1,
                L2 = File2,
                L3 = File3,
                L4 = File4)
  col = brewer.pal(4, "Set2")
  venn = venn.diagram(
    x = input,
    category.names = c("L1", "L2", "L3", "L4"),
    filename = paste(resultsDir, "BTA1", i, "_venn.png", sep = ""),
    fill = col,
    main = paste("BTA1", i, sep = "")
  )
  
}
