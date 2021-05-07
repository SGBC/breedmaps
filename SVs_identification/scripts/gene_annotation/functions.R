####################--------------- FUNCTIONS FOR GENE ANNOTATION OF CATTLE VCFs ---------------####################

# Function to read DELLY vcf files and extract the PRECISE flag. (vcfR was unable to handle the PRECISE/IMPRECISE flag)
read_BTA_delly_vcf = function(path) {
  vcf = read.table(
    file = path,
    quote = "",
    sep = "\t",
    header = FALSE,
    stringsAsFactors = F
  )
  # fix = V8, Meta = V9, gt = V10
  renamed = vcf %>% dplyr::rename(
    CHROM = V1,
    POS = V2,
    ID = V3,
    REF = V4,
    ALT = V5,
    QUAL = V6,
    FILTER = V7,
    INFO = V8,
    FORMAT = V9,
    GT = V10
  )
  tmp.split = data.frame(stringr::str_split_fixed(renamed$INFO, ";", n = Inf))
  info = data.frame(matrix(
    data = NA,
    nrow = dim(tmp.split)[1],
    ncol = dim(tmp.split)[2] + 1
  ))
  rownames(info) = renamed$ID
  colnames(info) = c(
    "PRECISE",
    "IMPRECISE",
    "SVTYPE",
    "SVMETHOD",
    "CHR2",
    "END",
    "PE",
    "MAPQ",
    "CT",
    "CIPOS",
    "CIEND",
    "SRMAPQ",
    "INSLEN",
    "HOMLEN",
    "SR",
    "SRQ",
    "CONSENSUS",
    "CE"
  )
  for (row in 1:dim(info)[1]) {
    for (col in 1:dim(info)[2]) {
      if (is.null(tmp.split[row, col])) {
        
      }
      else if (tmp.split[row, col] == "") {
        
      }
      else if (tmp.split[row, col] == "PRECISE") {
        info[row, "PRECISE"] = TRUE
        info[row, "IMPRECISE"] = FALSE
      }
      else if (tmp.split[row, col] == "IMPRECISE") {
        info[row, "PRECISE"] = FALSE
        info[row, "IMPRECISE"] = TRUE
      }
      else {
        split = stringr::str_split_fixed(tmp.split[row, col], "=", n = 2)
        #info[row, split[1]] = tmp.split[row, col]
        info[row, split[1]] = split[[2]]
      }
    }
  }
  id_column = rownames_to_column(info, var = "ID")
  result = inner_join(renamed, id_column, by = "ID")
  cleaned = result %>% dplyr::select(
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILTER,
    PRECISE,
    SVTYPE,
    CHR2,
    END,
    PE,
    MAPQ,
    CT,
    CIPOS,
    CIEND,
    SRMAPQ,
    INSLEN,
    HOMLEN,
    SR,
    SRQ,
    CONSENSUS,
    CE
  )
  return (cleaned)
}


read_RDC_delly_vcf = function(path) {
  vcf = read.table(
    file = path,
    quote = "",
    sep = "\t",
    header = FALSE,
    stringsAsFactors = F
  )
  # fix = V8, Meta = V9, gt = V10
  renamed = vcf %>% dplyr::rename(
    CHROM = V1,
    POS = V2,
    ID = V3,
    REF = V4,
    ALT = V5,
    QUAL = V6,
    FILTER = V7,
    INFO = V8,
    FORMAT = V9,
    GT = V10
  )
  tmp.split = data.frame(stringr::str_split_fixed(renamed$INFO, ";", n = Inf))
  info = data.frame(matrix(
    data = NA,
    nrow = dim(tmp.split)[1],
    ncol = dim(tmp.split)[2] + 1
  ))
  rownames(info) = renamed$ID
  colnames(info) = c(
    "PRECISE",
    "IMPRECISE",
    "SVTYPE",
    "SVMETHOD",
    "CHR2",
    "END",
    "PE",
    "MAPQ",
    "CT",
    "CIPOS",
    "CIEND",
    "SRMAPQ",
    "INSLEN",
    "HOMLEN",
    "SR",
    "SRQ",
    "CONSENSUS",
    "CE"
  )
  for (row in 1:dim(info)[1]) {
    for (col in 1:dim(info)[2]) {
      if (is.null(tmp.split[row, col])) {
        
      }
      else if (tmp.split[row, col] == "") {
        
      }
      else if (tmp.split[row, col] == "PRECISE") {
        info[row, "PRECISE"] = TRUE
        info[row, "IMPRECISE"] = FALSE
      }
      else if (tmp.split[row, col] == "IMPRECISE") {
        info[row, "PRECISE"] = FALSE
        info[row, "IMPRECISE"] = TRUE
      }
      else {
        split = stringr::str_split_fixed(tmp.split[row, col], "=", n = 2)
        #info[row, split[1]] = tmp.split[row, col]
        info[row, split[1]] = split[[2]]
      }
    }
  }
  info = info %>% dplyr::select(
    PRECISE,
    SVTYPE,
    CHR2,
    END,
    PE,
    MAPQ,
    CT,
    CIPOS,
    CIEND,
    SRMAPQ,
    INSLEN,
    HOMLEN,
    SR,
    SRQ,
    CONSENSUS,
    CE,
    SVLEN,
    POS2
  )
  id_column = rownames_to_column(info, var = "ID")
  result = inner_join(renamed, id_column, by = "ID")
  cleaned = result %>% dplyr::select(
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILTER,
    PRECISE,
    SVTYPE,
    CHR2,
    END,
    PE,
    MAPQ,
    CT,
    CIPOS,
    CIEND,
    SRMAPQ,
    INSLEN,
    HOMLEN,
    SR,
    SRQ,
    CONSENSUS,
    CE,
    SVLEN,
    POS2
  )
  return (cleaned)
}

findoverlap_dataframe = function(var_range, gene_range, var_df, gene_df){
  overlap = suppressWarnings(data.frame(findOverlaps(var_range, gene_range, maxgap = 0)))
  # The output from findOverlaps needs to be reformatted to include the geneIDs and variantIDs we need to join
  overlap$ID1 = var_df$ID1[overlap$queryHits]
  overlap$ID2 = gene_df$ID2[overlap$subjectHits]
  overlap = overlap %>% dplyr::select(ID1, ID2)
  joined_vcf  = inner_join(overlap, var_df, by = "ID1")
  joined_genes =  inner_join(joined_vcf,  gene_df, by = "ID2")
}

