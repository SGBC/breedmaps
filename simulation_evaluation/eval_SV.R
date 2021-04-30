# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("RSVSim")
library(RSVSim)
library("vcfR")

#########################################################################
# To load VCF files but still this package is not that useful for converting a VCF file into dataframe. 
vcf <- read.vcfR("C:/Users/Itachi Bal/Desktop/pipeline_eval/sim5_100_180_SV.vcf")
#########################################################################
head(vcf)




#-------------- Describe the paths----------------------#
truefiles <- "C:/Users/Itachi Bal/Desktop/pipeline_eval/"
vcffiles <- "C:/Users/Itachi Bal/Desktop/pipeline_eval/"


#-----------------Variants called by the tool-----------#
#querySVs <- read.table(paste(vcffiles,"sim5.bed", sep=""), header = FALSE, sep = "\t", stringsAsFactors=FALSE, fill=TRUE)
#querySVs

#-----------------True simulated variants------------#
simSVs1 <- read.csv(paste(truefiles,"deletions (1).csv" , sep=""), sep ="\t")

simSVs2 <- read.csv(paste(truefiles,"inversions (1).csv", sep=""), sep = "\t")

# The RSVsim has a method called "simulateSV", it produces true variant simualted data(CSV files)
# based on the type of SVs that we have installed into the reference genome.

#simSVs_1 <- simSVs2[,c(-7)]
#names(simSVs_1)[6] <- "BpSeq"

#names(simSVs1)
#names(simSVs2)

simSVs <- dplyr::bind_rows(simSVs1,simSVs2)
#simSVs <- rbind(simSVs1,simSVs2,simSVs3)

querySVs <- read.csv("C:/Users/Itachi Bal/Desktop/pipeline_eval/sim6.csv")


overlaps <- compareSV(querySVs , simSVs, tol=1000)

plot(overlaps)

data <- read.table(text="Name	Precise_variants	Imprecise_variants
BTA125_S1_L1_SV.vcf	68 	1474
BTA125_S1_L2_SV.vcf	74	1508
BTA125_S1_L3_SV.vcf	67	1433
BTA126_L1_SV.vcf	83	1574
BTA127_L1_SV.vcf	68	1356
BTA127_L2_SV.vcf	71	1353
BTA127_L3_SV.vcf	51	1259
BTA127_L4_SV.vcf	55	1277
BTA128_L3_SV.vcf	89	1497
BTA129_L1_SV.vcf	59	1310
BTA129_L2_SV.vcf	63	1254
BTA129_L3_SV.vcf	47	1156
BTA129_L4_SV.vcf	65	1210
BTA130_L1_SV.vcf	97	1553
BTA130_L2_SV.vcf	110	1642
BTA130_L3_SV.vcf	87	1524
BTA130_L4_SV.vcf	90	1452
BTA131_L1_SV.vcf	88	1502
BTA131_L2_SV.vcf	68	1615
BTA131_L3_SV.vcf	68	1394
BTA131_L4_SV.vcf	66	1475
BTA132_L1_SV.vcf	51	1198
BTA132_L2_SV.vcf	43	1231
BTA132_L3_SV.vcf	50	1123
BTA132_L4_SV.vcf	60	1155
", header = T)


data
hist(data$Precise_variants, xlab =  )

plot(data$Imprecise_variants, xlab= "Sample_genomes", ylab= "No. of SVs per sample", title("Imprecise variants"))


data2 <- read.table(text="BND	DEL DUP INS INV
10028	18624	3473	2	4086
", header = T)


hist(as.numeric(data2))











data %>% Precise_variants










