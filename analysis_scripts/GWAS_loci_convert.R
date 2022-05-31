#
# GWAS_loci_convert.R


library(biomaRt)
library(tidyverse)

New_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs"

###########################
### LD Block Conversion ###
###########################

ld38 <- read.table(file = paste0(CNVs_Folder, "/Berisa.EUR.hg38.bed"), header = F, sep = "\t", quote = "", stringsAsFactors = F, col.names = c("Chr", "LD_start", "LD_end"))
head(ld38)
ld38$Chr <- as.numeric(gsub(pattern = "chr", replacement = "", x = ld38$Chr))
filter(ld38, Chr == 1 & LD_start>10000000 & LD_end<20000000)
ld38 <- rbind(ld38, c(1, 12719464, 14565015)) # manually added LD block from nearest endpoints
filter(ld38, Chr == 7 & LD_start>140000000 & LD_end<148000000)
ld38 <- rbind(ld38, c(7, 141526757, 142959223)) # manually added LD block from nearest endpoints
filter(ld38, Chr == 17 & LD_start>32000000 & LD_end<45000000)
ld38 <- rbind(ld38, c(17, 36141651, 38653091)) # manually added LD block from nearest endpoints

GWAS_HRC <- read.table(file = paste0(New_Folder, "/GWAS_sig_loci_all_non_mucinous_hg38.bed"), header = F, sep = "\t", stringsAsFactors = F)
names(GWAS_HRC) <- c("Chr", "GWAS_start", "GWAS_end", "GWAS_pvalue")
head(GWAS_HRC)

GWAS_HRC$Chr <- gsub(pattern = "chr", replacement = "", x = GWAS_HRC$Chr)
head(GWAS_HRC)
GWAS_HRC <- mutate(GWAS_HRC, LD_start = NA, LD_end = NA)
head(GWAS_HRC)

for (SNP in 1:nrow(GWAS_HRC)) {
  ldsub <- filter(ld38, Chr==GWAS_HRC[SNP,1])
  head(ldsub)
  for(ldrow in 1:nrow(ldsub)) {
    if ((ldsub[ldrow,2] < GWAS_HRC[SNP,3]) && (ldsub[ldrow,3] >= GWAS_HRC[SNP,3])) {
      GWAS_HRC[SNP,5] <- ldsub[ldrow,2]
      GWAS_HRC[SNP,6] <- ldsub[ldrow,3]
    }
  }
}

summary(GWAS_HRC$LD_end)

head(GWAS_HRC)

filter(GWAS_HRC, is.na(LD_start))
nrow(unique(dplyr::select(GWAS_HRC, Chr, LD_start, LD_end)))
#filter(GWAS_HRC, is.na(LD_end))
write.table(x = GWAS_HRC, file = paste0(New_Folder, "/GWAS_Loci_HRC_coord_all_non_mucinous.txt"), row.names = F, quote = F, col.names = T, sep = "\t")

GWAS_HRC_unique <- unique(dplyr::select(GWAS_HRC, chrom = Chr, start = LD_start, end = LD_end))
write.table(x = GWAS_HRC_unique, file = paste0(New_Folder, "/GWAS_Loci_HRC_LD_unique_all_non_mucinous.txt"), row.names = F, quote = F, col.names = T, sep = "\t")



#################################
### LD Block Conversion HGSOC ###
#################################

ld38 <- read.table(file = paste0(CNVs_Folder, "/Berisa.EUR.hg38.bed"), header = F, sep = "\t", quote = "", stringsAsFactors = F, col.names = c("Chr", "LD_start", "LD_end"))
head(ld38)
ld38$Chr <- as.numeric(gsub(pattern = "chr", replacement = "", x = ld38$Chr))
filter(ld38, Chr == 1 & LD_start>10000000 & LD_end<20000000)
ld38 <- rbind(ld38, c(1, 12719464, 14565015)) # manually added LD block from nearest endpoints
filter(ld38, Chr == 7 & LD_start>140000000 & LD_end<148000000)
ld38 <- rbind(ld38, c(7, 141526757, 142959223)) # manually added LD block from nearest endpoints
filter(ld38, Chr == 17 & LD_start>32000000 & LD_end<45000000)
ld38 <- rbind(ld38, c(17, 36141651, 38653091)) # manually added LD block from nearest endpoints

GWAS_HRC_HGSOC <- read.table(file = paste0(New_Folder, "/GWAS_sig_loci_serous_hg_extra_hg38.bed"), header = F, sep = "\t", stringsAsFactors = F)
names(GWAS_HRC_HGSOC) <- c("Chr", "GWAS_start", "GWAS_end", "GWAS_pvalue")
head(GWAS_HRC_HGSOC)

GWAS_HRC_HGSOC$Chr <- gsub(pattern = "chr", replacement = "", x = GWAS_HRC_HGSOC$Chr)
head(GWAS_HRC_HGSOC)
GWAS_HRC_HGSOC <- mutate(GWAS_HRC_HGSOC, LD_start = NA, LD_end = NA)
head(GWAS_HRC_HGSOC)

for (SNP in 1:nrow(GWAS_HRC_HGSOC)) {
  ldsub <- filter(ld38, Chr==GWAS_HRC_HGSOC[SNP,1])
  head(ldsub)
  for(ldrow in 1:nrow(ldsub)) {
    if ((ldsub[ldrow,2] < GWAS_HRC_HGSOC[SNP,3]) && (ldsub[ldrow,3] >= GWAS_HRC_HGSOC[SNP,3])) {
      GWAS_HRC_HGSOC[SNP,5] <- ldsub[ldrow,2]
      GWAS_HRC_HGSOC[SNP,6] <- ldsub[ldrow,3]
    }
  }
}

summary(GWAS_HRC_HGSOC$LD_end)

head(GWAS_HRC_HGSOC)

filter(GWAS_HRC_HGSOC, is.na(LD_start))
nrow(unique(dplyr::select(GWAS_HRC_HGSOC, Chr, LD_start, LD_end)))
write.table(x = GWAS_HRC_HGSOC, file = paste0(New_Folder, "/GWAS_Loci_HRC_coord_serous_hg_extra.txt"), row.names = F, quote = F, col.names = T, sep = "\t")

GWAS_HRC_HGSOC_unique <- unique(dplyr::select(GWAS_HRC_HGSOC, chrom = Chr, start = LD_start, end = LD_end))
write.table(x = GWAS_HRC_HGSOC_unique, file = paste0(New_Folder, "/GWAS_Loci_HRC_LD_unique_serous_hg_extra.txt"), row.names = F, quote = F, col.names = T, sep = "\t")
