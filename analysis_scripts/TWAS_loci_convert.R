#
# TWAS_loci_convert.R

library(biomaRt)
library(tidyverse)

New_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs"


TWAS2018 <- read.table(file = paste0(CNVs_Folder, "/TWAS_Loci_2018_simple.txt"), header = T, sep = "\t", stringsAsFactors = F)
head(TWAS2018)
TWAS2018[grep(pattern = "OBFC1", x = TWAS2018$Gene), 1] <- "STN1"
TWAS2019 <- read.table(file = paste0(CNVs_Folder, "/TWAS_Loci_2019_simple.txt"), header = T, sep = "\t", stringsAsFactors = F)
head(TWAS2019)
TWAS2019[grep(pattern = "KIAA1267", x = TWAS2019$Gene), 1] <- "KANSL1"

TWAS2018bp <- read.table(file =  paste0(CNVs_Folder, "/TWAS_Loci_2018_bp.txt"), header = T, sep = "\t", stringsAsFactors = F)
head(TWAS2018bp)
TWAS2018bp[grep(pattern = "OBFC1", x = TWAS2018bp$Gene), 1] <- "STN1"
TWAS2019bp <- read.table(file =  paste0(CNVs_Folder, "/TWAS_Loci_2019_bp.txt"), header = T, sep = "\t", stringsAsFactors = F)
head(TWAS2019bp)
TWAS2019bp[grep(pattern = "KIAA1267", x = TWAS2019bp$Gene), 1] <- "KANSL1"


TWAS2018
unique(TWAS2018$Gene) %>% length()

TWAS2019
unique(TWAS2019$Gene) %>% length()

unique(c(TWAS2018$Gene, TWAS2019$Gene))
unique(c(TWAS2018$Gene, TWAS2019$Gene)) %>% length() # 43 including RP11
unique(c(TWAS2018$Gene, TWAS2019$Gene)) %>% grep(pattern = "RP11", ignore.case = T, value = T)
unique(c(TWAS2018$Gene, TWAS2019$Gene)) %>% grep(pattern = "RP11", ignore.case = T, value = T) %>% length() # 6 RP11
# 43-6 # 37 no-RP11 genes

listEnsembl()
listDatasets(useEnsembl("ensembl"))
hmart_gene <- useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")
listDatasets(useEnsembl("snp"))
hmart_snp <- useEnsembl(biomart="snp", dataset = "hsapiens_snp")



TWAS2018gene <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'), filters = 'external_gene_name', values = TWAS2018$Gene, mart = hmart_gene)
head(TWAS2018gene)
TWAS2018gene <- filter(TWAS2018gene, chromosome_name %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"))
length(unique(TWAS2018$Gene)) ## 34
length(unique(TWAS2018gene$external_gene_name)) ## 27, now 28 with OBFC1 fixed
setdiff(unique(TWAS2018$Gene), unique(TWAS2018gene$external_gene_name))
# the missing genes are "RP11-798G7.8"  "RP11-105N13.4" "RP11-669E14.6" "RP11-138C9.1"  "RP11-6N17.6"   "RP11-403A21.1" "OBFC1"  (OBFC1 now STN1) 


TWAS2018snp <- getBM(attributes = c('chr_name', 'chrom_start', 'chrom_end', 'refsnp_id'), filters = 'snp_filter', values = TWAS2018$rsID, mart = hmart_snp)
head(TWAS2018snp)
TWAS2018snp <- filter(TWAS2018snp, chr_name %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"))
length(unique(TWAS2018$rsID)) ## 13
length(unique(TWAS2018snp$refsnp_id)) ## 13

TWAS2018gene <- rename(TWAS2018gene, Chr = chromosome_name, Gene_start = start_position, Gene_end = end_position, Gene = external_gene_name)
TWAS2018snp <- rename(TWAS2018snp, Chr = chr_name, SNP_start = chrom_start, SNP_end = chrom_end, rsID = refsnp_id)

TWAS2018coord <- merge(merge(TWAS2018, TWAS2018gene), TWAS2018snp)

write.table(x = TWAS2018coord, file = "TWAS_Loci_2018_coord.txt", row.names = F, quote = F, col.names = T, sep = "\t")



TWAS2019gene <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'), filters = 'external_gene_name', values = TWAS2019$Gene, mart = hmart_gene)
head(TWAS2019gene)
TWAS2019gene <- filter(TWAS2019gene, chromosome_name %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"))
length(unique(TWAS2019$Gene)) ## 18
length(unique(TWAS2019gene$external_gene_name)) ## 17, now 18 with KIAA1267 fixed
setdiff(unique(TWAS2019$Gene), unique(TWAS2019gene$external_gene_name))
# the missing gene is "KIAA1267"



TWAS2019snp <- getBM(attributes = c('chr_name', 'chrom_start', 'chrom_end', 'refsnp_id'), filters = 'snp_filter', values = TWAS2019$rsID, mart = hmart_snp)
head(TWAS2019snp)
TWAS2019snp <- filter(TWAS2019snp, chr_name %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"))
length(unique(TWAS2019$rsID)) ## 12
length(unique(TWAS2019snp$refsnp_id)) ## 12


head(TWAS2019)
head(TWAS2019bp)
head(TWAS2019gene)
head(TWAS2019snp)
TWAS2019gene <- rename(TWAS2019gene, Chr = chromosome_name, Gene_start = start_position, Gene_end = end_position, Gene = external_gene_name)
TWAS2019snp <- rename(TWAS2019snp, Chr = chr_name, SNP_start = chrom_start, SNP_end = chrom_end, rsID = refsnp_id)

TWAS2019coord <- merge(merge(TWAS2019, TWAS2019gene), TWAS2019snp)
head(TWAS2019coord)

write.table(x = TWAS2019coord, file = "TWAS_Loci_2019_coord.txt", row.names = F, quote = F, col.names = T, sep = "\t")



###########################
### LD Block Conversion ###
###########################

ld38 <- read.table(file = "Berisa.EUR.hg38.bed", header = F, sep = "\t", quote = "", stringsAsFactors = F, col.names = c("Chr", "LD_start", "LD_end"))
head(ld38)
ld38$Chr <- as.numeric(gsub(pattern = "chr", replacement = "", x = ld38$Chr))
filter(ld38, Chr == 1 & LD_start>10000000 & LD_end<20000000)
ld38 <- rbind(ld38, c(1, 12719464, 14565015)) # manually added LD block from nearest endpoints
filter(ld38, Chr == 7 & LD_start>140000000 & LD_end<148000000)
ld38 <- rbind(ld38, c(7, 141526757, 142959223)) # manually added LD block from nearest endpoints
filter(ld38, Chr == 17 & LD_start>32000000 & LD_end<45000000)
ld38 <- rbind(ld38, c(17, 36141651, 38653091)) # manually added LD block from nearest endpoints

GWAS_HRC <- read.table(file = "sig_all_non_mucinous_SNPs_hg38.txt", header = T, sep = "\t", stringsAsFactors = F)
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
write.table(x = GWAS_HRC, file = "GWAS_Loci_HRC_coord.txt", row.names = F, quote = F, col.names = T, sep = "\t")
nrow(GWAS_HRC)
nrow(unique(GWAS_HRC))
nrow(select(GWAS_HRC, Chr, LD_start, LD_end))
nrow(unique(select(GWAS_HRC, Chr, LD_start, LD_end))) # 35
write.table(x = unique(select(GWAS_HRC, Chr, LD_start, LD_end)), file = "GWAS_Loci_HRC_LD_unique.txt", row.names = F, quote = F, col.names = T, sep = "\t")

head(TWAS2018coord)
TWAS2018coord <- mutate(TWAS2018coord, LD_start = NA, LD_end = NA, Gene_LD_start = NA, Gene_LD_end = NA)
for (SNP in 1:nrow(TWAS2018coord)) {
  ldsub <- filter(ld38, Chr==TWAS2018coord[SNP,2])
  head(ldsub)
  for(ldrow in 1:nrow(ldsub)) {
    if ((ldsub[ldrow,2] < TWAS2018coord[SNP,8]) && (ldsub[ldrow,3] >= TWAS2018coord[SNP,8])) {
      TWAS2018coord[SNP,9] <- ldsub[ldrow,2]
      TWAS2018coord[SNP,10] <- ldsub[ldrow,3]
    }
    if ((ldsub[ldrow,2] <= TWAS2018coord[SNP,5]) && (ldsub[ldrow,3] >= TWAS2018coord[SNP,6])) {
      TWAS2018coord[SNP,11] <- ldsub[ldrow,2]
      TWAS2018coord[SNP,12] <- ldsub[ldrow,3]
    }
  }
}



filter(TWAS2018coord, is.na(LD_start))
filter(TWAS2018coord, is.na(Gene_LD_start))
write.table(x = TWAS2018coord, file = "TWAS_Loci_2018_coord.txt", row.names = F, quote = F, col.names = T, sep = "\t")
nrow(TWAS2018coord)
nrow(unique(TWAS2018coord))
nrow(select(TWAS2018coord, Chr, LD_start, LD_end)) # 28
nrow(unique(select(TWAS2018coord, Chr, LD_start, LD_end))) # 12
nrow(unique(select(TWAS2018coord, Chr, Gene_start, Gene_end))) # 28
nrow(unique(select(TWAS2018coord, Chr, Gene_LD_start, Gene_LD_end))) # 13
nrow(unique(select(TWAS2018coord, Chr, Gene_start, Gene_end, LD_start, LD_end))) # 28
write.table(x = unique(select(TWAS2018coord, Chr, Gene_start, Gene_end)), file = "TWAS_Loci_2018_Gene_unique.txt", row.names = F, quote = F, col.names = T, sep = "\t")
write.table(x = unique(select(TWAS2018coord, Chr, Gene_LD_start, Gene_LD_end)), file = "TWAS_Loci_2018_Gene_LD_unique.txt", row.names = F, quote = F, col.names = T, sep = "\t")
write.table(x = unique(select(TWAS2018coord, Chr, LD_start, LD_end)), file = "TWAS_Loci_2018_LD_unique.txt", row.names = F, quote = F, col.names = T, sep = "\t")




head(TWAS2019coord)
TWAS2019coord <- mutate(TWAS2019coord, LD_start = NA, LD_end = NA, Gene_LD_start = NA, Gene_LD_end = NA)
for (SNP in 1:nrow(TWAS2019coord)) {
  ldsub <- filter(ld38, Chr==TWAS2019coord[SNP,2])
  head(ldsub)
  for(ldrow in 1:nrow(ldsub)) {
    #print(paste("ldsub[ldrow,2] = ", ldsub[ldrow,2]))
    #print(paste("GWAS_HRC[SNP,3] = ", GWAS_HRC[SNP,3]))
    #print(paste("ldsub[ldrow,3] = ", ldsub[ldrow,3]))
    if ((ldsub[ldrow,2] < TWAS2019coord[SNP,8]) && (ldsub[ldrow,3] >= TWAS2019coord[SNP,8])) {
      TWAS2019coord[SNP,9] <- ldsub[ldrow,2]
      TWAS2019coord[SNP,10] <- ldsub[ldrow,3]
    }
    if ((ldsub[ldrow,2] <= TWAS2019coord[SNP,5]) && (ldsub[ldrow,3] >= TWAS2019coord[SNP,6])) {
      TWAS2019coord[SNP,11] <- ldsub[ldrow,2]
      TWAS2019coord[SNP,12] <- ldsub[ldrow,3]
    }
  }
}

filter(TWAS2019coord, is.na(LD_start))
filter(TWAS2019coord, is.na(Gene_LD_start))
write.table(x = TWAS2019coord, file = "TWAS_Loci_2019_coord.txt", row.names = F, quote = F, col.names = T, sep = "\t")
nrow(TWAS2019coord)
nrow(unique(TWAS2019coord))
nrow(select(TWAS2019coord, Chr, LD_start, LD_end)) # 19
nrow(unique(select(TWAS2019coord, Chr, LD_start, LD_end))) # 7
nrow(unique(select(TWAS2019coord, Chr, Gene_start, Gene_end))) # 18
nrow(unique(select(TWAS2019coord, Chr, Gene_LD_start, Gene_LD_end))) # 5
nrow(unique(select(TWAS2019coord, Chr, Gene_start, Gene_end, LD_start, LD_end))) # 18
write.table(x = unique(select(TWAS2019coord, Chr, Gene_start, Gene_end)), file = "TWAS_Loci_2019_Gene_unique.txt", row.names = F, quote = F, col.names = T, sep = "\t")
write.table(x = unique(select(TWAS2019coord, Chr, Gene_LD_start, Gene_LD_end)), file = "TWAS_Loci_2019_Gene_LD_unique.txt", row.names = F, quote = F, col.names = T, sep = "\t")
write.table(x = unique(select(TWAS2019coord, Chr, LD_start, LD_end)), file = "TWAS_Loci_2019_LD_unique.txt", row.names = F, quote = F, col.names = T, sep = "\t")

