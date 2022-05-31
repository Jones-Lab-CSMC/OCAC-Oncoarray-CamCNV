# PaperFigures.R

library(dplyr)
library(qqman)

New_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs"

tablecutoff = 5e-4

overallprobes <- read.table(file = paste0(New_Folder, "/", "ocac_overall_cnv_results_v7_all_probes_hg38.txt"), sep = "\t", header = T)
nrow(overallprobes)
overallp <- 0.05/nrow(overallprobes)
overallprobes <- mutate(overallprobes, Chr = as.numeric(gsub(pattern = "chr", replacement = "", x = gsub(pattern = "X", replacement = "23", x = Chr_numeric))))
table(overallprobes$Chr)
max(-log10(overallprobes$pvalue))

hgsocprobes <- read.table(file = paste0(New_Folder, "/", "ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt"), sep = "\t", header = T)
nrow(hgsocprobes)
hgsocp <- 0.05/nrow(hgsocprobes)
hgsocprobes <- mutate(hgsocprobes, Chr = as.numeric(gsub(pattern = "chr", replacement = "", x = gsub(pattern = "X", replacement = "23", x = Chr_numeric))))
table(hgsocprobes$Chr)
max(-log10(hgsocprobes$pvalue))


jpeg(filename = paste0(New_Folder, "/Figures/Final_Paper/", "Probe_Results_Manhattans", ".jpg"), width = 6, height = 6, units = "in", res = 1250)
par(mfrow=c(2,1))
manhattan(x = overallprobes, chr = "Chr", bp = "position_b37", p = "pvalue", snp = NA, main = "Overall Single Probes", col = c(rgb(204, 0, 0, maxColorValue = 255), 
                                                    rgb(255, 153, 153, maxColorValue = 255)), 
          suggestiveline = -log10(overallp), genomewideline = -log10(tablecutoff), ylim = c(0, 12), cex = 1.1)
manhattan(x = hgsocprobes, chr = "Chr", bp = "position_b37", p = "pvalue", snp = NA, main = "HGSOC Single Probes", col = c(rgb(102, 0, 204, maxColorValue = 255), 
                                                     rgb(204, 153, 255, maxColorValue = 255)), 
          suggestiveline = -log10(hgsocp), genomewideline = -log10(tablecutoff), ylim = c(0, 12), cex = 1.1)
dev.off()





all_non_mucinous <- read.table(paste0(CNVs_Folder, "/lift_all_non_mucinous_HRC_stats_hg38.txt"), 
                               na.strings = c("", "NA", "-99"), header = T, fill = T, stringsAsFactors = F)
data1 <- all_non_mucinous %>% dplyr::select(snp, chrchromosome, position, pvalue) %>% 
  dplyr::rename(SNP = snp, CHR = chrchromosome, BP = position, P = pvalue) %>% na.omit()
data1$CHR <- as.numeric(gsub("chr", "", data1$CHR))

endometrioid <- read.table(paste0(CNVs_Folder, "/lift_endometrioid_HRC_stats_hg38.txt"), 
                           na.strings = c("", "NA", "-99"), header = T, fill = T, stringsAsFactors = F)
data2 <- endometrioid %>% dplyr::select(snp, chrchromosome, position, pvalue) %>% 
  dplyr::rename(SNP = snp, CHR = chrchromosome, BP = position, P = pvalue) %>% na.omit()
data2$CHR <- as.numeric(gsub("chr", "", data2$CHR))

serous_hg_extra <- read.table(paste0(CNVs_Folder, "/lift_serous_hg_extra_HRC_stats_hg38.txt"), 
                              na.strings = c("", "NA", "-99"), header = T, fill = T, stringsAsFactors = F)
data3 <- serous_hg_extra %>% dplyr::select(snp, chrchromosome, position, pvalue) %>% 
  dplyr::rename(SNP = snp, CHR = chrchromosome, BP = position, P = pvalue) %>% na.omit()
data3$CHR <- as.numeric(gsub("chr", "", data3$CHR))

length(unique(all_non_mucinous$snp))
length(unique(endometrioid$snp))
length(unique(serous_hg_extra$snp))
length(unique(c(all_non_mucinous$snp, endometrioid$snp, serous_hg_extra$snp)))
desiredp = 0.05/23960
jpeg(filename = paste0("Figures/Final_Paper/", "HRC_tagSNPs_hg38_stacked_manhattans", ".jpg"), width = 6,
     height = 6, units = "in", res = 1250)
par(mfrow=c(3,1))
par(mfrow=c(2,1))
manhattan(data1, main = "All Non Mucinous", col = c(rgb(204, 0, 0, maxColorValue = 255), 
                                                    rgb(255, 153, 153, maxColorValue = 255)), 
          suggestiveline = -log10(desiredp), genomewideline = F, ylim = c(0, 70), cex = 1.1)
manhattan(data3, main = "High Grade Serous", col = c(rgb(102, 0, 204, maxColorValue = 255), 
                                                     rgb(204, 153, 255, maxColorValue = 255)), 
          suggestiveline = -log10(desiredp), genomewideline = F, ylim = c(0, 70), cex = 1.1)
dev.off()

filter(data1, CHR==10 & P<1e-4)
filter(data1, P<desiredp)
filter(data2, P<desiredp)
filter(data3, P<desiredp)

subdata1 <- all_non_mucinous %>% dplyr::select(snp, chrchromosome, position, pvalue, min_individual_p) %>% 
  dplyr::rename(SNP = snp, CHR = chrchromosome, BP = position, P = pvalue, MINP = min_individual_p) %>% na.omit()
subdata1$CHR <- as.numeric(gsub("chr", "", subdata1$CHR))
head(subdata1)

nrow(filter(subdata1, P<desiredp))
nrow(filter(subdata1, MINP<desiredp))
nrow(filter(subdata1, P<desiredp & MINP<desiredp))
filter(subdata1, P>desiredp & MINP<desiredp)
filter(subdata1, P<desiredp & MINP<desiredp)
min(filter(subdata1, CHR==17)$P)
filter(subdata1, P == min(filter(subdata1, CHR==17)$P))
min(filter(subdata1, CHR==9)$P)
filter(subdata1, P == min(filter(subdata1, CHR==9)$P))

subdata2 <- endometrioid %>% dplyr::select(snp, chrchromosome, position, pvalue, min_individual_p) %>% 
  dplyr::rename(SNP = snp, CHR = chrchromosome, BP = position, P = pvalue, MINP = min_individual_p) %>% na.omit()
subdata2$CHR <- as.numeric(gsub("chr", "", subdata2$CHR))
head(subdata2)

nrow(filter(subdata2, P<desiredp))
nrow(filter(subdata2, MINP<desiredp))
nrow(filter(subdata2, P<desiredp & MINP<desiredp))
filter(subdata2, P>desiredp & MINP<desiredp)
min(filter(subdata2, CHR==17)$P)
filter(subdata2, P == min(filter(subdata2, CHR==17)$P))
min(filter(subdata2, CHR==9)$P)
filter(subdata2, P == min(filter(subdata2, CHR==9)$P))



subdata3 <- serous_hg_extra %>% dplyr::select(snp, chrchromosome, position, pvalue, min_individual_p) %>% 
  dplyr::rename(SNP = snp, CHR = chrchromosome, BP = position, P = pvalue, MINP = min_individual_p) %>% na.omit()
subdata3$CHR <- as.numeric(gsub("chr", "", subdata3$CHR))
head(subdata3)

nrow(filter(subdata3, P<desiredp))
nrow(filter(subdata3, MINP<desiredp))
nrow(filter(subdata3, P<desiredp & MINP<desiredp))
filter(subdata3, P>desiredp & MINP<desiredp)
min(filter(subdata3, CHR==17)$P)
filter(subdata3, P == min(filter(subdata3, CHR==17)$P))
min(filter(subdata3, CHR==9)$P)
filter(subdata3, P == min(filter(subdata3, CHR==9)$P))



all_non_mucinous_CNVs <- read.table(paste0(CNVs_Folder, "/lift_all_non_mucinous_CNVs_SNPs_hg38.txt"), 
                                    na.strings = c("", "NA", "-99"), header = T, fill = T, stringsAsFactors = F)
serous_hg_extra_CNVs <- read.table(paste0(CNVs_Folder, "/lift_serous_hg_extra_CNVs_SNPs_hg38.txt"), 
                                    na.strings = c("", "NA", "-99"), header = T, fill = T, stringsAsFactors = F)

sig_all_non_mucinous_CNVs <- filter(all_non_mucinous_CNVs, pvalue < desiredp)
sig_all_non_mucinous_CNVs
unique(sig_all_non_mucinous_CNVs$CNP)
