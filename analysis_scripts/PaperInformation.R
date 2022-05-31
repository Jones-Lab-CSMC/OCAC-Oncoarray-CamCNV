# PaperInformation.R

library(tidyverse)

New_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs"

phenotypes <- read.table(file = paste0(New_Folder, "/oncid_phenotypes_selected_samples_20200910.txt"), header = T, sep = "\t")
overallcnvs <- read.table(file = paste0(New_Folder, "/ocac_overall_cnvcalls_hg38.txt"), header = T)
hgsoccnvs <- read.table(file = paste0(New_Folder, "/ocac_hgsoc_cnvcalls_hg38.txt"), header = T)

overallprobes <- read.table(file = paste0(New_Folder, "/ocac_overall_cnv_results_v7_all_probes_hg38.txt"), header = T, sep = "\t")
hgsocprobes <- read.table(file = paste0(New_Folder, "/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt"), header = T, sep = "\t")
nrow(unique(select(overallprobes, Chr_numeric, position_b37)))
nrow(unique(select(overallprobes, Chr_numeric, position_b37, del_dup))) # 57413 tests
0.05/nrow(unique(select(overallprobes, Chr_numeric, position_b37, del_dup))) #8.708829e-07
nrow(unique(select(hgsocprobes, Chr_numeric, position_b37, del_dup)))
0.05/nrow(unique(select(hgsocprobes, Chr_numeric, position_b37, del_dup))) #8.57942e-07
summary(overallprobes$cnv_count)
summary(hgsocprobes$cnv_count)
nrow(unique(select(filter(overallprobes, del_dup=="del"), Chr_numeric, position_b37)))
nrow(unique(select(filter(overallprobes, del_dup=="dup"), Chr_numeric, position_b37)))
nrow(unique(select(filter(hgsocprobes, del_dup=="del"), Chr_numeric, position_b37)))
nrow(unique(select(filter(hgsocprobes, del_dup=="dup"), Chr_numeric, position_b37)))


nrow(unique(select(overallprobes, Chr_numeric, position_b37, del_dup))) - nrow(unique(select(overallprobes, Chr_numeric, position_b37)))
# 5515 both
nrow(unique(select(filter(overallprobes, del_dup=="del"), Chr_numeric, position_b37))) - (nrow(unique(select(overallprobes, Chr_numeric, position_b37, del_dup))) - nrow(unique(select(overallprobes, Chr_numeric, position_b37))))
# 20981 dels only
nrow(unique(select(filter(overallprobes, del_dup=="dup"), Chr_numeric, position_b37))) - (nrow(unique(select(overallprobes, Chr_numeric, position_b37, del_dup))) - nrow(unique(select(overallprobes, Chr_numeric, position_b37))))
# 25402 dups only

overallsig01cnvs <- read.table(file = paste0(New_Folder, "/ocac_overall_sig01_CNVs_hg38.txt"), header = T, sep="\t")
hgsocsig01cnvs <- read.table(file = paste0(New_Folder, "/ocac_hgsoc_sig01_CNVs_hg38.txt"), header = T, sep="\t")

overallcnvsids <- read.table(file = paste0(New_Folder, "/ocac_overall_cnvcalls_IDs_list_hg38.txt"), header = T)

brca1dels <- read.table(file = paste0(New_Folder, "/BRCA1_all_dels.txt"), header = T)
brca1dups <- read.table(file = paste0(New_Folder, "/BRCA1_all_dups.txt"), header = T)

tagSNPs <- read.table(file = paste0(CNVs_Folder, "/tagSNPs_CEULD_hg38.txt"), header = T)
all_non_mucinous <- read.table(paste0(CNVs_Folder,"/lift_all_non_mucinous_HRC_stats_hg38.txt"), 
                               na.strings = c("", "NA", "-99"), header = T, fill = T, stringsAsFactors = F)
endometrioid <- read.table(paste0(CNVs_Folder,"/lift_endometrioid_HRC_stats_hg38.txt"), 
                           na.strings = c("", "NA", "-99"), header = T, fill = T, stringsAsFactors = F)
serous_hg_extra <- read.table(paste0(CNVs_Folder,"/lift_serous_hg_extra_HRC_stats_hg38.txt"), 
                              na.strings = c("", "NA", "-99"), header = T, fill = T, stringsAsFactors = F)


overallcnvs <- merge(overallcnvs, select(phenotypes, OncID=onc_id, overall, serous_hg_extra), by = "OncID")
head(overallcnvs)

nrow(unique(select(overallcnvs, Chr, startpos, endpos))) # 60449
nrow(unique(select(overallcnvs, Chr, startpos, endpos, del_dup))) # 61459
nrow(unique(select(filter(overallcnvs, overall==0), Chr, startpos, endpos, del_dup)))
nrow(unique(select(filter(overallcnvs, overall==1), Chr, startpos, endpos, del_dup)))


nrow(overallprobes)
nrow(unique(select(overallprobes, Chr_numeric, position_b37)))
nrow(overallprobes) - nrow(unique(select(overallprobes, Chr_numeric, position_b37)))
nrow(filter(overallprobes, del_dup=="del"))
nrow(filter(overallprobes, del_dup=="dup"))

AllCNVsPerOncID <- data.frame(OncID = unique(phenotypes$onc_id), count=NA)
CNVsPerOncID <- data.frame(table(overallcnvs$OncID))
names(CNVsPerOncID) <- c("OncID", "count")
for (id in AllCNVsPerOncID$OncID) {
  if(id %in% CNVsPerOncID$OncID) {
    AllCNVsPerOncID[which(AllCNVsPerOncID$OncID==id), 'count'] <- CNVsPerOncID[which(CNVsPerOncID$OncID==id), 'count']
  }
}
summary(AllCNVsPerOncID$count)

overallcnvs <- unite(overallcnvs, col = CNVID, Chr, startpos, endpos, del_dup, remove = F)
OncIDsPerCNV <- data.frame(table(overallcnvs$CNVID))
names(OncIDsPerCNV) <- c("CNVID", "count")
summary(OncIDsPerCNV$count)
filter(OncIDsPerCNV, count >500)
nrow(filter(OncIDsPerCNV, count>1))

summary(brca1dels$length)
summary(brca1dups$length)


# number tagSNPs and tagged CNVs
length(unique(tagSNPs$CNP))
length(unique(tagSNPs$SNP))
tagSNPs <- unite(data = tagSNPs, col = "CHR_SNPPOS", CHR, SNPPOS, sep = ":", remove = F)
all_non_mucinous <- mutate(all_non_mucinous, CHR=gsub(pattern = "chr", replacement = "", x = chrchromosome)) %>% unite(., col = "CHR_SNPPOS", CHR, position, sep = ":", remove = F)
endometrioid <- mutate(endometrioid, CHR=gsub(pattern = "chr", replacement = "", x = chrchromosome)) %>% unite(., col = "CHR_SNPPOS", CHR, position, sep = ":", remove = F)
serous_hg_extra <- mutate(serous_hg_extra, CHR=gsub(pattern = "chr", replacement = "", x = chrchromosome)) %>% unite(., col = "CHR_SNPPOS", CHR, position, sep = ":", remove = F)
tagSNPslist <- unique(c(all_non_mucinous$CHR_SNPPOS, endometrioid$CHR_SNPPOS, serous_hg_extra$CHR_SNPPOS))


length(unique(c(all_non_mucinous$snp, endometrioid$snp, serous_hg_extra$snp))) # 23996
length(tagSNPslist) # 23960
nrow(unique(select(filter(tagSNPs, CHR_SNPPOS %in% tagSNPslist), CHR, START, END)))
length(unique(tagSNPs$CHR_SNPPOS))
nrow(unique(select(filter(tagSNPs, CHR_SNPPOS %in% tagSNPslist), CHR_SNPPOS)))

# number of unique bonferroni significant probes and cnvs
nrow(unique(select(filter(overallprobes, pvalue<8.71e-7), Chr_numeric, position_b37, del_dup)))
nrow(unique(select(filter(hgsocprobes, pvalue<8.56e-7), Chr_numeric, position_b37, del_dup)))

nrow(unique(select(filter(overallprobes, pvalue<8.71e-7), Chr_numeric, position_b37, del_dup)))
unique(filter(overallprobes, pvalue<8.71e-7)$position_b37)
unique(filter(hgsocprobes, pvalue<8.56e-7)$position_b37)
length(intersect(unique(filter(overallprobes, pvalue<8.71e-7)$position_b37),
          unique(filter(hgsocprobes, pvalue<8.56e-7)$position_b37)))
length(union(unique(filter(overallprobes, pvalue<8.71e-7)$position_b37),
                 unique(filter(hgsocprobes, pvalue<8.56e-7)$position_b37)))

nrow(unique(select(filter(overallsig01cnvs, pvalue<8.71e-7), Chr, startpos, endpos, del_dup)))
nrow(filter(overallsig01cnvs, pvalue<8.71e-7))

nrow(unique(select(filter(hgsocsig01cnvs, pvalue<8.56e-7), Chr, startpos, endpos, del_dup)))
nrow(filter(hgsocsig01cnvs, pvalue<8.56e-7))
nrow(filter(filter(hgsocsig01cnvs, pvalue<8.56e-7), del_dup==del_dup))

controlBRCA1 <- c(rep.int(x = 1, times = 17), rep.int(x = 0, times = (17306-17))) # 17/17306
caseBRCA1 <- c(rep.int(x = 1, times = 105), rep.int(x = 0, times = (13071-105))) # 105/13071
HGSOCBRCA1 <- c(rep.int(x = 1, times = 93), rep.int(x = 0, times = (8679-93))) # 93/8679
t.test(x = controlBRCA1, y = caseBRCA1)
t.test(x = controlBRCA1, y = HGSOCBRCA1)

OverallGroup <- data.frame(BRCA = c(rep(x = "CNV", times = 17), rep(x = "None", times = (17306-17)), rep(x = "CNV", times = 105), rep(x = "None", times = (13071-105))), CaseCtrl = c(rep(x = "Ctrl", times = 17306), rep(x = "Case", times = 13071)))
HGSOCGroup <- data.frame(BRCA = c(rep(x = "CNV", times = 17), rep(x = "None", times = (17306-17)), rep(x = "CNV", times = 93), rep(x = "None", times = (8679-93))), HGSOC = c(rep(x = "Ctrl", times = 17306), rep(x = "Case", times = 8679)))
table(OverallGroup)
table(HGSOCGroup)

chisq.test(table(OverallGroup))
chisq.test(table(HGSOCGroup))
prop.test(table(OverallGroup))
prop.test(t(table(OverallGroup)))
prop.test(t(table(HGSOCGroup)))
prop.test(t(table(OverallGroup)))$p.value
prop.test(t(table(HGSOCGroup)))$p.value
prop.test((table(OverallGroup)))$p.value
prop.test((table(HGSOCGroup)))$p.value
chisq.test(t(table(OverallGroup)))$p.value
chisq.test(t(table(HGSOCGroup)))$p.value


