#!/usr/local/bin/Rscript
#
# FunciVar_Submit_sig05_CNVs_Overall_GWAS_TWAS_Loci_NewQC.R
#


# library(plyr)
library(GenomicRanges)
library(BiocGenerics)
library(ggplot2)
library(foreach)
library(doMC)
library(dplyr)
#https://www.r-bloggers.com/parallel-r-loops-for-windows-and-linux/
as.vector <- BiocGenerics::as.vector


New_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs"

source(paste0(CNVs_Folder, "/", "funciVar_biofeaturebackground.R"))

#import hg38 H3K27Ac or other mark
source(paste0(CNVs_Folder, "/", "funciVar_importBED_functions.R"))
# includes function importBEDs(markname)


#import seqlengths
hg38_sizes <- read.delim(paste0(CNVs_Folder, "/", "hg38.chrom.sizes.limited.txt"), header = F)
hg38_seqlengths <- hg38_sizes$V2
names(hg38_seqlengths) <- hg38_sizes$V1
hg38_chroms <- hg38_sizes$V1


GWAS = read.delim(paste0(New_Folder, "/", "GWAS_Loci_HRC_LD_unique_all_non_mucinous.txt"), header = T, stringsAsFactors = F, col.names = c("chrom", "start", "end"))

GWAS$chrom = paste0("chr", GWAS$chrom)
GWAS = mutate(GWAS, sample="GWAS_LD", state="GWAS_LD")
GWAS_noBRCA1 <- filter(GWAS, !(chrom=="chr17"&start==41743558&end==43694719))
GWAS_noBRCA1 <- mutate(GWAS_noBRCA1, sample="GWAS_noBRCA1_LD", state="GWAS_noBRCA1_LD")
gr_GWAS <- makeGRangesFromDataFrame(GWAS, keep.extra.columns = T, seqinfo = hg38_seqlengths)
seqlengths(gr_GWAS) <- replace(seqlengths(gr_GWAS), intersect(names(seqlengths(gr_GWAS)), names(hg38_seqlengths)), hg38_seqlengths[names(seqlengths(gr_GWAS))])
gr_GWAS_noBRCA1 <- makeGRangesFromDataFrame(GWAS_noBRCA1, keep.extra.columns = T, seqinfo = hg38_seqlengths)

TWAS2018gene = read.delim(paste0(CNVs_Folder, "/", "TWAS_Loci_2018_Gene_unique.txt"), header = T, stringsAsFactors = F, col.names = c("chrom", "start", "end"))
TWAS2018gene$chrom = paste0("chr", TWAS2018gene$chrom)
TWAS2018gene = mutate(TWAS2018gene, sample="TWAS_2018_Gene", state="TWAS_2018_Gene")
gr_TWAS_2018_Gene <- makeGRangesFromDataFrame(TWAS2018gene, keep.extra.columns = T, seqinfo = hg38_seqlengths)

TWAS2018geneld = read.delim(paste0(CNVs_Folder, "/", "TWAS_Loci_2018_Gene_LD_unique.txt"), header = T, stringsAsFactors = F, col.names = c("chrom", "start", "end"))
TWAS2018geneld$chrom = paste0("chr", TWAS2018geneld$chrom)
TWAS2018geneld = mutate(TWAS2018geneld, sample="TWAS_2018_Gene_LD", state="TWAS_2018_Gene_LD")
gr_TWAS_2018_Gene_LD <- makeGRangesFromDataFrame(TWAS2018geneld, keep.extra.columns = T, seqinfo = hg38_seqlengths)

TWAS2018ld = read.delim(paste0(CNVs_Folder, "/", "TWAS_Loci_2018_LD_unique.txt"), header = T, stringsAsFactors = F, col.names = c("chrom", "start", "end"))
TWAS2018ld$chrom = paste0("chr", TWAS2018ld$chrom)
TWAS2018ld = mutate(TWAS2018ld, sample="TWAS_2018_LD", state="TWAS_2018_LD")
gr_TWAS_2018_LD <- makeGRangesFromDataFrame(TWAS2018ld, keep.extra.columns = T, seqinfo = hg38_seqlengths)

TWAS2019gene = read.delim(paste0(CNVs_Folder, "/", "TWAS_Loci_2019_Gene_unique.txt"), header = T, stringsAsFactors = F, col.names = c("chrom", "start", "end"))
TWAS2019gene$chrom = paste0("chr", TWAS2019gene$chrom)
TWAS2019gene = mutate(TWAS2019gene, sample="TWAS_2019_Gene", state="TWAS_2019_Gene")
gr_TWAS_2019_Gene <- makeGRangesFromDataFrame(TWAS2019gene, keep.extra.columns = T, seqinfo = hg38_seqlengths)




TWAS2019geneld = read.delim(paste0(CNVs_Folder, "/", "TWAS_Loci_2019_Gene_LD_unique.txt"), header = T, stringsAsFactors = F, col.names = c("chrom", "start", "end"))
TWAS2019geneld$chrom = paste0("chr", TWAS2019geneld$chrom)
TWAS2019geneld = mutate(TWAS2019geneld, sample="TWAS_2019_Gene_LD", state="TWAS_2019_Gene_LD")
gr_TWAS_2019_Gene_LD <- makeGRangesFromDataFrame(TWAS2019geneld, keep.extra.columns = T, seqinfo = hg38_seqlengths)

TWAS2019ld = read.delim(paste0(CNVs_Folder, "/", "TWAS_Loci_2019_LD_unique.txt"), header = T, stringsAsFactors = F, col.names = c("chrom", "start", "end"))
TWAS2019ld$chrom = paste0("chr", TWAS2019ld$chrom)
TWAS2019ld = mutate(TWAS2019ld, sample="TWAS_2019_LD", state="TWAS_2019_LD")
gr_TWAS_2019_LD <- makeGRangesFromDataFrame(TWAS2019ld, keep.extra.columns = T, seqinfo = hg38_seqlengths)


# instructions in https://www.simoncoetzee.com/bioc2017.html

#converting variants into GRanges

fore <- read.delim("ocac_overall_sig05_CNVs_hg38_unique.txt", header = T, stringsAsFactors = F)

GRfore <- with(fore, GRanges(Chr, IRanges(startpos, endpos), cnv_type=del_dup, cnv_id=Chr0100000))
seqlengths(GRfore) <- replace(seqlengths(GRfore), intersect(names(seqlengths(GRfore)), names(hg38_seqlengths)), hg38_seqlengths[names(seqlengths(GRfore))])

back <- read.delim("ocac_overall_sig05_CNVs_hg38_sim1000_new.txt", header = T, stringsAsFactors = F, na.strings = "NA")

GRback <- with(back, GRanges(Chr, IRanges(startpos, endpos), cnv_type=del_dup, cnv_id=Chr0100000))
seqlengths(GRback) <- replace(seqlengths(GRback), intersect(names(seqlengths(GRback)), names(hg38_seqlengths)), hg38_seqlengths[names(seqlengths(GRback))])
GRback
GRback <- trim(GRback)

list_GR <- list(fg=GRfore, bg=GRback)

rm(back)
rm(GRback)


gr_GWAS
gr_TWAS_2018_Gene
gr_TWAS_2018_LD
gr_TWAS_2019_Gene
gr_TWAS_2019_LD
gr_TWAS_Gene <- union(gr_TWAS_2018_Gene, gr_TWAS_2019_Gene)
mcols(gr_TWAS_Gene) <- DataFrame(sample = rep("TWAS_Gene", length(gr_TWAS_Gene)), state = rep("TWAS_Gene", length(gr_TWAS_Gene)))
gr_TWAS_Gene_LD <- union(gr_TWAS_2018_Gene_LD, gr_TWAS_2019_Gene_LD)
mcols(gr_TWAS_Gene_LD) <- DataFrame(sample = rep("TWAS_Gene_LD", length(gr_TWAS_Gene_LD)), state = rep("TWAS_Gene_LD", length(gr_TWAS_Gene_LD)))
gr_TWAS_LD <- union(gr_TWAS_2018_LD, gr_TWAS_2019_LD)
mcols(gr_TWAS_LD) <- DataFrame(sample = rep("TWAS_LD", length(gr_TWAS_LD)), state = rep("TWAS_LD", length(gr_TWAS_LD)))
gr_TWAS <- union(gr_TWAS_Gene, gr_TWAS_LD)
mcols(gr_TWAS) <- DataFrame(sample = rep("TWAS", length(gr_TWAS)), state = rep("TWAS", length(gr_TWAS)))
gr_GWAS_TWAS <- list(gr_GWAS, gr_GWAS_noBRCA1, gr_TWAS_2018_Gene, gr_TWAS_2018_Gene_LD, gr_TWAS_2018_LD, gr_TWAS_2019_Gene, gr_TWAS_2019_LD, gr_TWAS_2019_Gene_LD, gr_TWAS_Gene, gr_TWAS_Gene_LD, gr_TWAS_LD, gr_TWAS)


GWAS_enrichment <- CalculateEnrichment(list_GR, gr_GWAS)
gwas_overlaps <- SetOverlaps(list_GR$fg, gr_GWAS)
gwas_cnvs_ld <- mergeByOverlaps(gwas_overlaps, gr_GWAS)
gwas_cnvs_ld_df = as.data.frame(DataFrame(Chr = gsub(pattern = "chr", replacement = "", x = as.character(seqnames(gwas_cnvs_ld$gwas_overlaps))), LD_start = start(gwas_cnvs_ld$gr_GWAS), LD_end = end(gwas_cnvs_ld$gr_GWAS), CNV_start = start(gwas_cnvs_ld$gwas_overlaps), CNV_end = end(gwas_cnvs_ld$gwas_overlaps), CNV_type = gwas_cnvs_ld$cnv_type))
head(gwas_cnvs_ld_df)


GWAS_loci_LD <- read.table(file = paste0(New_Folder, "/GWAS_Loci_HRC_coord_all_non_mucinous.txt"), header = T)
head(GWAS_loci_LD)

GWAS_loci_LD_sort <- GWAS_loci_LD[with(GWAS_loci_LD, order(GWAS_pvalue)),]
GWAS_loci_LD_unique <- GWAS_loci_LD_sort[!duplicated(GWAS_loci_LD_sort[c('LD_start', 'LD_end')]),]
GWAS_loci_LD_unique <- GWAS_loci_LD_unique[with(GWAS_loci_LD_unique, order(Chr, LD_start)),]
head(GWAS_loci_LD_unique)


GWAS_loci_LD_unique_overlap <- mutate(GWAS_loci_LD_unique, Chr_Pos = paste(Chr, GWAS_end, sep = "_"))
head(GWAS_loci_LD_unique_overlap)

GWAS_coords <- read.table(file = paste0(New_Folder, "/GWAS_Loci_HRC_coord_all_non_mucinous.txt"), header = T, sep = "\t")
GWAS_coords <- mutate(GWAS_coords, Chr_Pos = paste(Chr, GWAS_end, sep = "_"))
GWAS_coords <- mutate(GWAS_coords, Color = "0,0,0,")
head(GWAS_coords)

# black: 0,0,0,
# purple: 102,0,204

for (row in 1:nrow(GWAS_coords)) {
  if (GWAS_coords[row, 'Chr_Pos'] %in% GWAS_loci_LD_unique_overlap$Chr_Pos) {
    GWAS_coords[row, 'Color'] = "102,0,204,"
  }
}
table(GWAS_coords$Color)

GWAS_hits_bed <- data.frame(paste0("chr", GWAS_coords$Chr), GWAS_coords$GWAS_start, GWAS_coords$GWAS_end, GWAS_coords$GWAS_pvalue, 0, ".", GWAS_coords$GWAS_start, GWAS_coords$GWAS_end, GWAS_coords$Color)
head(GWAS_hits_bed)
write.table(x = GWAS_hits_bed, file = "HRC_GWAS_hits_hg38_sig05_NewQC.bed", sep = "\t", col.names = F, row.names = F, quote = F)

# must manually re-add track title:
# made bed of GWAS hits
# ~/Documents/OncoArrayCNVs/UCSC_hg38_BED_Files/HRC_GWAS_hits_hg38.bed
# manually added track name: track name="HRC_GWAS_hits_hg38" description="All Genome-Wide Significance GWAS Hits for HRC Imputation, Black is p<5e-8, Purple is Most Significant GWAS Hit in an LD Block, name is p-value" itemRgb="on"







all_cnvs_ld_gwasSNP <- merge(gwas_cnvs_ld_df, GWAS_loci_LD_unique, by = c("Chr", "LD_start", "LD_end"))
write.table(x = all_cnvs_ld_gwasSNP, file = "CNVs_at_GWAS_Loci_with_LD_SNP_coords_sig05_NewQC.txt", quote = F, sep = "\t", row.names = F, col.names = T)
only_ld_gwasSNP <- unique(select(all_cnvs_ld_gwasSNP, -CNV_start, -CNV_end, -CNV_type))
head(only_ld_gwasSNP)
# first row, chr17  41743558  43694719, contains all BRCA1 CNVs
write.table(x = only_ld_gwasSNP, file = "CNVs_at_GWAS_Loci_only_LD_SNP_coords_sig05_NewQC.txt", quote = F, sep = "\t", row.names = F, col.names = T)

gwas_overlaps_df <- mutate(as.data.frame(mcols(gwas_overlaps)), chr = as.character(seqnames(gwas_overlaps)), start = start(gwas_overlaps), end = end(gwas_overlaps))
gwas_overlap_cnvs <- select(filter(gwas_overlaps_df, GWAS_LD==1), chr, start, end, cnv_type, cnv_id)
head(gwas_overlap_cnvs)
write.table(x = gwas_overlap_cnvs, file = "CNVs_at_GWAS_Loci_sig05_NewQC.txt", quote = F, sep = "\t", row.names = F, col.names = T)



registerDoMC(2)
list_enrichments <- foreach(i=1:length(gr_GWAS_TWAS)) %dopar% {
  cur_gr <- gr_GWAS_TWAS[[i]]
  temp_enrichment <- NULL
  temp_enrichment <- CalculateEnrichment(list_GR, cur_gr)
  return(temp_enrichment)
}

enrichment = data.frame()
for(single_enrichment in list_enrichments) {
  enrichment <- rbind(enrichment, single_enrichment)
}

write.table(x = enrichment, file = "GWAS_TWAS_enrichment_results_sig05_NewQC.txt", quote = F, sep = "\t", row.names = F, col.names = T)


enrichment_short <- select(enrichment, sample, fg.ratio, bg.ratio, probability, difference, fg.success, fg.total, bg.success, bg.total, significant)

write.table(x = enrichment_short, file = "GWAS_TWAS_enrichment_results_short_sig05_NewQC.txt", quote = F, sep = "\t", row.names = F, col.names = T)

