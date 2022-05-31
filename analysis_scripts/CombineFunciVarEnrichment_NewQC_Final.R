#!/usr/local/bin/Rscript
#
#
# CombineFunciVarEnrichment_NewQC_Final.R
#
#


library(dplyr)
library(GenomicRanges)
library(BiocGenerics)
library(ggplot2)

New_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs"

source(paste0(CNVs_Folder, "/funciVar_biofeaturebackground.R"))

setwd(New_Folder)



##############################################
### Overall ChromHMM Non-Coding Enrichment ###
##############################################

load("E1_ChromHMM_overall_sig05_noncoding_enrichment.RData")
E1_ChromHMM_overall_sig05_noncoding_enrichment <- enrichment
rownames(E1_ChromHMM_overall_sig05_noncoding_enrichment) <- paste0(E1_ChromHMM_overall_sig05_noncoding_enrichment$sample, "_", E1_ChromHMM_overall_sig05_noncoding_enrichment$state)

load("E2_ChromHMM_overall_sig05_noncoding_enrichment.RData")
E2_ChromHMM_overall_sig05_noncoding_enrichment <- enrichment
rownames(E2_ChromHMM_overall_sig05_noncoding_enrichment) <- paste0(E2_ChromHMM_overall_sig05_noncoding_enrichment$sample, "_", E2_ChromHMM_overall_sig05_noncoding_enrichment$state)

load("E3_ChromHMM_overall_sig05_noncoding_enrichment.RData")
E3_ChromHMM_overall_sig05_noncoding_enrichment <- enrichment
rownames(E3_ChromHMM_overall_sig05_noncoding_enrichment) <- paste0(E3_ChromHMM_overall_sig05_noncoding_enrichment$sample, "_", E3_ChromHMM_overall_sig05_noncoding_enrichment$state)

load("E4_ChromHMM_overall_sig05_noncoding_enrichment.RData")
E4_ChromHMM_overall_sig05_noncoding_enrichment <- enrichment
rownames(E4_ChromHMM_overall_sig05_noncoding_enrichment) <- paste0(E4_ChromHMM_overall_sig05_noncoding_enrichment$sample, "_", E4_ChromHMM_overall_sig05_noncoding_enrichment$state)

load("E5_ChromHMM_overall_sig05_noncoding_enrichment.RData")
E5_ChromHMM_overall_sig05_noncoding_enrichment <- enrichment
rownames(E5_ChromHMM_overall_sig05_noncoding_enrichment) <- paste0(E5_ChromHMM_overall_sig05_noncoding_enrichment$sample, "_", E5_ChromHMM_overall_sig05_noncoding_enrichment$state)

load("E6_ChromHMM_overall_sig05_noncoding_enrichment.RData")
E6_ChromHMM_overall_sig05_noncoding_enrichment <- enrichment
rownames(E6_ChromHMM_overall_sig05_noncoding_enrichment) <- paste0(E6_ChromHMM_overall_sig05_noncoding_enrichment$sample, "_", E6_ChromHMM_overall_sig05_noncoding_enrichment$state)

load("E7_ChromHMM_overall_sig05_noncoding_enrichment.RData")
E7_ChromHMM_overall_sig05_noncoding_enrichment <- enrichment
rownames(E7_ChromHMM_overall_sig05_noncoding_enrichment) <- paste0(E7_ChromHMM_overall_sig05_noncoding_enrichment$sample, "_", E7_ChromHMM_overall_sig05_noncoding_enrichment$state)

load("E8_ChromHMM_overall_sig05_noncoding_enrichment.RData")
E8_ChromHMM_overall_sig05_noncoding_enrichment <- enrichment
rownames(E8_ChromHMM_overall_sig05_noncoding_enrichment) <- paste0(E7_ChromHMM_overall_sig05_noncoding_enrichment$sample, "_", E8_ChromHMM_overall_sig05_noncoding_enrichment$state)

rm(enrichment)

all_ChromHMM_overall_sig05_noncoding_enrichment <- rbind(E1_ChromHMM_overall_sig05_noncoding_enrichment, E2_ChromHMM_overall_sig05_noncoding_enrichment, E3_ChromHMM_overall_sig05_noncoding_enrichment, E4_ChromHMM_overall_sig05_noncoding_enrichment, E5_ChromHMM_overall_sig05_noncoding_enrichment, E6_ChromHMM_overall_sig05_noncoding_enrichment, E7_ChromHMM_overall_sig05_noncoding_enrichment)
all_ChromHMM_overall_sig05_noncoding_enrichment_8state <- rbind(E1_ChromHMM_overall_sig05_noncoding_enrichment, E2_ChromHMM_overall_sig05_noncoding_enrichment, E3_ChromHMM_overall_sig05_noncoding_enrichment, E4_ChromHMM_overall_sig05_noncoding_enrichment, E5_ChromHMM_overall_sig05_noncoding_enrichment, E6_ChromHMM_overall_sig05_noncoding_enrichment, E7_ChromHMM_overall_sig05_noncoding_enrichment, E8_ChromHMM_overall_sig05_noncoding_enrichment)

etomarktype <- data.frame(state=c("E1", "E2", "E3", "E4", "E5", "E6", "E7"), mark=factor(c("Weak Promoter", "Active Promoter", "Active Region", "Active Enhancer", "Weak Enhancer", "Insulator", "Transcribed"), ordered = TRUE, levels = c("Weak Promoter", "Active Promoter", "Active Region", "Active Enhancer", "Weak Enhancer", "Insulator", "Transcribed")))
all_ChromHMM_overall_sig05_noncoding_enrichment <- merge(x = all_ChromHMM_overall_sig05_noncoding_enrichment, y = etomarktype, by = "state", all.x = TRUE)
all_ChromHMM_overall_sig05_noncoding_enrichment <- mutate(.data = all_ChromHMM_overall_sig05_noncoding_enrichment, histotype=sub(pattern = "_consensus", replacement = "", x = sample))
write.table(x = all_ChromHMM_overall_sig05_noncoding_enrichment, file = "funcivar_ChromHMM_overall_sig05_noncoding_enrichment_cnvs.txt", quote = FALSE, sep = "\t", row.names = F, col.names = T)

etomarktype_8state <- data.frame(state=c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8"), mark=factor(c("Weak Promoter", "Active Promoter", "Active Region", "Active Enhancer", "Weak Enhancer", "Insulator", "Transcribed", "Low Signal"), ordered = TRUE, levels = c("Weak Promoter", "Active Promoter", "Active Region", "Active Enhancer", "Weak Enhancer", "Insulator", "Transcribed", "Low Signal")))
all_ChromHMM_overall_sig05_noncoding_enrichment_8state <- merge(x = all_ChromHMM_overall_sig05_noncoding_enrichment_8state, y = etomarktype_8state, by = "state", all.x = TRUE)
all_ChromHMM_overall_sig05_noncoding_enrichment_8state <- mutate(.data = all_ChromHMM_overall_sig05_noncoding_enrichment_8state, histotype=sub(pattern = "_consensus", replacement = "", x = sample))


# formatted but normal font

variant.enrichment = all_ChromHMM_overall_sig05_noncoding_enrichment
graphtitle = NULL
jpgname = "ChromHMM_overall_sig05_noncoding_enrichment_celltype_formatted.jpg"
filepath = "Figures/Bubble_Plots"
plot_height = 3.8
plot_width = 7.5

filter(all_ChromHMM_overall_sig05_noncoding_enrichment, difference==min(all_ChromHMM_overall_sig05_noncoding_enrichment$difference))

all_ChromHMM_overall_sig05_noncoding_enrichment$sample
variant.enrichment$histotype <- factor(variant.enrichment$histotype, levels = c("EEC", "IOSE", "FT", "HGSOC", "LGSOC", "CCOC", "MOC"), labels = c("EEC", "IOSE", "FT", "HGSOC", "LGSOC", "CCOC", "MOC"))



formattedplot <- 
  ggplot(data = variant.enrichment, aes(x = histotype, y = mark)) +
  geom_point(aes(size = probability, color = difference)) +
  scale_y_discrete(limits = rev(levels(variant.enrichment$mark))) +
  geom_point(data = variant.enrichment[variant.enrichment$significant == T,], color = "black", aes(shape = "a", stroke = 2, size = probability)) +
  scale_shape_discrete(solid=F, name = "Significant", labels = "true") +
  theme_bw() +
  scale_color_distiller(direction = -1, palette = "RdBu", limits = c(-1,1)*max(abs(variant.enrichment$difference))) + 
  theme(legend.position="right", legend.key.size = unit(.15, "inches"),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(1)),
        strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = rel(1.4)), 
        strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 1, size = rel(1.8)),
        strip.background = element_blank(),
        panel.spacing = unit(.1, "lines")) +
  ggtitle(graphtitle) + 
  labs(x = "Histotype Group", y = "Biofeature", size = "Probability", color = "Difference")

formattedplot

ggsave(filename = jpgname,
       plot = formattedplot, #+ coord_fixed(),
       device = "jpeg",
       width = plot_width,
       height = plot_height,
       dpi = 600,
       path = filepath)


############################################
### HGSOC ChromHMM Non-Coding Enrichment ###
############################################

load("E1_ChromHMM_hgsoc_sig05_noncoding_enrichment.RData")
E1_ChromHMM_hgsoc_sig05_noncoding_enrichment <- enrichment
rownames(E1_ChromHMM_hgsoc_sig05_noncoding_enrichment) <- paste0(E1_ChromHMM_hgsoc_sig05_noncoding_enrichment$sample, "_", E1_ChromHMM_hgsoc_sig05_noncoding_enrichment$state)

load("E2_ChromHMM_hgsoc_sig05_noncoding_enrichment.RData")
E2_ChromHMM_hgsoc_sig05_noncoding_enrichment <- enrichment
rownames(E2_ChromHMM_hgsoc_sig05_noncoding_enrichment) <- paste0(E2_ChromHMM_hgsoc_sig05_noncoding_enrichment$sample, "_", E2_ChromHMM_hgsoc_sig05_noncoding_enrichment$state)

load("E3_ChromHMM_hgsoc_sig05_noncoding_enrichment.RData")
E3_ChromHMM_hgsoc_sig05_noncoding_enrichment <- enrichment
rownames(E3_ChromHMM_hgsoc_sig05_noncoding_enrichment) <- paste0(E3_ChromHMM_hgsoc_sig05_noncoding_enrichment$sample, "_", E3_ChromHMM_hgsoc_sig05_noncoding_enrichment$state)

load("E4_ChromHMM_hgsoc_sig05_noncoding_enrichment.RData")
E4_ChromHMM_hgsoc_sig05_noncoding_enrichment <- enrichment
rownames(E4_ChromHMM_hgsoc_sig05_noncoding_enrichment) <- paste0(E4_ChromHMM_hgsoc_sig05_noncoding_enrichment$sample, "_", E4_ChromHMM_hgsoc_sig05_noncoding_enrichment$state)

load("E5_ChromHMM_hgsoc_sig05_noncoding_enrichment.RData")
E5_ChromHMM_hgsoc_sig05_noncoding_enrichment <- enrichment
rownames(E5_ChromHMM_hgsoc_sig05_noncoding_enrichment) <- paste0(E5_ChromHMM_hgsoc_sig05_noncoding_enrichment$sample, "_", E5_ChromHMM_hgsoc_sig05_noncoding_enrichment$state)

load("E6_ChromHMM_hgsoc_sig05_noncoding_enrichment.RData")
E6_ChromHMM_hgsoc_sig05_noncoding_enrichment <- enrichment
rownames(E6_ChromHMM_hgsoc_sig05_noncoding_enrichment) <- paste0(E6_ChromHMM_hgsoc_sig05_noncoding_enrichment$sample, "_", E6_ChromHMM_hgsoc_sig05_noncoding_enrichment$state)

load("E7_ChromHMM_hgsoc_sig05_noncoding_enrichment.RData")
E7_ChromHMM_hgsoc_sig05_noncoding_enrichment <- enrichment
rownames(E7_ChromHMM_hgsoc_sig05_noncoding_enrichment) <- paste0(E7_ChromHMM_hgsoc_sig05_noncoding_enrichment$sample, "_", E7_ChromHMM_hgsoc_sig05_noncoding_enrichment$state)

rm(enrichment)

all_ChromHMM_hgsoc_sig05_noncoding_enrichment <- rbind(E1_ChromHMM_hgsoc_sig05_noncoding_enrichment, E2_ChromHMM_hgsoc_sig05_noncoding_enrichment, E3_ChromHMM_hgsoc_sig05_noncoding_enrichment, E4_ChromHMM_hgsoc_sig05_noncoding_enrichment, E5_ChromHMM_hgsoc_sig05_noncoding_enrichment, E6_ChromHMM_hgsoc_sig05_noncoding_enrichment, E7_ChromHMM_hgsoc_sig05_noncoding_enrichment)

head(all_ChromHMM_hgsoc_sig05_noncoding_enrichment)

etomarktype <- data.frame(state=c("E1", "E2", "E3", "E4", "E5", "E6", "E7"), mark=factor(c("Weak Promoter", "Active Promoter", "Active Region", "Active Enhancer", "Weak Enhancer", "Insulator", "Transcribed"), ordered = TRUE, levels = c("Weak Promoter", "Active Promoter", "Active Region", "Active Enhancer", "Weak Enhancer", "Insulator", "Transcribed")))
all_ChromHMM_hgsoc_sig05_noncoding_enrichment <- merge(x = all_ChromHMM_hgsoc_sig05_noncoding_enrichment, y = etomarktype, by = "state", all.x = TRUE)
all_ChromHMM_hgsoc_sig05_noncoding_enrichment <- mutate(.data = all_ChromHMM_hgsoc_sig05_noncoding_enrichment, histotype=sub(pattern = "_consensus", replacement = "", x = sample))
write.table(x = all_ChromHMM_hgsoc_sig05_noncoding_enrichment, file = "funcivar_ChromHMM_hgsoc_sig05_noncoding_enrichment_cnvs.txt", quote = FALSE, sep = "\t", row.names = F, col.names = T)


# formatted but normal font

variant.enrichment = all_ChromHMM_hgsoc_sig05_noncoding_enrichment
graphtitle = NULL
jpgname = "ChromHMM_hgsoc_sig05_noncoding_enrichment_celltype_formatted.jpg"
filepath = "Figures/Bubble_Plots"
plot_height = 3.8
plot_width = 7.5

variant.enrichment$histotype <- factor(variant.enrichment$histotype, levels = c("EEC", "IOSE", "FT", "HGSOC", "LGSOC", "CCOC", "MOC"), labels = c("EEC", "IOSE", "FT", "HGSOC", "LGSOC", "CCOC", "MOC"))


formattedplot <- 
  ggplot(data = variant.enrichment, aes(x = histotype, y = mark)) +
  geom_point(aes(size = probability, color = difference)) +
  scale_y_discrete(limits = rev(levels(variant.enrichment$mark))) +
  geom_point(data = variant.enrichment[variant.enrichment$significant == T,], color = "black", aes(shape = "a", stroke = 2, size = probability)) +
  scale_shape_discrete(solid=F, name = "Significant", labels = "true") +
  theme_bw() +
  scale_color_distiller(direction = -1, palette = "RdBu", limits = c(-1,1)*max(abs(variant.enrichment$difference))) + 
  theme(legend.position="right", legend.key.size = unit(.15, "inches"),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(1)),
        strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = rel(1.4)),
        strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 1, size = rel(1.8)),
        strip.background = element_blank(),
        panel.spacing = unit(.1, "lines")) +
  ggtitle(graphtitle) + 
  labs(x = "Histotype Group", y = "Biofeature", size = "Probability", color = "Difference")

formattedplot

ggsave(filename = jpgname,
       plot = formattedplot, #+ coord_fixed(),
       device = "jpeg",
       width = plot_width,
       height = plot_height,
       dpi = 600,
       path = filepath)



#################################
### Make New Enrichment Table ###
#################################


# sample	AnalysisGroup	fg.ratio	bg.ratio	probability	difference	fg.success	fg.total	bg.success	bg.total	significant
all_ChromHMM_overall_sig05_noncoding_enrichment
all_ChromHMM_hgsoc_sig05_noncoding_enrichment

combinedSuppTable <- rbind(mutate(all_ChromHMM_overall_sig05_noncoding_enrichment, set = "Overall"), mutate(all_ChromHMM_hgsoc_sig05_noncoding_enrichment, set = "HGSOC")) %>% select(., set, state, sample, fg.ratio, bg.ratio, probability, difference, fg.success, fg.total, bg.success, bg.total, significant)

combinedSuppTable[which(combinedSuppTable$state=="E1"), 'state'] <- "Weak_Promoter"
combinedSuppTable[which(combinedSuppTable$state=="E2"), 'state'] <- "Active_Promoter"
combinedSuppTable[which(combinedSuppTable$state=="E3"), 'state'] <- "Active_Region"
combinedSuppTable[which(combinedSuppTable$state=="E4"), 'state'] <- "Active_Enhancer"
combinedSuppTable[which(combinedSuppTable$state=="E5"), 'state'] <- "Weak_Enhancer"
combinedSuppTable[which(combinedSuppTable$state=="E6"), 'state'] <- "Insulator"
combinedSuppTable[which(combinedSuppTable$state=="E7"), 'state'] <- "Transcribed"

write.table(x = combinedSuppTable, file = "Tables/Funcivar_ChromHMM_Enrichment_Table.txt", sep = "\t", row.names = F, quote = F, col.names = T)
