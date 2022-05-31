#!/usr/local/bin/Rscript
#
# FunciVar_Submit_NonCoding_sig05_CNVs_HGSOC_ChromHMM.R
#

print("successfully entered R")
library(optparse)

my_options = list(
  make_option(c("-m", "--mark"), action="store", default=NA, type='character',
              help="A variable storing the name of the mark you're using.")
  #http://www.cureffi.org/2014/01/15/running-r-batch-mode-linux/
  #https://cran.r-project.org/web/packages/optparse/vignettes/optparse.html
)
opt = parse_args(OptionParser(option_list=my_options))

print("parsed args")
if(!is.na(opt$mark)) {
  print(paste0("you inputted the mark: ", opt$mark))

  
  library(plyr)
  library(GenomicRanges)
  library(BiocGenerics)
  library(ggplot2)
  as.vector <- BiocGenerics::as.vector
  library(foreach)
  library(doParallel)
  
  New_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
  CNVs_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs"
  
  source(paste0(CNVs_Folder, "/", "funciVar_biofeaturebackground.R"))
  
  #import hg38 H3K27Ac or other mark
  source(paste0(New_Folder, "/", "funciVar_importBED_functions_chromHMM.R"))
  # includes function importBEDs(markname)
  
  setwd(New_Folder)
  
  #import seqlengths
  hg38_sizes <- read.delim(paste0(CNVs_Folder, "/", "hg38.chrom.sizes.limited.txt"), header = F)
  hg38_seqlengths <- hg38_sizes$V2
  names(hg38_seqlengths) <- hg38_sizes$V1
  hg38_chroms <- hg38_sizes$V1
  
  mark_features <- importBEDs_chromHMM(opt$mark, hg38_seqlengths, hg38_chroms)
  
  samples <- unique(mark_features$sample)
  
  # instructions in https://www.simoncoetzee.com/bioc2017.html
  
  #converting variants into GRanges
  fore <- read.delim(paste0(New_Folder, "/", "ocac_hgsoc_sig05_noncoding_CNVs_hg38.txt"), header = T, stringsAsFactors = F)
  
  GRfore <- with(fore, GRanges(Chr, IRanges(startpos, endpos), cnv_type=del_dup, cnv_id=Chr0100000))
  seqlengths(GRfore) <- replace(seqlengths(GRfore), intersect(names(seqlengths(GRfore)), names(hg38_seqlengths)), hg38_seqlengths[names(seqlengths(GRfore))])
  
  back <- read.delim(paste0(New_Folder, "/", "ocac_hgsoc_sig05_noncoding_CNVs_hg38_sim2000.txt"), header = T, stringsAsFactors = F, na.strings = "NA")
  
  GRback <- with(back, GRanges(Chr, IRanges(startpos, endpos), cnv_type=del_dup, cnv_id=Chr0100000))
  seqlengths(GRback) <- replace(seqlengths(GRback), intersect(names(seqlengths(GRback)), names(hg38_seqlengths)), hg38_seqlengths[names(seqlengths(GRback))])
  
  
  list_GR <- list(fg=GRfore, bg=GRback)
  
  rm(back)
  rm(GRback)
  
  if (length(samples)<4) {
    registerDoParallel(length(samples))
  } else {
    registerDoParallel(4)
  }
  list_enrichments <- foreach(i=1:length(samples)) %dopar% {
    cur_sample <- samples[i]
    cur_features <- subset(mark_features, mark_features$sample==cur_sample)
    temp_enrichment <- NULL
    temp_enrichment <- CalculateEnrichment(list_GR, cur_features)
    return(temp_enrichment)
  }
  
  enrichment = data.frame()
  for(single_enrichment in list_enrichments) {
    enrichment <- rbind(enrichment, single_enrichment)
  }
  
  enrichment <- mutate(enrichment, state=opt$mark)
  
  save(enrichment, file = paste0(New_Folder, "/", opt$mark, "_ChromHMM_hgsoc_sig05_noncoding_enrichment.RData"))
  
  
  
} else {
  cat("you didn't specify mark\n", file=stderr()) # print error messages to stderr
}