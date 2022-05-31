#
# funciVar_importBED_functions_chromHMM.R
#

# importBEDs_chromHMM(opt$mark, hg38_seqlengths, hg38_chroms)
importBEDs_chromHMM <- function(markname = NULL, my_seqlengths = NULL, my_chroms = NULL) {
  if(is.null(markname)) {
    print("Input a mark name.")
    exit(1)
  }
  # mark name will be E1 to E7
  
  library(foreach)
  library(doParallel)
  library(dplyr)
  library(GenomicRanges)
  
  CNVs_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs"
  New_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
  
  allbeds = list.files(path = paste0(New_Folder, "/consensus_chromHMM_states/state_", markname), pattern = "*_hg38.bed$")
  allGRs <- vector(mode="list", length=length(allbeds))
  
  registerDoParallel(4)
  list_allGRs <- foreach(i=1:length(allbeds)) %dopar% {
    bedname = allbeds[[i]]
    sample_name = gsub(paste0("_", markname, ".*"), "", bedname)
    bedfile = read.delim(paste0(New_Folder, "/consensus_chromHMM_states/state_", markname, "/", bedname), header = F, stringsAsFactors = F)
    bedfile = bedfile[,1:3]
    names(bedfile) = c("chrom", "start", "end")
    if(!is.null(my_chroms)) {
      bedfile = bedfile[bedfile$chrom %in% my_chroms,] 
    }
    bedfile = mutate(bedfile, sample=sample_name, state=markname)
    gr <- makeGRangesFromDataFrame(bedfile, keep.extra.columns = T)
    if(!is.null(my_seqlengths)) {
      seqlengths(gr) <- replace(seqlengths(gr), intersect(names(seqlengths(gr)), names(my_seqlengths)), my_seqlengths[names(seqlengths(gr))])
    }
    print(paste0("processed ", bedname))
    return(gr)
  }
  
  print("processed all BEDs")
  
  allGRs = GRanges()
  for(single_gr in list_allGRs) {
    allGRs <- c(allGRs, single_gr)
  }
  return(allGRs)
}



