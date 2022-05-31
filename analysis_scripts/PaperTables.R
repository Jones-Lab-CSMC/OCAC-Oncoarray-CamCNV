# PaperTables.R

library(dplyr)
library(qqman)

New_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder <- "/Users/devriesa/Documents/OncoArrayCNVs"

phenotypes <- read.table(file = paste0(New_Folder, "/oncid_phenotypes_selected_samples_20200910.txt"), header = T, sep = "\t")
overallcnvs <- read.table(file = paste0(New_Folder, "/ocac_overall_cnvcalls_hg38.txt"), header = T)
hgsoccnvs <- read.table(file = paste0(New_Folder, "/ocac_hgsoc_cnvcalls_hg38.txt"), header = T)

table(hgsoccnvs$serous_hg_extra)

###########
### AGE ###
###########

summary(filter(phenotypes, refage<125)$refage)
summary(filter(filter(phenotypes, overall==0), refage<125)$refage)
summary(filter(filter(phenotypes, overall==1), refage<125)$refage)
summary(filter(filter(phenotypes, serous_hg_extra==1), refage<125)$refage)

t.test(x = filter(filter(phenotypes, overall==0), refage<125)$refage, y = filter(filter(phenotypes, overall==1), refage<125)$refage)
t.test(x = filter(filter(phenotypes, overall==0), refage<125)$refage, y = filter(filter(phenotypes, serous_hg_extra==1), refage<125)$refage)

################
### Segments ###
################

nrow(overallcnvs)

nrow(overallcnvs)
nrow(filter(overallcnvs, del_dup=="del"))
nrow(filter(overallcnvs, del_dup=="del"))/nrow(overallcnvs)
nrow(filter(overallcnvs, del_dup=="dup"))
nrow(filter(overallcnvs, del_dup=="dup"))/nrow(overallcnvs)

nrow(filter(overallcnvs, case_control==0))
nrow(filter(overallcnvs, case_control==1))

nrow(filter(filter(overallcnvs, case_control==0), del_dup=="del"))
nrow(filter(filter(overallcnvs, case_control==0), del_dup=="del"))/nrow(filter(overallcnvs, case_control==0))
nrow(filter(filter(overallcnvs, case_control==0), del_dup=="dup"))
nrow(filter(filter(overallcnvs, case_control==0), del_dup=="dup"))/nrow(filter(overallcnvs, case_control==0))

nrow(filter(filter(overallcnvs, case_control==1), del_dup=="del"))
nrow(filter(filter(overallcnvs, case_control==1), del_dup=="del"))/nrow(filter(overallcnvs, case_control==1))
nrow(filter(filter(overallcnvs, case_control==1), del_dup=="dup"))
nrow(filter(filter(overallcnvs, case_control==1), del_dup=="dup"))/nrow(filter(overallcnvs, case_control==1))

nrow(hgsoccnvs)
nrow(filter(hgsoccnvs, serous_hg_extra==0)) # this is different than above
nrow(filter(hgsoccnvs, serous_hg_extra==1))

nrow(filter(filter(hgsoccnvs, serous_hg_extra==1), del_dup=="del"))
nrow(filter(filter(hgsoccnvs, serous_hg_extra==1), del_dup=="del"))/nrow(filter(hgsoccnvs, serous_hg_extra==1))
nrow(filter(filter(hgsoccnvs, serous_hg_extra==1), del_dup=="dup"))
nrow(filter(filter(hgsoccnvs, serous_hg_extra==1), del_dup=="dup"))/nrow(filter(hgsoccnvs, serous_hg_extra==1))


# how many CNVs per person
controlids <- unique(filter(phenotypes, overall==0)$onc_id)
overallids <- unique(filter(phenotypes, overall==1)$onc_id)
hgsocids <- unique(filter(phenotypes, serous_hg_extra==1)$onc_id)

cnvsPerControl <- filter(overallcnvs, OncID %in% controlids) %>% group_by(OncID) %>% summarize(., count = n())
allcontrols <- data.frame(IDs = controlids, count = 0)
for(id in allcontrols$IDs) {
  if(id %in% cnvsPerControl$OncID) {
    allcontrols[which(allcontrols$IDs==id), 'count'] = cnvsPerControl[which(cnvsPerControl$OncID==id), 'count']
  }
}
summary(allcontrols$count)


cnvsPerOverall <- filter(overallcnvs, OncID %in% overallids) %>% group_by(OncID) %>% summarize(., count = n())
alloverall <- data.frame(IDs = overallids, count = 0)
for(id in alloverall$IDs) {
  if(id %in% cnvsPerOverall$OncID) {
    alloverall[which(alloverall$IDs==id), 'count'] = cnvsPerOverall[which(cnvsPerOverall$OncID==id), 'count']
  }
}
summary(alloverall$count)
t.test(x = allcontrols$count, y = alloverall$count)

cnvsPerHGSOC <- filter(overallcnvs, OncID %in% hgsocids) %>% group_by(OncID) %>% summarize(., count = n())
allhgsoc <- data.frame(IDs = hgsocids, count = 0)
for(id in allhgsoc$IDs) {
  if(id %in% cnvsPerHGSOC$OncID) {
    allhgsoc[which(allhgsoc$IDs==id), 'count'] = cnvsPerHGSOC[which(cnvsPerHGSOC$OncID==id), 'count']
  }
}
summary(allhgsoc$count)
t.test(x = allcontrols$count, y = allhgsoc$count)



### Dels ###
cnvsPerControlDels <- filter(overallcnvs, OncID %in% controlids) %>% filter(., del_dup=="del") %>% group_by(OncID) %>% summarize(., count = n())
allcontrolsdels <- data.frame(IDs = controlids, count = 0)
for(id in allcontrolsdels$IDs) {
  if(id %in% cnvsPerControlDels$OncID) {
    allcontrolsdels[which(allcontrolsdels$IDs==id), 'count'] = cnvsPerControlDels[which(cnvsPerControlDels$OncID==id), 'count']
  }
}
summary(allcontrolsdels$count)


cnvsPerOverallDels <- filter(overallcnvs, OncID %in% overallids) %>% filter(., del_dup=="del") %>% group_by(OncID) %>% summarize(., count = n())
alloveralldels <- data.frame(IDs = overallids, count = 0)
for(id in alloveralldels$IDs) {
  if(id %in% cnvsPerOverallDels$OncID) {
    alloveralldels[which(alloveralldels$IDs==id), 'count'] = cnvsPerOverallDels[which(cnvsPerOverallDels$OncID==id), 'count']
  }
}
summary(alloveralldels$count)
t.test(x = allcontrolsdels$count, y = alloveralldels$count)

cnvsPerHGSOCDels <- filter(overallcnvs, OncID %in% hgsocids) %>% filter(., del_dup=="del") %>% group_by(OncID) %>% summarize(., count = n())
allhgsocdels <- data.frame(IDs = hgsocids, count = 0)
for(id in allhgsocdels$IDs) {
  if(id %in% cnvsPerHGSOCDels$OncID) {
    allhgsocdels[which(allhgsocdels$IDs==id), 'count'] = cnvsPerHGSOCDels[which(cnvsPerHGSOCDels$OncID==id), 'count']
  }
}
summary(allhgsocdels$count)
t.test(x = allcontrolsdels$count, y = allhgsocdels$count)




### Dups ###
cnvsPerControlDups <- filter(overallcnvs, OncID %in% controlids) %>% filter(., del_dup=="dup") %>% group_by(OncID) %>% summarize(., count = n())
allcontrolsdups <- data.frame(IDs = controlids, count = 0)
for(id in allcontrolsdups$IDs) {
  if(id %in% cnvsPerControlDups$OncID) {
    allcontrolsdups[which(allcontrolsdups$IDs==id), 'count'] = cnvsPerControlDups[which(cnvsPerControlDups$OncID==id), 'count']
  }
}
summary(allcontrolsdups$count)


cnvsPerOverallDups <- filter(overallcnvs, OncID %in% overallids) %>% filter(., del_dup=="dup") %>% group_by(OncID) %>% summarize(., count = n())
alloveralldups <- data.frame(IDs = overallids, count = 0)
for(id in alloveralldups$IDs) {
  if(id %in% cnvsPerOverallDups$OncID) {
    alloveralldups[which(alloveralldups$IDs==id), 'count'] = cnvsPerOverallDups[which(cnvsPerOverallDups$OncID==id), 'count']
  }
}
summary(alloveralldups$count)
t.test(x = allcontrolsdups$count, y = alloveralldups$count)

cnvsPerHGSOCDups <- filter(overallcnvs, OncID %in% hgsocids) %>% filter(., del_dup=="dup") %>% group_by(OncID) %>% summarize(., count = n())
allhgsocdups <- data.frame(IDs = hgsocids, count = 0)
for(id in allhgsocdups$IDs) {
  if(id %in% cnvsPerHGSOCDups$OncID) {
    allhgsocdups[which(allhgsocdups$IDs==id), 'count'] = cnvsPerHGSOCDups[which(cnvsPerHGSOCDups$OncID==id), 'count']
  }
}
summary(allhgsocdups$count)
t.test(x = allcontrolsdups$count, y = allhgsocdups$count)



### CNV Length ###

summary(overallcnvs$length)
summary(filter(overallcnvs, OncID %in% controlids)$length)
summary(filter(overallcnvs, OncID %in% overallids)$length)
summary(filter(overallcnvs, OncID %in% hgsocids)$length)
t.test(x = filter(overallcnvs, OncID %in% controlids)$length, y = filter(overallcnvs, OncID %in% overallids)$length)
t.test(x = filter(overallcnvs, OncID %in% controlids)$length, y = filter(overallcnvs, OncID %in% hgsocids)$length)

# dels
summary(filter(overallcnvs, del_dup=="del")$length)
summary(filter(filter(overallcnvs, OncID %in% controlids), del_dup=="del")$length)
summary(filter(filter(overallcnvs, OncID %in% overallids), del_dup=="del")$length)
summary(filter(filter(overallcnvs, OncID %in% hgsocids), del_dup=="del")$length)
t.test(x = filter(filter(overallcnvs, OncID %in% controlids), del_dup=="del")$length, y = filter(filter(overallcnvs, OncID %in% overallids), del_dup=="del")$length)
t.test(x = filter(filter(overallcnvs, OncID %in% controlids), del_dup=="del")$length, y = filter(filter(overallcnvs, OncID %in% hgsocids), del_dup=="del")$length)

# dups
summary(filter(overallcnvs, del_dup=="dup")$length)
summary(filter(filter(overallcnvs, OncID %in% controlids), del_dup=="dup")$length)
summary(filter(filter(overallcnvs, OncID %in% overallids), del_dup=="dup")$length)
summary(filter(filter(overallcnvs, OncID %in% hgsocids), del_dup=="dup")$length)
t.test(x = filter(filter(overallcnvs, OncID %in% controlids), del_dup=="dup")$length, y = filter(filter(overallcnvs, OncID %in% overallids), del_dup=="dup")$length)
t.test(x = filter(filter(overallcnvs, OncID %in% controlids), del_dup=="dup")$length, y = filter(filter(overallcnvs, OncID %in% hgsocids), del_dup=="dup")$length)

table(phenotypes$casecon)
table(phenotypes$overall)
table(phenotypes$overall_all)
table(phenotypes$mucinous)
table(phenotypes$mucinous_all)
table(phenotypes$ser_lg_lmp)
table(phenotypes$serouslowgrade)

table(filter(phenotypes, set_country=="")$set_centre)
# table(filter(phenotypes, set_centre=="UKR")$set_country)
for(i in which(phenotypes$set_country=="")) {
  if(phenotypes[i, 'set_centre']=="BGS") {
    phenotypes[i, 'set_country']="UK"
  }
  if(phenotypes[i, 'set_centre']=="CAM") {
    phenotypes[i, 'set_country']="UK"
  }
  if(phenotypes[i, 'set_centre']=="EPC") {
    phenotypes[i, 'set_country']="Europe"
  }
  if(phenotypes[i, 'set_centre']=="MAY") {
    phenotypes[i, 'set_country']="USA"
  }
  if(phenotypes[i, 'set_centre']=="MOF") {
    phenotypes[i, 'set_country']="USA"
  }
  if(phenotypes[i, 'set_centre']=="NCO") {
    phenotypes[i, 'set_country']="USA"
  }
  if(phenotypes[i, 'set_centre']=="POL") {
    phenotypes[i, 'set_country']="Poland"
  }
  if(phenotypes[i, 'set_centre']=="RMH") {
    phenotypes[i, 'set_country']="UK"
  }
  if(phenotypes[i, 'set_centre']=="SEA") {
    phenotypes[i, 'set_country']="UK"
  }
  if(phenotypes[i, 'set_centre']=="TBO") {
    phenotypes[i, 'set_country']="USA"
  }
  if(phenotypes[i, 'set_centre']=="UKO") {
    phenotypes[i, 'set_country']="UK"
  }
  if(phenotypes[i, 'set_centre']=="UKR") {
    phenotypes[i, 'set_country']="UK"
  }
}
table(filter(phenotypes, set_country=="")$set_centre)
table(phenotypes$set_country)


table(filter(phenotypes, is.na(set_country))$set_centre)
table(filter(phenotypes, set_centre=="MOF")$set_country)
for(i in which(is.na(phenotypes$set_country))) {
  if(phenotypes[i, 'set_centre']=="MOF") {
    phenotypes[i, 'set_country']="USA"
  }
}
table(filter(phenotypes, is.na(set_country))$set_centre)
table(phenotypes$set_country)





newpheno <- select(phenotypes, Country=set_country, Control=overall, All=overall, LGSerous=ser_lg_lmp, HGSerous=serous_hg_extra, Mucinous=mucinous, Endometrioid=endometrioid, ClearCell=clearcell)

newpheno <- mutate(newpheno, Control=ifelse(test = Control==1, yes = 0, no = 1))
table(phenotypes$set_country)
table(newpheno$Country)
newpheno$Country <- factor(newpheno$Country)
levels(newpheno$Country)
sumpheno <- newpheno %>% dplyr::group_by(Country) %>% dplyr::summarise(., sum(Control), sum(All), sum(LGSerous, na.rm = T), sum(HGSerous, na.rm = T), sum(Mucinous, na.rm = T), sum(Endometrioid, na.rm = T), sum(ClearCell, na.rm = T))
names(sumpheno) <- c("Country", "Control", "All", "LGSerous", "HGSerous", "Mucinous", "Endometrioid", "ClearCell")
sumpheno <- mutate(sumpheno, Other=All-LGSerous-HGSerous-Mucinous-Endometrioid-ClearCell)
names(sumpheno) <- c("Country", "Control", "All", "LG serous", "HG serous", "Mucinous", "Endometrioid", "Clear cell", "Other")
sum(sumpheno$Control)
sum(sumpheno$All)

write.table(x = sumpheno, file = "Tables/CountrySampleCounts.txt", quote = F, sep = "\t", row.names = F, col.names = T)


# BRCA1 CNVs Distribution
BRCA1cnvs <- read.table(file = "BRCA1_CNVs_Info.txt", header = T, sep = "\t", na.strings = "NA")

BRCA1cnvs <- mutate(BRCA1cnvs, uniqueCNV=paste0(Chr, ":", Start, "-", End, "(", del_dup, ")"))
length(unique(BRCA1cnvs$uniqueCNV))
uniqueBRCA1 <- as.data.frame(table(BRCA1cnvs$uniqueCNV))
multiplePatientsBRCA1 <- filter(uniqueBRCA1, Freq!=1)
summary(as.numeric(table(BRCA1cnvs$uniqueCNV)))




#######################
### Authorship List ###
#######################

centercounts <- table(phenotypes$set_centre) %>% t() %>% t()
centercasecounts <- table(filter(phenotypes, overall==1)$set_centre) %>% t() %>% t()

center_allcounts <- data.frame(study = rownames(centercounts), samplecounts = centercounts[,1], casecounts = centercasecounts[,1])
View(center_allcounts)

write.table(x = center_allcounts, file = "set_centre_counts.txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)


table(select(phenotypes, set_centre, set_country)) %>% View()

for (country in unique(phenotypes$set_country)) {
  print(paste0("country ", country, " includes studies:"))
  print(paste0(unique(filter(phenotypes, set_country==country)$set_centre)))
}
