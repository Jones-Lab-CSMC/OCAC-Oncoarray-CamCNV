#!/usr/bin/env bash
#
#
# selectNumPatientsAndCNVsAtInputLocations.sh
#


New_Folder="/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder="/Users/devriesa/Documents/OncoArrayCNVs"


if [ $# -ne 3 ]
then
  echo "Usage: selectNumPatientsAndCNVsAtInputLocations.sh chrname startloc endloc (where chrname, startloc, and endloc are numeric)"
  exit 1
fi

chrname=$1
startloc=$2
endloc=$3


locbed=$(echo "scratch_"$chrname"_"$startloc"_"$endloc"_loc.bed")
scratchCNVsbed=$(echo "scratch_"$chrname"_"$startloc"_"$endloc"_overallcalls_CNVs.bed")
scratchCNVstxt=$(echo "scratch_"$chrname"_"$startloc"_"$endloc"_overallcalls_CNVs.txt")
scratchhgsocCNVsbed=$(echo "scratch_"$chrname"_"$startloc"_"$endloc"_hgsoccalls_CNVs.bed")
scratchhgsocCNVstxt=$(echo "scratch_"$chrname"_"$startloc"_"$endloc"_hgsoccalls_CNVs.txt")

################################
##### get location as .bed #####
################################

awk -v chr=$chrname -v start=$startloc -v end=$endloc 'BEGIN{OFS="\t"; print "chr"chr, start, end}' >| $locbed

head $locbed

###############################################
##### intersect probes/locations and CNVs #####
###############################################


bedtools intersect -wa -a $New_Folder/ocac_overall_cnvcalls_IDs_hg38.bed -b $locbed | awk 'BEGIN{OFS="\t"} !seen[$0]++' >| $scratchCNVsbed

bedtools intersect -wa -a $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.bed -b $locbed | awk 'BEGIN{OFS="\t"} !seen[$0]++' >| $scratchhgsocCNVsbed


awk 'BEGIN{OFS="\t"} {gsub(/;/, "\t")} {print}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.header $scratchCNVsbed | awk 'BEGIN{OFS="\t"} !seen[$0]++' >| $scratchCNVstxt


awk 'BEGIN{OFS="\t"} {gsub(/;/, "\t")} {print}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.header $scratchhgsocCNVsbed | awk 'BEGIN{OFS="\t"} !seen[$0]++' >| $scratchhgsocCNVstxt


#################
### All Cases ###
#################

##### get number overall cases in the gene #####
# print case dels
awk '($5=="del" && $9=="1"){print}' $scratchCNVstxt | awk -v numcasedels=0 'BEGIN{OFS=" "} (!seen[$6]++){numcasedels=numcasedels+1} END{print "total number overall case deletions is:", numcasedels}'

# print case dups
awk '($5=="dup" && $9=="1"){print}' $scratchCNVstxt | awk -v numcasedups=0 'BEGIN{OFS=" "} (!seen[$6]++){numcasedups=numcasedups+1} END{print "total number overall case duplications is:", numcasedups}'

# print case combined
awk '($9=="1"){print}' $scratchCNVstxt | awk -v numcaseboth=0 'BEGIN{OFS=" "} (!seen[$6]++){numcaseboth=numcaseboth+1} END{print "total number overall cases is:", numcaseboth}'


##### get number hgsoc cases in the gene #####

# print case dels
awk '($5=="del" && $9=="1"){print}' $scratchhgsocCNVstxt | awk -v numcasedels=0 'BEGIN{OFS=" "} (!seen[$6]++){numcasedels=numcasedels+1} END{print "total number hgsoc case deletions is:", numcasedels}'


# print case dups
awk '($5=="dup" && $9=="1"){print}' $scratchhgsocCNVstxt | awk -v numcasedups=0 'BEGIN{OFS=" "} (!seen[$6]++){numcasedups=numcasedups+1} END{print "total number overall case duplications is:", numcasedups}'


# print case combined
awk '($9=="1"){print}' $scratchhgsocCNVstxt | awk -v numcaseboth=0 'BEGIN{OFS=" "} (!seen[$6]++){numcaseboth=numcaseboth+1} END{print "total number hgsoc cases is:", numcaseboth}'



##### get number controls in the gene #####

# print control dels
awk '($5=="del" && $9=="0"){print}' $scratchCNVstxt | awk -v numcontroldels=0 'BEGIN{OFS=" "} (!seen[$6]++){numcontroldels=numcontroldels+1} END{print "total number control deletions is:", numcontroldels}'

# print control dups
awk '($5=="dup" && $9=="0"){print}' $scratchCNVstxt | awk -v numcontroldups=0 'BEGIN{OFS=" "} (!seen[$6]++){numcontroldups=numcontroldups+1} END{print "total number control duplications is:", numcontroldups}'


# print control combined
awk '($9=="0"){print}' $scratchCNVstxt | awk -v numcontrolboth=0 'BEGIN{OFS=" "} (!seen[$6]++){numcontrolboth=numcontrolboth+1} END{print "total number controls is:", numcontrolboth}'



rm $locbed
rm $scratchCNVsbed
rm $scratchCNVstxt
rm $scratchhgsocCNVsbed
rm $scratchhgsocCNVstxt
