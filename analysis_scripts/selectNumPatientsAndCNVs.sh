#!/usr/bin/env bash
#
#
# selectNumPatientsAndCNVs.sh
#

#pwd

New_Folder="/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder="/Users/devriesa/Documents/OncoArrayCNVs"


#this uses the by country CNVs and probes

if [ $# -ne 3 ]
then
  echo "Usage: selectNumPatientsAndCNVs.sh genename cnvtype analysistype (where genename is \"BRCA1\", cnvtype is \"del\", analysistype is \"overall\")"
  exit 1
fi

genename=$1
cnvtype=$2
analysistype=$3
#awk -v genesymbol=$genename '{print "You inputted the gene "genesymbol}'
#echo "You inputted the gene "$genename

genelocbed=$(echo "scratch_"$genename"_"$cnvtype"_"$analysistype"_loc.bed")
scratchCNVsbed=$(echo "scratch_"$genename"_"$cnvtype"_"$analysistype"_overallcalls_CNVs.bed")
scratchCNVstxt=$(echo "scratch_"$genename"_"$cnvtype"_"$analysistype"_overallcalls_CNVs.txt")
scratchhgsocCNVsbed=$(echo "scratch_"$genename"_"$cnvtype"_"$analysistype"_hgsoccalls_CNVs.bed")
scratchhgsocCNVstxt=$(echo "scratch_"$genename"_"$cnvtype"_"$analysistype"_hgsoccalls_CNVs.txt")



#########################################
##### get location at relevant gene #####
#########################################

awk -v genesymbol=$genename 'BEGIN{OFS="\t"} $4==genesymbol{print}' $CNVs_Folder/Biomart_All_Genes_Hg38.bed >| $genelocbed


###############################################
##### intersect probes/locations and CNVs #####
###############################################


bedtools intersect -wa -a $New_Folder/ocac_overall_cnvcalls_IDs_hg38.bed -b $genelocbed | awk 'BEGIN{OFS="\t"} !seen[$0]++' >| $scratchCNVsbed

bedtools intersect -wa -a $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.bed -b $genelocbed | awk 'BEGIN{OFS="\t"} !seen[$0]++' >| $scratchhgsocCNVsbed


awk 'BEGIN{OFS="\t"} {gsub(/;/, "\t")} {print}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.header $scratchCNVsbed | awk 'BEGIN{OFS="\t"} !seen[$0]++' >| $scratchCNVstxt

awk 'BEGIN{OFS="\t"} {gsub(/;/, "\t")} {print}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.header $scratchhgsocCNVsbed | awk 'BEGIN{OFS="\t"} !seen[$0]++' >| $scratchhgsocCNVstxt


##### get number overall cases in the gene #####
# print case dels
if [ $cnvtype = "del" ] && [ $analysistype = "overall" ]; then
  #statements
  awk '($5=="del" && $9=="1"){print}' $scratchCNVstxt | awk -v numcasedels=0 '(!seen[$6]++){numcasedels=numcasedels+1} END{print numcasedels}'
fi


# print case dups
if [ $cnvtype = "dup" ] && [ $analysistype = "overall" ]; then
  #statements
  awk '($5=="dup" && $9=="1"){print}' $scratchCNVstxt | awk -v numcasedups=0 '(!seen[$6]++){numcasedups=numcasedups+1} END{print numcasedups}'
fi


##### get number hgsoc cases in the gene #####


# print case dels
if [ $cnvtype = "del" ] && [ $analysistype = "hgsoc" ]; then
  #statements
  awk '($5=="del" && $9=="1"){print}' $scratchhgsocCNVstxt | awk -v numcasedels=0 '(!seen[$6]++){numcasedels=numcasedels+1} END{print numcasedels}'
fi


# print case dups
if [ $cnvtype = "dup" ] && [ $analysistype = "hgsoc" ]; then
  #statements
  awk '($5=="dup" && $9=="1"){print}' $scratchhgsocCNVstxt | awk -v numcasedups=0 '(!seen[$6]++){numcasedups=numcasedups+1} END{print numcasedups}'
fi





##### get number controls in the gene #####

# print control dels
if [ $cnvtype = "del" ] && [ $analysistype = "control" ]; then
  #statements
  awk '($5=="del" && $9=="0"){print}' $scratchCNVstxt | awk -v numcontroldels=0 '(!seen[$6]++){numcontroldels=numcontroldels+1} END{print numcontroldels}'
fi

# print control dups
if [ $cnvtype = "dup" ] && [ $analysistype = "control" ]; then
  #statements
  awk '($5=="dup" && $9=="0"){print}' $scratchCNVstxt | awk -v numcontroldups=0 '(!seen[$6]++){numcontroldups=numcontroldups+1} END{print numcontroldups}'
fi

# print control combined
if [ $cnvtype = "both" ] && [ $analysistype = "control" ]; then
  #statements
  awk '($9=="0"){print}' $scratchCNVstxt | awk -v numcontrolboth=0 '(!seen[$6]++){numcontrolboth=numcontrolboth+1} END{print numcontrolboth}'
fi


##### get number combined cases (dels and dups, "both") in the gene #####

# print case combined
if [ $cnvtype = "both" ] && [ $analysistype = "overall" ]; then
  #statements
  awk '($9=="1"){print}' $scratchCNVstxt | awk -v numcaseboth=0 '(!seen[$6]++){numcaseboth=numcaseboth+1} END{print numcaseboth}'
fi


# print hgsoc combined
if [ $cnvtype = "both" ] && [ $analysistype = "hgsoc" ]; then
  #statements
  awk '($9=="1"){print}' $scratchhgsocCNVstxt | awk -v numhgsocboth=0 '(!seen[$6]++){numhgsocboth=numhgsocboth+1} END{print numhgsocboth}'
fi


# done!


rm $genelocbed
rm $scratchCNVsbed
rm $scratchCNVstxt
rm $scratchhgsocCNVsbed
rm $scratchhgsocCNVstxt
