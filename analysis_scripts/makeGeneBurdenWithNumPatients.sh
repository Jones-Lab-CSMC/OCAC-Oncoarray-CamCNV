#!/usr/bin/env bash
#
#
# makeGeneBurdenWithNumPatients.sh
#

New_Folder="/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder="/Users/devriesa/Documents/OncoArrayCNVs"


awk 'BEGIN{OFS="\t"} (NR==1){gsub("\r", ""); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "cnv_call"; next} (FNR==NR){gsub("\r", ""); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "del"; next} (FNR!=NR && FNR!=1){gsub("\r", ""); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "dup"; next}' $New_Folder/raml_results_overall_dels.txt $New_Folder/raml_results_overall_dups.txt >| $New_Folder/raml_results_overall.txt

awk 'BEGIN{OFS="\t"} (NR==1){gsub("\r", ""); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "cnv_call"; next} (FNR==NR){gsub("\r", ""); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "del"; next} (FNR!=NR && FNR!=1){gsub("\r", ""); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "dup"; next}' $New_Folder/raml_results_hgsoc_dels.txt $New_Folder/raml_results_hgsoc_dups.txt >| $New_Folder/raml_results_hgsoc.txt

wc -l $New_Folder/raml_results_hgsoc.txt
wc -l $New_Folder/raml_results_hgsoc_both.txt
wc -l $New_Folder/raml_results_hgsoc_dels.txt
wc -l $New_Folder/raml_results_hgsoc_dups.txt
wc -l $New_Folder/raml_results_hgsoc_both.txt

wc -l $New_Folder/raml_results_overall.txt
for (( i = 1; i <= 3955; i++ )); do
  if [[ $i -eq 1 ]]; then
    awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, "NumCases", "NumControls", $16, $10, "AnalysisSet"; exit }' $New_Folder/raml_results_overall.txt
  fi
  gene=$(awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1; exit }' $New_Folder/raml_results_overall.txt)
  if [[ $gene == *\"* ]]
  then
    awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, "NA", "NA", $16, $10, "overall"; exit }' $New_Folder/raml_results_overall.txt
  elif [[ $i -ne 1 ]]
  then
    deldup=$(awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $16; exit }' $New_Folder/raml_results_overall.txt)
    casenum=$( $New_Folder/selectNumPatientsAndCNVs.sh $gene $deldup "overall" )
    ctrlnum=$( $New_Folder/selectNumPatientsAndCNVs.sh $gene $deldup "control" )
    awk -v row=$i -v cases=$casenum -v ctrls=$ctrlnum 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, cases, ctrls, $16, $10, "overall"; exit }' $New_Folder/raml_results_overall.txt
  fi
done >| $New_Folder/Tables/raml_results_overall_counts.txt

wc -l $New_Folder/raml_results_overall_both.txt
for (( i = 1; i <= 3897; i++ )); do
  if [[ $i -eq 1 ]]; then
    awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, "NumCases", "NumControls", "cnv_call", $10, "AnalysisSet"; exit }' $New_Folder/raml_results_overall_both.txt
  fi
  gene=$(awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1; exit }' $New_Folder/raml_results_overall_both.txt)
  if [[ $gene == *\"* ]]
  then
    awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, "NA", "NA", "both", $10, "overall"; exit }' $New_Folder/raml_results_overall_both.txt
  elif [[ $i -ne 1 ]]
  then
    casenum=$( $New_Folder/selectNumPatientsAndCNVs.sh $gene "both" "overall" )
    ctrlnum=$( $New_Folder/selectNumPatientsAndCNVs.sh $gene "both" "control" )
    awk -v row=$i -v cases=$casenum -v ctrls=$ctrlnum 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, cases, ctrls, "both", $10, "overall"; exit }' $New_Folder/raml_results_overall_both.txt
  fi
done >| $New_Folder/Tables/raml_results_overall_counts_both.txt

wc -l $New_Folder/raml_results_hgsoc.txt
for (( i = 1; i <= 3481; i++ )); do
  if [[ $i -eq 1 ]]; then
    awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, "NumCases", "NumControls", $16, $10, "AnalysisSet"; exit }' $New_Folder/raml_results_hgsoc.txt
  fi
  gene=$(awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1; exit }' $New_Folder/raml_results_hgsoc.txt)
  if [[ $gene == *\"* ]]
  then
    awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, "NA", "NA", $16, $10, "hgsoc"; exit }' $New_Folder/raml_results_hgsoc.txt
  elif [[ $i -ne 1 ]]
  then
    deldup=$(awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $16; exit }' $New_Folder/raml_results_hgsoc.txt)
    casenum=$( $New_Folder/selectNumPatientsAndCNVs.sh $gene $deldup "hgsoc" )
    ctrlnum=$( $New_Folder/selectNumPatientsAndCNVs.sh $gene $deldup "control" )
    awk -v row=$i -v cases=$casenum -v ctrls=$ctrlnum 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, cases, ctrls, $16, $10, "hgsoc"; exit }' $New_Folder/raml_results_hgsoc.txt
  fi
done >| $New_Folder/Tables/raml_results_hgsoc_counts.txt


wc -l $New_Folder/raml_results_hgsoc_both.txt
wc -l $New_Folder/raml_results_hgsoc_both.txt
for (( i = 1; i <= 3437; i++ )); do
  if [[ $i -eq 1 ]]; then
    awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, "NumCases", "NumControls", "cnv_call", $10, "AnalysisSet"; exit }' $New_Folder/raml_results_hgsoc_both.txt
  fi
  gene=$(awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1; exit }' $New_Folder/raml_results_hgsoc_both.txt)
  if [[ $gene == *\"* ]]
  then
    awk -v row=$i 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, "NA", "NA", "both", $10, "hgsoc"; exit }' $New_Folder/raml_results_hgsoc_both.txt
  elif [[ $i -ne 1 ]]
  then
    casenum=$( $New_Folder/selectNumPatientsAndCNVs.sh $gene "both" "hgsoc" )
    ctrlnum=$( $New_Folder/selectNumPatientsAndCNVs.sh $gene "both" "control" )
    awk -v row=$i -v cases=$casenum -v ctrls=$ctrlnum 'BEGIN{OFS="\t"} NR==row { gsub("\r", ""); print $1, $2, $3, cases, ctrls, "both", $10, "hgsoc"; exit }' $New_Folder/raml_results_hgsoc_both.txt
  fi
done >| $New_Folder/Tables/raml_results_hgsoc_counts_both.txt



######################################################
### Combining Overall and HGSOC into sorted tables ###
######################################################

awk 'BEGIN{OFS="\t"} (FNR==1||$1=="Gene"){next} {print}' $New_Folder/Tables/raml_results_overall_counts.txt $New_Folder/Tables/raml_results_hgsoc_counts.txt | sort -k 7n | awk 'BEGIN{OFS="\t"; print "Gene", "AMLstat", "factor1", "NumCases", "NumControls", "cnv_call", "AMLpvalue", "AnalysisSet"} {print}' >| $New_Folder/Tables/raml_results_overall_hgsoc_counts_sorted.txt

awk '($1=="Gene"){print}' $New_Folder/Tables/raml_results_overall_hgsoc_counts_sorted.txt | wc -l
awk '($1=="Gene"){print}' $New_Folder/Tables/raml_results_overall_counts.txt | wc -l
awk '($1=="Gene"){print}' $New_Folder/Tables/raml_results_hgsoc_counts.txt | wc -l # ISSUE FOR LATER



awk 'BEGIN{OFS="\t"} (FNR==1){next} {print}' $New_Folder/Tables/raml_results_overall_counts_both.txt $New_Folder/Tables/raml_results_hgsoc_counts_both.txt | sort -k 7n | awk 'BEGIN{OFS="\t"; print "Gene", "AMLstat", "factor1", "NumCases", "NumControls", "cnv_call", "AMLpvalue", "AnalysisSet"} {print}' >| $New_Folder/Tables/raml_results_overall_hgsoc_counts_both_sorted.txt


awk 'BEGIN{OFS="\t"} (FNR==1){next} ($7<0.002){print}' $New_Folder/Tables/raml_results_overall_hgsoc_counts_sorted.txt $New_Folder/Tables/raml_results_overall_hgsoc_counts_both_sorted.txt | sort -k 7n | awk 'BEGIN{OFS="\t"; print "Gene", "AMLstat", "factor1", "NumCases", "NumControls", "cnv_call", "AMLpvalue", "AnalysisSet"} {print}' >| $New_Folder/Tables/raml_results_overall_hgsoc_counts_del_dups_both_sorted_top.txt

awk 'BEGIN{OFS="\t"} (FNR==1){next} ($7<0.002){print}' $New_Folder/Tables/raml_results_overall_hgsoc_counts_sorted.txt | sort -k 7n | awk 'BEGIN{OFS="\t"; print "Gene", "AMLstat", "factor1", "NumCases", "NumControls", "cnv_call", "AMLpvalue", "AnalysisSet"} {print}' >| $New_Folder/Tables/raml_results_overall_hgsoc_counts_del_dups_sorted_top.txt
