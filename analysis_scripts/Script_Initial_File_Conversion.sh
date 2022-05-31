# Script_Initial_File_Conversion.sh

New_Folder="/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder="/Users/devriesa/Documents/OncoArrayCNVs"

########################################
### Liftover CNV Calls (Overall Set) ###
########################################
echo "Liftover CNV Calls (Overall Set)"


awk 'BEGIN{FS = ","; OFS = "\t"} (FNR!=1 && FNR==NR){sub(/23/, "X", $6); print "chr"$6, $7-1, $8, $1, "del", $2, $3, $4, $5, $9, $10, $11} (FNR!=1 && FNR!=NR){sub(/23/, "X", $6); print "chr"$6, $7-1, $8, $1, "dup", $2, $3, $4, $5, $9, $10, $11}' $New_Folder/ocac_overall_export_v7_dels.csv $New_Folder/ocac_overall_export_v7_dups.csv >| $New_Folder/ocac_overall_cnvcalls_hg19.bedPlus
awk 'BEGIN{FS = ","; OFS = "\t"} (FNR==1){print $6, $7, $8, $1, "del_dup", $2, $3, $4, $5, $9, $10, $11}' $New_Folder/ocac_overall_export_v7_dels.csv >| $New_Folder/ocac_overall_cnvcalls_hg19.header

liftOver -bedPlus=4 -tab $New_Folder/ocac_overall_cnvcalls_hg19.bedPlus $CNVs_Folder/hg19ToHg38.over.chain $New_Folder/ocac_overall_cnvcalls_hg38.bedPlus $New_Folder/ocac_overall_cnvcalls_hg19tohg38_fail.txt

awk 'BEGIN{FS="\t"; OFS="\t"} !($1 ~ /alt/){gsub(/\r/, ""); print}' $New_Folder/ocac_overall_cnvcalls_hg19.header $New_Folder/ocac_overall_cnvcalls_hg38.bedPlus >| $New_Folder/ocac_overall_cnvcalls_hg38.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR>1){ print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13 }' $New_Folder/ocac_overall_cnvcalls_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_hg38.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {sub(/chr/, ""); print $1, $2, $3, $5}' $New_Folder/ocac_overall_cnvcalls_hg38.txt | awk 'BEGIN{FS="\t"; OFS="\t"} (!seen[$1";"$2";"$3]++){print $1, $2, $3, $4}' | awk 'BEGIN{FS="\t"; OFS="\t"} (NR==1){print $1, $2, $3, $1"0"(100000+NR-1), $4} (NR>1){print "chr"$1, $2, $3, $1"0"(100000+NR-1), $4}' >| $New_Folder/ocac_overall_cnvcalls_IDs_list_hg38.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==NR){ids[$1";"$2";"$3]=$4; next} (FNR!=NR){print $0, ids[$1";"$2";"$3]}' $New_Folder/ocac_overall_cnvcalls_IDs_list_hg38.txt $New_Folder/ocac_overall_cnvcalls_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print; exit}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_IDs_hg38.header
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR!=1){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_IDs_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$9==0){print}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_controls_hg38.txt
awk 'BEGIN{FS="\t"; OFS="\t"} ($9==0){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_controls_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$9==1){print}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_cases_hg38.txt
awk 'BEGIN{FS="\t"; OFS="\t"} ($9==1){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_cases_hg38.bed



bedtools intersect -v -loj -wa -a $New_Folder/ocac_overall_cnvcalls_controls_hg38.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed >| $New_Folder/ocac_overall_cnvcalls_controls_noncoding_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.header $New_Folder/ocac_overall_cnvcalls_controls_noncoding_hg38.bed >| $New_Folder/ocac_overall_cnvcalls_controls_noncoding_hg38.txt
bedtools intersect -v -loj -wa -a $New_Folder/ocac_overall_cnvcalls_cases_hg38.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed >| $New_Folder/ocac_overall_cnvcalls_cases_noncoding_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.header $New_Folder/ocac_overall_cnvcalls_cases_noncoding_hg38.bed >| $New_Folder/ocac_overall_cnvcalls_cases_noncoding_hg38.txt



awk 'BEGIN{FS="\t"; OFS="\t"} (FNR>1){seen[$13]++} (FNR>1 && $9==0){ctrl[$13]++; overall_case[$13]=overall_case[$13]+0} (FNR>1 && $9==1){overall_case[$13]++; ctrl[$13]=ctrl[$13]+0} END{ print "Chr0100000", "total_count", "ctrl_count", "overall_case_count"; PROCINFO["sorted_in"] = "@ind_num_asc"; for(id in seen) {print id, seen[id], ctrl[id], overall_case[id]}}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_IDs_hg38_OncIDcounts.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==NR){total[$1]=$2; ctrl[$1]=$3; overall_case[$1]=$4} (FNR!=NR && !seen[$13]++){print $1, $2, $3, $13, $5, $12, $10, total[$13], ctrl[$13], overall_case[$13]}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38_OncIDcounts.txt $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_IDs_hg38_unique.txt


awk 'BEGIN{FS="\t"; OFS="\t"} FNR==1{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38_unique.txt >| $New_Folder/ocac_overall_cnvcalls_IDs_hg38_unique.header
awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38_unique.txt >| $New_Folder/ocac_overall_cnvcalls_IDs_hg38_unique.bed




######################################
### Liftover CNV Calls (HGSOC Set) ###
######################################
echo "Liftover CNV Calls (HGSOC Set)"

awk 'BEGIN{sum=0; rows=0} (FNR!=1){sum+=$5; rows+=1} END{print sum/rows,rows}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt # average age
awk 'BEGIN{sum=0; rows=0} (FNR!=1&&$7==0){sum+=$5; rows+=1} END{print sum/rows, rows}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt # average control age
awk 'BEGIN{sum=0; rows=0} (FNR!=1&&$7==1){sum+=$5; rows+=1} END{print sum/rows, rows}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt # average age for overall
awk 'BEGIN{sum=0; rows=0} (FNR!=1&&$8==1){sum+=$5; rows+=1} END{print sum/rows, rows}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt # average age for hgsoc

awk '{print $8}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt | sort | uniq -c
awk '{print $9}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt | sort | uniq -c

awk 'BEGIN{OFS="\t"; hgsoc["OncID"]="serous_hg_extra"} (FNR==NR && ($8==0 || $8==1 || NR==1)){hgsoc[$2]=$8; next} (FNR!=NR && ($6 in hgsoc)){print $1, $2, $3, $4, $5, $6, $7, $8, hgsoc[$6], $10, $11, $12, $13}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt $New_Folder/ocac_overall_cnvcalls_hg38.txt | awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/\t$/, ""); print}' >| $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt



head -n 3 $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt | cat -t


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR>1){ print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13 }' $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13 }' $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_hg38.header


awk 'BEGIN{FS="\t"; OFS="\t"} {sub(/chr/, ""); print $1, $2, $3, $5}' $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt | awk 'BEGIN{FS="\t"; OFS="\t"} (!seen[$1";"$2";"$3]++){print $1, $2, $3, $4}' | awk 'BEGIN{FS="\t"; OFS="\t"} (NR==1){print $1, $2, $3, $1"0"(100000+NR-1), $4} (NR>1){print "chr"$1, $2, $3, $1"0"(100000+NR-1), $4}' >| $New_Folder/ocac_hgsoc_cnvcalls_IDs_list_hg38.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==NR){ids[$1";"$2";"$3]=$4; next} (FNR!=NR){print $0, ids[$1";"$2";"$3]}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_list_hg38.txt $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print; exit}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.header
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR!=1){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$9==0){print}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_controls_hg38.txt
awk 'BEGIN{FS="\t"; OFS="\t"} ($9==0){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_controls_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$9==1){print}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_cases_hg38.txt
awk 'BEGIN{FS="\t"; OFS="\t"} ($9==1){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_cases_hg38.bed




awk 'BEGIN{FS="\t"; OFS="\t"} (FNR>1){seen[$13]++} (FNR>1 && $9==0){ctrl[$13]++; hgsoc_case[$13]=hgsoc_case[$13]+0} (FNR>1 && $9==1){hgsoc_case[$13]++; ctrl[$13]=ctrl[$13]+0} END{ print "Chr0100000", "total_count", "ctrl_count", "hgsoc_case_count"; PROCINFO["sorted_in"] = "@ind_num_asc"; for(id in seen) {print id, seen[id], ctrl[id], hgsoc_case[id]}}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_OncIDcounts.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==NR){total[$1]=$2; ctrl[$1]=$3; hgsoc_case[$1]=$4} (FNR!=NR && !seen[$13]++){print $1, $2, $3, $13, $5, $12, $10, total[$13], ctrl[$13], hgsoc_case[$13]}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_OncIDcounts.txt $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_unique.txt

awk 'BEGIN{FS="\t"; OFS="\t"} FNR==1{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_unique.txt >| $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_unique.header
awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_unique.txt >| $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_unique.bed



######################################
### Liftover CNV Probe Information ###
######################################
echo "Liftover CNV Probe Information"

# from $New_Folder/ocac_overall_cnv_results_v7.xlsx I manually made:
# $New_Folder/ocac_overall_cnv_results_v7_dels_probes_all.txt
# $New_Folder/ocac_overall_cnv_results_v7_dups_probes_all.txt
# $New_Folder/ocac_overall_cnv_results_v7_dels_regions.txt
# $New_Folder/ocac_overall_cnv_results_v7_dups_regions.txt

# # from $New_Folder/ocac_overall_hgsoc_cnv_results_v7_all_probes.xlsx I manually made:
# $New_Folder/ocac_overall_dels_cnv_results_v7_all_probes.txt
# $New_Folder/ocac_overall_dups_cnv_results_v7_all_probes.txt
# $New_Folder/ocac_hgsoc_dels_cnv_results_v7_all_probes.txt
# $New_Folder/ocac_hgsoc_dups_cnv_results_v7_all_probes.txt


awk 'BEGIN{FS = "\t"; OFS = "\t"} (FNR!=1 && FNR==NR){sub(/23/, "X", $5); print "chr"$5, $6-1, $6, $1, "del", $3, $4, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16} (FNR!=1 && FNR!=NR){sub(/23/, "X", $6); print "chr"$4, $5-1, $5, $1, "dup", $2, $3, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' $New_Folder/ocac_overall_dels_cnv_results_v7_all_probes.txt $New_Folder/ocac_overall_dups_cnv_results_v7_all_probes.txt >| $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg19.bedPlus
awk 'BEGIN{FS = "\t"; OFS = "\t"} (FNR==1){sub(/23/, "X", $5); print $5, "position_b37_0base", $6, $1, "del_dup", $3, $4, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' $New_Folder/ocac_overall_dels_cnv_results_v7_all_probes.txt >| $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg19.header

liftOver -bedPlus=4 -tab $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg19.bedPlus $CNVs_Folder/hg19ToHg38.over.chain $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.bedPlus $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg19tohg38_fail.txt

awk 'BEGIN{FS="\t"; OFS="\t"} !($1 ~ /alt/){print}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg19.header $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.bedPlus >| $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt


awk 'BEGIN{FS = "\t"; OFS = "\t"} (FNR!=1 && FNR==NR && $5!=""){sub(/23/, "X", $5); print "chr"$5, $6-1, $6, $1, "del", $3, $4, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16} (FNR!=1 && FNR!=NR && $5!=""){sub(/23/, "X", $5); print "chr"$5, $6-1, $6, $1, "dup", $3, $4, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' $New_Folder/ocac_hgsoc_dels_cnv_results_v7_all_probes.txt $New_Folder/ocac_hgsoc_dups_cnv_results_v7_all_probes.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg19.bedPlus

awk 'BEGIN{FS = "\t"; OFS = "\t"} (FNR==1){sub(/23/, "X", $5); print $5, "position_b37_0base", $6, $1, "del_dup", $3, $4, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' $New_Folder/ocac_hgsoc_dels_cnv_results_v7_all_probes.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg19.header

liftOver -bedPlus=4 -tab $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg19.bedPlus $CNVs_Folder/hg19ToHg38.over.chain $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.bedPlus $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg19tohg38_fail.txt

awk 'BEGIN{FS="\t"; OFS="\t"} !($1 ~ /alt/){print}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg19.header $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.bedPlus >| $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt




########################################
### Subset to Significant CNV Probes ###
########################################
echo "Subset to Significant CNV Probes"

# Bonferroni
# head -n 3 $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt
# wc -l $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt #57432 without header
# # overall is p < 8.71x10^-7
# head -n 3 $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt
# wc -l $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt #58382 without header
# hgsoc is p < 8.56x10^-7

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR == 1 || $16 < 0.01){print}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_overall_cnv_results_v7_sig01_probes_hg38.txt

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR == 1 || $16 < 0.05){print}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.txt

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR>1){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16";"$17}' $New_Folder/ocac_overall_cnv_results_v7_sig01_probes_hg38.txt >| $New_Folder/ocac_overall_cnv_results_v7_sig01_probes_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print $0; exit}' $New_Folder/ocac_overall_cnv_results_v7_sig01_probes_hg38.txt >| $New_Folder/ocac_overall_cnv_results_v7_sig01_probes_hg38.header

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR>1){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16";"$17}' $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.txt >| $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print $0; exit}' $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.txt >| $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.header

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR == 1 || $16 < 0.01){print}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_sig01_probes_hg38.txt

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR == 1 || $16 < 0.05){print}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.txt

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR>1){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16";"$17}' $New_Folder/ocac_hgsoc_cnv_results_v7_sig01_probes_hg38.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_sig01_probes_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print $0; exit}' $New_Folder/ocac_hgsoc_cnv_results_v7_sig01_probes_hg38.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_sig01_probes_hg38.header

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR>1){print $1, $2, $3, $4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16";"$17}' $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.bed
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print $0; exit}' $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.header



#########################################################
### Intersect Sig Probes with CNV Calls (Overall Set) ###
#########################################################
echo "Intersect Sig Probes with CNV Calls (Overall Set)"

bedtools intersect -loj -wb -wa -a $New_Folder/ocac_overall_cnv_results_v7_sig01_probes_hg38.bed -b $New_Folder/ocac_overall_cnvcalls_IDs_hg38_unique.bed >| $New_Folder/ocac_overall_sig01_CNVs_hg38_notab.txt

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_overall_sig01_CNVs_hg38_notab.txt >| $New_Folder/ocac_overall_sig01_CNVs_hg38_nohead.txt


cat $New_Folder/ocac_overall_cnv_results_v7_sig01_probes_hg38.header $New_Folder/ocac_overall_cnvcalls_IDs_hg38_unique.header | gsed '$!{:a;N;s/\r//;s/\n/\t/;ta}' >| $New_Folder/ocac_overall_sig01_CNVs_hg38.header

cat $New_Folder/ocac_overall_sig01_CNVs_hg38.header $New_Folder/ocac_overall_sig01_CNVs_hg38_nohead.txt >| $New_Folder/ocac_overall_sig01_CNVs_hg38_mismatch.txt


awk '(FNR==1||$18=="."){print $0}' $New_Folder/ocac_overall_sig01_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_overall_sig01_CNVs_hg38_nomatch.txt

rm $New_Folder/ocac_overall_sig01_CNVs_hg38_notab.txt
rm $New_Folder/ocac_overall_sig01_CNVs_hg38_nohead.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$5==$22){print $0}' $New_Folder/ocac_overall_sig01_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_overall_sig01_CNVs_hg38.txt


awk 'BEGIN{FS="\t"; OFS="\t"} NR == FNR { if (!($21 in min) || $16 < min[$21]) { min[$21] = $16; count[$21]++; } next; } $16 == min[$21] { if (FNR==1) { print $0, "num_sig_probes"; } else { print $0, count[$21]; } }' $New_Folder/ocac_overall_sig01_CNVs_hg38.txt $New_Folder/ocac_overall_sig01_CNVs_hg38.txt | awk 'BEGIN{FS="\t"; OFS="\t"} !seen[$21]++' >| $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.txt

awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16}' $New_Folder/ocac_overall_sig01_CNVs_hg38.txt >| $New_Folder/ocac_overall_sig01_CNVs_hg38.bed
# name is Chr0100000;del_dup;total_count;ctrl_count;overall_case_count;position_b37;OR;pvalue

awk 'BEGIN{FS="\t"; OFS="\t"} FNR==1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16";"$28}' $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.txt >| $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.bedheader
awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16";"$28}' $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.txt >| $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.bed
# name is Chr0100000;del_dup;total_count;ctrl_count;overall_case_count;position_b37;OR;pvalue;num_sig_probes


bedtools intersect -loj -wb -wa -a $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.bed -b $New_Folder/ocac_overall_cnvcalls_IDs_hg38_unique.bed >| $New_Folder/ocac_overall_sig05_CNVs_hg38_notab.txt

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_overall_sig05_CNVs_hg38_notab.txt >| $New_Folder/ocac_overall_sig05_CNVs_hg38_nohead.txt


cat $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.header $New_Folder/ocac_overall_cnvcalls_IDs_hg38_unique.header | gsed '$!{:a;N;s/\r//;s/\n/\t/;ta}' >| $New_Folder/ocac_overall_sig05_CNVs_hg38.header


cat $New_Folder/ocac_overall_sig05_CNVs_hg38.header $New_Folder/ocac_overall_sig05_CNVs_hg38_nohead.txt >| $New_Folder/ocac_overall_sig05_CNVs_hg38_mismatch.txt

awk '(FNR==1||$18=="."){print $0}' $New_Folder/ocac_overall_sig05_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_overall_sig05_CNVs_hg38_nomatch.txt

rm $New_Folder/ocac_overall_sig05_CNVs_hg38_notab.txt
rm $New_Folder/ocac_overall_sig05_CNVs_hg38_nohead.txt


# confirming dels match dels and dups match dups with $5==$22
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$5==$22){print $0}' $New_Folder/ocac_overall_sig05_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_overall_sig05_CNVs_hg38.txt

awk 'BEGIN{FS="\t"; OFS="\t"} NR == FNR { if (!($21 in min) || $16 < min[$21]) { min[$21] = $16; count[$21]++; } next; } $16 == min[$21] { if (FNR==1) { print $0, "num_sig_probes"; } else { print $0, count[$21]; } }' $New_Folder/ocac_overall_sig05_CNVs_hg38.txt $New_Folder/ocac_overall_sig05_CNVs_hg38.txt | awk 'BEGIN{FS="\t"; OFS="\t"} !seen[$21]++' >| $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.txt



awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16}' $New_Folder/ocac_overall_sig05_CNVs_hg38.txt >| $New_Folder/ocac_overall_sig05_CNVs_hg38.bed
# name is Chr0100000;del_dup;total_count;ctrl_count;overall_case_count;position_b37;OR;pvalue


awk 'BEGIN{FS="\t"; OFS="\t"} FNR==1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16";"$28}' $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.txt >| $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.bedheader
awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16";"$28}' $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.txt >| $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.bed
# name is Chr0100000;del_dup;total_count;ctrl_count;overall_case_count;position_b37;OR;pvalue;num_sig_probes




#######################################################
### Intersect Sig Probes with CNV Calls (HGSOC Set) ###
#######################################################
echo "Intersect Sig Probes with CNV Calls (HGSOC Set)"

bedtools intersect -loj -wb -wa -a $New_Folder/ocac_hgsoc_cnv_results_v7_sig01_probes_hg38.bed -b $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_unique.bed >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_notab.txt

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_notab.txt >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_nohead.txt


cat $New_Folder/ocac_hgsoc_cnv_results_v7_sig01_probes_hg38.header $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_unique.header | gsed '$!{:a;N;s/\r//;s/\n/\t/;ta}' >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38.header

cat $New_Folder/ocac_hgsoc_sig01_CNVs_hg38.header $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_nohead.txt >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_mismatch.txt

awk '(FNR==1||$18=="."){print $0}' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_nomatch.txt

rm $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_notab.txt
rm $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_nohead.txt

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$5==$22){print $0}' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38.txt

awk 'BEGIN{FS="\t"; OFS="\t"} NR == FNR { if (!($21 in min) || $16 < min[$21]) { min[$21] = $16; count[$21]++; } next; } $16 == min[$21] { if (FNR==1) { print $0, "num_sig_probes"; } else { print $0, count[$21]; } }' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38.txt $New_Folder/ocac_hgsoc_sig01_CNVs_hg38.txt | awk 'BEGIN{FS="\t"; OFS="\t"} !seen[$21]++' >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.txt


awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16}' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38.txt >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38.bed
# name is Chr0100000;del_dup;total_count;ctrl_count;hgsoc_case_count;position_b37;OR;pvalue

awk 'BEGIN{FS="\t"; OFS="\t"} FNR==1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16";"$28}' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.txt >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.bedheader
awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16";"$28}' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.txt >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.bed
# name is Chr0100000;del_dup;total_count;ctrl_count;hgsoc_case_count;position_b37;OR;pvalue;num_sig_probes


bedtools intersect -loj -wb -wa -a $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.bed -b $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_unique.bed >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_notab.txt

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_notab.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_nohead.txt


cat $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.header $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_unique.header | gsed '$!{:a;N;s/\r//;s/\n/\t/;ta}' >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38.header


cat $New_Folder/ocac_hgsoc_sig05_CNVs_hg38.header $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_nohead.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_mismatch.txt


awk '(FNR==1||$18=="."){print $0}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_nomatch.txt


rm $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_notab.txt
rm $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_nohead.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$5==$22){print $0}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38.txt


awk 'BEGIN{FS="\t"; OFS="\t"} NR == FNR { if (!($21 in min) || $16 < min[$21]) { min[$21] = $16; count[$21]++; } next; } $16 == min[$21] { if (FNR==1) { print $0, "num_sig_probes"; } else { print $0, count[$21]; } }' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38.txt $New_Folder/ocac_hgsoc_sig05_CNVs_hg38.txt | awk 'BEGIN{FS="\t"; OFS="\t"} !seen[$21]++' >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.txt



awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38.bed
# name is Chr0100000;del_dup;total_count;ctrl_count;hgsoc_case_count;position_b37;OR;pvalue

awk 'BEGIN{FS="\t"; OFS="\t"} FNR==1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16";"$28}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.bedheader
awk 'BEGIN{FS="\t"; OFS="\t"} FNR>1{print $18, $19, $20, $21";"$22";"$25";"$26";"$27";"$3";"$12";"$16";"$28}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.bed
# name is Chr0100000;del_dup;total_count;ctrl_count;hgsoc_case_count;position_b37;OR;pvalue;num_sig_probes



##################################################
##### Get Non-Coding CNVs Only (Overall Set) #####
##################################################
echo "Get Non-Coding CNVs Only (Overall Set)"

bedtools intersect -v -loj -wa -a $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed >| $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38.bed


awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.bedheader $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38.bed >| $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38.txt


bedtools intersect -v -loj -wa -a $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed >| $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38.bed


awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.bedheader $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38.bed >| $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38.txt



####################################################
##### Non-Coding CNVs (Overall Set) Simulation #####
####################################################

echo "Non-Coding CNVs (Overall Set) Simulation"

##### Sig 01 #####

:> $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_forsim2000.bed

for i in {1..2000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38.bed >> $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_forsim2000.bed; done

bedtools shuffle -i $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_forsim2000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_sim2000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_sim2000_all.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_sim2000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.bedheader $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_sim2000.bed >| $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_sim2000.txt


rm $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_forsim2000.bed
rm $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_sim2000_all.bed
rm $New_Folder/ocac_overall_sig01_noncoding_CNVs_hg38_sim2000.bed



##### Sig 05 #####

:> $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_forsim2000.bed

for i in {1..2000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38.bed >> $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_forsim2000.bed; done

bedtools shuffle -i $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_forsim2000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_sim2000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_sim2000_all.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_sim2000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.bedheader $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_sim2000.bed >| $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_sim2000.txt


rm $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_forsim2000.bed
rm $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_sim2000_all.bed
rm $New_Folder/ocac_overall_sig05_noncoding_CNVs_hg38_sim2000.bed





################################################
##### Get Non-Coding CNVs Only (HGSOC Set) #####
################################################
echo "Get Non-Coding CNVs Only (HGSOC Set)"

bedtools intersect -v -loj -wa -a $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed >| $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38.bed


awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.bedheader $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38.bed >| $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38.txt


bedtools intersect -v -loj -wa -a $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed >| $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38.bed


awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.bedheader $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38.bed >| $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38.txt




##################################################
##### Non-Coding CNVs (HGSOC Set) Simulation #####
##################################################

echo "Non-Coding CNVs (HGSOC Set) Simulation"

##### Sig 01 #####

:> $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_forsim2000.bed

for i in {1..2000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38.bed >> $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_forsim2000.bed; done

bedtools shuffle -i $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_forsim2000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_sim2000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_sim2000_all.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_sim2000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.bedheader $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_sim2000.bed >| $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_sim2000.txt


rm $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_forsim2000.bed
rm $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_sim2000_all.bed
rm $New_Folder/ocac_hgsoc_sig01_noncoding_CNVs_hg38_sim2000.bed



##### Sig 05 #####

:> $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_forsim2000.bed

for i in {1..2000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38.bed >> $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_forsim2000.bed; done

bedtools shuffle -i $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_forsim2000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_sim2000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_sim2000_all.bed -b $CNVs_Folder/Coding_hg38/All_Coding_hg38.bed $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_sim2000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.bedheader $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_sim2000.bed >| $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_sim2000.txt


rm $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_forsim2000.bed
rm $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_sim2000_all.bed
rm $New_Folder/ocac_hgsoc_sig05_noncoding_CNVs_hg38_sim2000.bed


#############################################
##### All CNVs (Overall Set) Simulation #####
#############################################

echo "All CNVs (Overall Set) Simulation"

##### Sig 01 #####

:> $New_Folder/ocac_overall_sig01_CNVs_hg38_forsim2000.bed

for i in {1..2000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.bed >> $New_Folder/ocac_overall_sig01_CNVs_hg38_forsim2000.bed; done

bedtools shuffle -i $New_Folder/ocac_overall_sig01_CNVs_hg38_forsim2000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_overall_sig01_CNVs_hg38_sim2000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_overall_sig01_CNVs_hg38_sim2000_all.bed -b $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_overall_sig01_CNVs_hg38_sim2000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_overall_sig01_CNVs_hg38_unique.bedheader $New_Folder/ocac_overall_sig01_CNVs_hg38_sim2000.bed >| $New_Folder/ocac_overall_sig01_CNVs_hg38_sim2000.txt


rm $New_Folder/ocac_overall_sig01_CNVs_hg38_forsim2000.bed
rm $New_Folder/ocac_overall_sig01_CNVs_hg38_sim2000_all.bed
rm $New_Folder/ocac_overall_sig01_CNVs_hg38_sim2000.bed


##### Sig 05 #####

:> $New_Folder/ocac_overall_sig05_CNVs_hg38_forsim2000.bed

for i in {1..2000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.bed >> $New_Folder/ocac_overall_sig05_CNVs_hg38_forsim2000.bed; done

bedtools shuffle -i $New_Folder/ocac_overall_sig05_CNVs_hg38_forsim2000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_overall_sig05_CNVs_hg38_sim2000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_overall_sig05_CNVs_hg38_sim2000_all.bed -b $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_overall_sig05_CNVs_hg38_sim2000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.bedheader $New_Folder/ocac_overall_sig05_CNVs_hg38_sim2000.bed >| $New_Folder/ocac_overall_sig05_CNVs_hg38_sim2000.txt


rm $New_Folder/ocac_overall_sig05_CNVs_hg38_forsim2000.bed
rm $New_Folder/ocac_overall_sig05_CNVs_hg38_sim2000_all.bed
rm $New_Folder/ocac_overall_sig05_CNVs_hg38_sim2000.bed


##### Sig 05 (1000x) #####

:> $New_Folder/ocac_overall_sig05_CNVs_hg38_forsim1000.bed

for i in {1..1000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.bed >> $New_Folder/ocac_overall_sig05_CNVs_hg38_forsim1000.bed; done

bedtools shuffle -i $New_Folder/ocac_overall_sig05_CNVs_hg38_forsim1000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_overall_sig05_CNVs_hg38_sim1000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_overall_sig05_CNVs_hg38_sim1000_all.bed -b $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_overall_sig05_CNVs_hg38_sim1000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.bedheader $New_Folder/ocac_overall_sig05_CNVs_hg38_sim1000.bed >| $New_Folder/ocac_overall_sig05_CNVs_hg38_sim1000_new.txt


rm $New_Folder/ocac_overall_sig05_CNVs_hg38_forsim1000.bed
rm $New_Folder/ocac_overall_sig05_CNVs_hg38_sim1000_all.bed
rm $New_Folder/ocac_overall_sig05_CNVs_hg38_sim1000.bed



#############################################
##### All CNVs (HGSOC Set) Simulation #####
#############################################

echo "All CNVs (HGSOC Set) Simulation"

##### Sig 01 #####

:> $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_forsim2000.bed

for i in {1..2000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.bed >> $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_forsim2000.bed; done

bedtools shuffle -i $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_forsim2000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_sim2000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_sim2000_all.bed -b $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_sim2000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_unique.bedheader $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_sim2000.bed >| $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_sim2000.txt


rm $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_forsim2000.bed
rm $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_sim2000_all.bed
rm $New_Folder/ocac_hgsoc_sig01_CNVs_hg38_sim2000.bed


##### Sig 05 #####

:> $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_forsim2000.bed

for i in {1..2000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.bed >> $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_forsim2000.bed; done

bedtools shuffle -i $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_forsim2000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim2000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim2000_all.bed -b $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim2000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.bedheader $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim2000.bed >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim2000.txt


rm $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_forsim2000.bed
rm $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim2000_all.bed
rm $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim2000.bed

##### Sig 05 (1000x) #####

:> $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_forsim1000.bed

for i in {1..1000};do awk 'BEGIN{FS="\t"; OFS="\t"} { print $1, $2, $3, $4 }' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.bed >> $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_forsim1000.bed; done

bedtools shuffle -i $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_forsim1000.bed -g $CNVs_Folder/hg38.chrom.sizes.limited.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim1000_all.bed

bedtools intersect -v -loj -wa -a $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim1000_all.bed -b $CNVs_Folder/gap.bed $CNVs_Folder/centromeres.bed >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim1000.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t")} {print} ' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.bedheader $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim1000.bed >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim1000_new.txt

rm $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_forsim1000.bed
rm $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim1000_all.bed
rm $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_sim1000.bed

### END
