#!/usr/bin/env bash
#
#
# getGWASloci_allnonmucinous.sh
#

GWAS_Folder="/home/devriesa/common/OncoArrayCNVs/CombinedOncoArrayFiles"


awk 'BEGIN{OFS="\t"} (FNR==1||$11<5e-8){print}' $GWAS_Folder/lift_all_non_mucinous_converted_hg38.txt >| $GWAS_Folder/GWAS_sig_loci_all_non_mucinous_hg38.txt

awk 'BEGIN{OFS="\t"} (FNR!=1){print $1, $2, $3, $11}' $GWAS_Folder/GWAS_sig_loci_all_non_mucinous_hg38.txt >| $GWAS_Folder/GWAS_sig_loci_all_non_mucinous_hg38.bed


bedtools merge -i $GWAS_Folder/GWAS_sig_loci_all_non_mucinous_hg38.bed -d 250000 | wc -l
# 32


awk 'BEGIN{OFS="\t"} (FNR==1||$11<5e-8){print}' $GWAS_Folder/lift_serous_hg_extra_converted_hg38.txt >| $GWAS_Folder/GWAS_sig_loci_serous_hg_extra_hg38.txt

awk 'BEGIN{OFS="\t"} (FNR!=1){print $1, $2, $3, $11}' $GWAS_Folder/GWAS_sig_loci_serous_hg_extra_hg38.txt >| $GWAS_Folder/GWAS_sig_loci_serous_hg_extra_hg38.bed


bedtools merge -i $GWAS_Folder/GWAS_sig_loci_serous_hg_extra_hg38.bed -d 250000 | wc -l
#33
