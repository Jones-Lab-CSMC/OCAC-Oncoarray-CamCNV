# Script_Get_Bonferroni_CNVs.sh

New_Folder="/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder="/Users/devriesa/Documents/OncoArrayCNVs"

# Overall Significant SNPs

awk 'BEGIN{OFS="\t"} (FNR>1&&FNR<9){print $1, $2-1, $2, $6}' $New_Folder/Tables/ocac_overall_minp_regions_top_hg38.txt >| $New_Folder/bonferroni_sig_SNPs_overall.bed

bedtools intersect -loj -wb -wa -a $New_Folder/bonferroni_sig_SNPs_overall.bed -b $New_Folder/ocac_overall_cnvcalls_hg38.txt | awk '($4==$9){print}' | wc -l # 145

bedtools intersect -loj -wb -wa -a $New_Folder/bonferroni_sig_SNPs_overall.bed -b $New_Folder/ocac_overall_cnvcalls_hg38.txt | awk '($4==$9){print $5, $6, $7, $9}' | sort | uniq | wc -l # 16

# HGSOC Significant SNPs

awk 'BEGIN{OFS="\t"} (FNR>1&&FNR<8){print $1, $2-1, $2, $6}' $New_Folder/Tables/ocac_hgsoc_minp_regions_top_hg38.txt >| $New_Folder/bonferroni_sig_SNPs_hgsoc.bed

bedtools intersect -loj -wb -wa -a $New_Folder/bonferroni_sig_SNPs_hgsoc.bed -b $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt | awk '($4==$9){print}' | wc -l # 117

bedtools intersect -loj -wb -wa -a $New_Folder/bonferroni_sig_SNPs_hgsoc.bed -b $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt | awk '($4==$9){print $5, $6, $7, $9}' | sort | uniq | wc -l # 16
