#!/usr/bin/env bash
#
#
# makeResultsTablesHg38.sh
#

New_Folder="/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder="/Users/devriesa/Documents/OncoArrayCNVs"

#####################################
### Merging CNV Regions (Overall) ###
#####################################
echo "Merging CNV Regions (Overall)"


awk 'BEGIN{OFS="\t"} ($5=="del"){print $1, $2, $3}' $New_Folder/ocac_overall_cnvcalls_hg38.txt | sort -k1,1 -k2,2n >| $New_Folder/ocac_overall_cnvcalls_dels_hg38.sorted.bed

bedtools merge -d 1 -i $New_Folder/ocac_overall_cnvcalls_dels_hg38.sorted.bed -c 1,2,3 -o count,collapse,collapse > $New_Folder/ocac_overall_cnvcalls_dels_merged_hg38.bed


awk 'BEGIN{OFS="\t"} ($5=="dup"){print $1, $2, $3}' $New_Folder/ocac_overall_cnvcalls_hg38.txt | sort -k1,1 -k2,2n >| $New_Folder/ocac_overall_cnvcalls_dups_hg38.sorted.bed

bedtools merge -d 1 -i $New_Folder/ocac_overall_cnvcalls_dups_hg38.sorted.bed -c 1,2,3 -o count,collapse,collapse > $New_Folder/ocac_overall_cnvcalls_dups_merged_hg38.bed



###################################
### Merging CNV Regions (HGSOC) ###
###################################
echo "Merging CNV Regions (HGSOC)"


awk 'BEGIN{OFS="\t"} ($5=="del"){print $1, $2, $3}' $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt | sort -k1,1 -k2,2n >| $New_Folder/ocac_hgsoc_cnvcalls_dels_hg38.sorted.bed

bedtools merge -d 1 -i $New_Folder/ocac_hgsoc_cnvcalls_dels_hg38.sorted.bed -c 1,2,3 -o count,collapse,collapse > $New_Folder/ocac_hgsoc_cnvcalls_dels_merged_hg38.bed


awk 'BEGIN{OFS="\t"} ($5=="dup"){print $1, $2, $3}' $New_Folder/ocac_hgsoc_cnvcalls_hg38.txt | sort -k1,1 -k2,2n >| $New_Folder/ocac_hgsoc_cnvcalls_dups_hg38.sorted.bed

bedtools merge -d 1 -i $New_Folder/ocac_hgsoc_cnvcalls_dups_hg38.sorted.bed -c 1,2,3 -o count,collapse,collapse > $New_Folder/ocac_hgsoc_cnvcalls_dups_merged_hg38.bed



###################################################
### Intersect Merged CNVs with Probes (Overall) ###
###################################################
echo "Intersect Merged CNVs with Probes (Overall)"


awk 'BEGIN{OFS="\t"} ($5=="del"){gsub("\r", ""); print $1, $2, $3, $5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16";"$17}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_overall_cnv_results_v7_del_probes_hg38.bed
awk 'BEGIN{OFS="\t"} ($5=="dup"){gsub("\r", ""); print $1, $2, $3, $5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16";"$17}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_overall_cnv_results_v7_dup_probes_hg38.bed
awk 'BEGIN{OFS="\t"} (FNR==1){gsub("\r", ""); print $1, $2, $3, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.header


bedtools intersect -loj -wb -wa -a $New_Folder/ocac_overall_cnvcalls_dels_merged_hg38.bed -b $New_Folder/ocac_overall_cnv_results_v7_del_probes_hg38.bed >| $New_Folder/ocac_overall_dels_probes_CNVs_hg38_notab.txt
awk 'BEGIN{OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_overall_dels_probes_CNVs_hg38_notab.txt > $New_Folder/ocac_overall_dels_probes_CNVs_hg38_nohead.txt
awk 'BEGIN{OFS="\t"; print "chr_window", "start_maxwindow", "end_maxwindow", "num_CNVs_window", "all_startpos", "all_endpos"} FNR==1{print $0; exit}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.header | gsed '$!{:a;N;s/\n/\t/;ta}' | awk 'BEGIN{OFS="\t"} {print}' > $New_Folder/ocac_overall_dels_probes_CNVs_hg38.header
cat $New_Folder/ocac_overall_dels_probes_CNVs_hg38.header $New_Folder/ocac_overall_dels_probes_CNVs_hg38_nohead.txt > $New_Folder/ocac_overall_dels_probes_CNVs_hg38.txt

bedtools intersect -loj -wb -wa -a $New_Folder/ocac_overall_cnvcalls_dups_merged_hg38.bed -b $New_Folder/ocac_overall_cnv_results_v7_dup_probes_hg38.bed >| $New_Folder/ocac_overall_dups_probes_CNVs_hg38_notab.txt
awk 'BEGIN{OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_overall_dups_probes_CNVs_hg38_notab.txt > $New_Folder/ocac_overall_dups_probes_CNVs_hg38_nohead.txt
awk 'BEGIN{OFS="\t"; print "chr_window", "start_maxwindow", "end_maxwindow", "num_CNVs_window", "all_startpos", "all_endpos"} FNR==1{print $0; exit}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.header | gsed '$!{:a;N;s/\n/\t/;ta}' | awk 'BEGIN{OFS="\t"} {print}' > $New_Folder/ocac_overall_dups_probes_CNVs_hg38.header
cat $New_Folder/ocac_overall_dups_probes_CNVs_hg38.header $New_Folder/ocac_overall_dups_probes_CNVs_hg38_nohead.txt > $New_Folder/ocac_overall_dups_probes_CNVs_hg38.txt



#################################################
### Intersect Merged CNVs with Probes (HGSOC) ###
#################################################
# Intersect Merged CNVs with Probes (HGSOC)

awk 'BEGIN{OFS="\t"} ($5=="del"){gsub("\r", ""); print $1, $2, $3, $5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16";"$17}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_del_probes_hg38.bed
awk 'BEGIN{OFS="\t"} ($5=="dup"){gsub("\r", ""); print $1, $2, $3, $5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16";"$17}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_dup_probes_hg38.bed
awk 'BEGIN{OFS="\t"} (FNR==1){gsub("\r", ""); print $1, $2, $3, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.header


bedtools intersect -loj -wb -wa -a $New_Folder/ocac_hgsoc_cnvcalls_dels_merged_hg38.bed -b $New_Folder/ocac_hgsoc_cnv_results_v7_del_probes_hg38.bed >| $New_Folder/ocac_hgsoc_dels_probes_CNVs_hg38_notab.txt
awk 'BEGIN{OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_hgsoc_dels_probes_CNVs_hg38_notab.txt > $New_Folder/ocac_hgsoc_dels_probes_CNVs_hg38_nohead.txt
awk 'BEGIN{OFS="\t"; print "chr_window", "start_maxwindow", "end_maxwindow", "num_CNVs_window", "all_startpos", "all_endpos"} FNR==1{print $0; exit}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.header | gsed '$!{:a;N;s/\n/\t/;ta}' | awk 'BEGIN{OFS="\t"} {print}' > $New_Folder/ocac_hgsoc_dels_probes_CNVs_hg38.header
cat $New_Folder/ocac_hgsoc_dels_probes_CNVs_hg38.header $New_Folder/ocac_hgsoc_dels_probes_CNVs_hg38_nohead.txt > $New_Folder/ocac_hgsoc_dels_probes_CNVs_hg38.txt

bedtools intersect -loj -wb -wa -a $New_Folder/ocac_hgsoc_cnvcalls_dups_merged_hg38.bed -b $New_Folder/ocac_hgsoc_cnv_results_v7_dup_probes_hg38.bed >| $New_Folder/ocac_hgsoc_dups_probes_CNVs_hg38_notab.txt
awk 'BEGIN{OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_hgsoc_dups_probes_CNVs_hg38_notab.txt > $New_Folder/ocac_hgsoc_dups_probes_CNVs_hg38_nohead.txt
awk 'BEGIN{OFS="\t"; print "chr_window", "start_maxwindow", "end_maxwindow", "num_CNVs_window", "all_startpos", "all_endpos"} FNR==1{print $0; exit}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.header | gsed '$!{:a;N;s/\n/\t/;ta}' | awk 'BEGIN{OFS="\t"} {print}' > $New_Folder/ocac_hgsoc_dups_probes_CNVs_hg38.header
cat $New_Folder/ocac_hgsoc_dups_probes_CNVs_hg38.header $New_Folder/ocac_hgsoc_dups_probes_CNVs_hg38_nohead.txt > $New_Folder/ocac_hgsoc_dups_probes_CNVs_hg38.txt



############################################
### Get Min P Value Per Window (Overall) ###
############################################
echo "Get Min P Value Per Window (Overall)"

awk 'BEGIN{OFS="\t"} (NR==FNR){ if (!($1";"$2";"$3 in min) || ($21 < min[$1";"$2";"$3])){min[$1";"$2";"$3] = $21; count[$1";"$2";"$3]++;} next;} ($21 == min[$1";"$2";"$3]){ if (FNR==1){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, "num_minp_probes";} else {print $0, count[$1";"$2";"$3];}}' $New_Folder/ocac_overall_dels_probes_CNVs_hg38.txt $New_Folder/ocac_overall_dels_probes_CNVs_hg38.txt >| $New_Folder/ocac_overall_dels_minp_probes_CNVs_hg38.txt
awk 'BEGIN{OFS="\t"} (NR==FNR){ if (!($1";"$2";"$3 in min) || ($21 < min[$1";"$2";"$3])){min[$1";"$2";"$3] = $21; count[$1";"$2";"$3]++;} next;} ($21 == min[$1";"$2";"$3]){ if (FNR==1){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, "num_minp_probes";} else {print $0, count[$1";"$2";"$3];}}' $New_Folder/ocac_overall_dups_probes_CNVs_hg38.txt $New_Folder/ocac_overall_dups_probes_CNVs_hg38.txt >| $New_Folder/ocac_overall_dups_minp_probes_CNVs_hg38.txt

awk '(FNR==1&&FNR!=NR){next} ($7!="."){print}' $New_Folder/ocac_overall_dels_minp_probes_CNVs_hg38.txt $New_Folder/ocac_overall_dups_minp_probes_CNVs_hg38.txt >| $New_Folder/ocac_overall_minp_probes_CNVs_hg38.txt



##########################################
### Get Min P Value Per Window (HGSOC) ###
##########################################
echo "Get Min P Value Per Window (HGSOC)"


awk 'BEGIN{OFS="\t"} (NR==FNR){ if (!($1";"$2";"$3 in min) || ($21 < min[$1";"$2";"$3])){min[$1";"$2";"$3] = $21; count[$1";"$2";"$3]++;} next;} ($21 == min[$1";"$2";"$3]){ if (FNR==1){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, "num_minp_probes";} else {print $0, count[$1";"$2";"$3];}}' $New_Folder/ocac_hgsoc_dels_probes_CNVs_hg38.txt $New_Folder/ocac_hgsoc_dels_probes_CNVs_hg38.txt >| $New_Folder/ocac_hgsoc_dels_minp_probes_CNVs_hg38.txt
awk 'BEGIN{OFS="\t"} (NR==FNR){ if (!($1";"$2";"$3 in min) || ($21 < min[$1";"$2";"$3])){min[$1";"$2";"$3] = $21; count[$1";"$2";"$3]++;} next;} ($21 == min[$1";"$2";"$3]){ if (FNR==1){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, "num_minp_probes";} else {print $0, count[$1";"$2";"$3];}}' $New_Folder/ocac_hgsoc_dups_probes_CNVs_hg38.txt $New_Folder/ocac_hgsoc_dups_probes_CNVs_hg38.txt >| $New_Folder/ocac_hgsoc_dups_minp_probes_CNVs_hg38.txt


awk '(FNR==1&&FNR!=NR){next} ($7!="."){print}' $New_Folder/ocac_hgsoc_dels_minp_probes_CNVs_hg38.txt $New_Folder/ocac_hgsoc_dups_minp_probes_CNVs_hg38.txt >| $New_Folder/ocac_hgsoc_minp_probes_CNVs_hg38.txt



###################################
### Make Final Tables (Overall) ###
###################################
echo "Make Final Tables (Overall)"


awk 'BEGIN{OFS = "\t"} (FNR==1){gsub("\r", ""); print "Chromosome", "SNP position", "CNV region start", "CNV region end", "SNP p-value", "CNV type", "SNP variant type", "Genes at SNP", "SNP EAF controls", "SNP EAF cases", "SNP log(OR)", "SNP SE", "SNP Wald", "SNP LRT", "Num CNVs at SNP"} (FNR!=1){gsub("\r", ""); print $1, $9, $2, $3, $21, $10, $11, $12, $13, $14, $17, $18, $19, $20, $22}' $New_Folder/ocac_overall_minp_probes_CNVs_hg38.txt >| $New_Folder/Tables/ocac_overall_minp_regions_probes_CNVs_hg38.txt

awk 'BEGIN{OFS="\t"} (FNR==1){gsub("\r", ""); print "Chr", "SNP position", "CNV region start", "CNV region end", "SNP p-value", "CNV type", "SNP variant type", "Genes at SNP", "SNP EAF controls", "SNP EAF cases", "SNP log(OR)"} ($5<0.0005){gsub("\r", ""); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' $New_Folder/Tables/ocac_overall_minp_regions_probes_CNVs_hg38.txt >| $New_Folder/Tables/ocac_overall_minp_regions_top_hg38.txt


#################################
### Make Final Tables (HGSOC) ###
#################################
echo "Make Final Tables (HGSOC)"

awk 'BEGIN{OFS = "\t"} (FNR==1){gsub("\r", ""); print "Chr", "SNP position", "CNV region start", "CNV region end", "SNP p-value", "CNV type", "SNP variant type", "Genes at SNP", "SNP EAF controls", "SNP EAF cases", "SNP log(OR)", "SNP SE", "SNP Wald", "SNP LRT", "Num CNVs at SNP"} (FNR!=1){gsub("\r", ""); print $1, $9, $2, $3, $21, $10, $11, $12, $13, $14, $17, $18, $19, $20, $22}' $New_Folder/ocac_hgsoc_minp_probes_CNVs_hg38.txt >| $New_Folder/Tables/ocac_hgsoc_minp_regions_probes_CNVs_hg38.txt

awk 'BEGIN{OFS="\t"} (FNR==1){gsub("\r", ""); print "Chr", "SNP position", "CNV region start", "CNV region end", "SNP p-value", "CNV type", "SNP variant type", "Genes at SNP", "SNP EAF controls", "SNP EAF cases", "SNP log(OR)"} ($5<0.0005){gsub("\r", ""); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' $New_Folder/Tables/ocac_hgsoc_minp_regions_probes_CNVs_hg38.txt >| $New_Folder/Tables/ocac_hgsoc_minp_regions_top_hg38.txt



###########################
### Make Segments Table ###
###########################


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print $7, $8, $13, $12, $11, $19, $16}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$2=="237121"){print}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt
awk 'BEGIN{OFS="\t"; FS="\t"; print "chr", "start", "end", "del_dup", "histotype", "OncID", "ocac_study", "country", "probes", "segmean", "length"} (FNR==NR && $7==0){histotype[$2]="Control"; next} (FNR==NR && $8==1){histotype[$2]="HGSOC"; next} (FNR==NR && $13==1){histotype[$2]="LGSOC"; next} (FNR==NR && $12==1){histotype[$2]="Clear_Cell"; next} (FNR==NR && $11==1){histotype[$2]="Endometrioid"; next} (FNR==NR && $19==1){histotype[$2]="Mucinous"; next} (FNR==NR && $16==1){histotype[$2]="Other_Case"; next} (FNR!=NR && FNR!=1){print $1, $2, $3, $5, histotype[$6], $6, $7, $8, $10, $11, $3-$2}' $New_Folder/oncid_phenotypes_selected_samples_20200910.txt $New_Folder/ocac_overall_cnvcalls_hg38.txt >| $New_Folder/Tables/ocac_overall_cnvcalls_with_histotypes_hg38.txt

awk 'BEGIN{FS="\t"}($6=="237121"){print}' $New_Folder/Tables/ocac_overall_cnvcalls_with_histotypes_hg38.txt



###########
### END ###
###########
