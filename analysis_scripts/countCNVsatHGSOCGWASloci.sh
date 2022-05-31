#!/usr/bin/env bash
#
#
# countCNVsatHGSOCGWASloci.sh

New_Folder="/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder="/Users/devriesa/Documents/OncoArrayCNVs"


awk 'BEGIN{OFS="\t"} (FNR!=1){print "chr"$1, $2, $3, $1":"$2"-"$3}' $New_Folder/CNVs_at_GWAS_Loci_only_LD_SNP_coords_sig05_NewQC_hgsoc.txt >| $New_Folder/CNVs_at_GWAS_Loci_only_loc_hgsoc.bed


awk 'BEGIN{FS="\t"; OFS="\t"} {print $18, $19, $20, $21, $22, $16, $17, $25, $26, $27}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_cnvs_forcount.txt

# ";"$6";"$7 removed for non-unique
awk 'BEGIN{OFS="\t"} (FNR!=1){print $1, $2, $3, $4";"$5";"$8";"$9";"$10}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_cnvs_forcount.txt >| $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_cnvs_forcount.bed

#  "pvalue", "cnv_count", removed because non-unique
bedtools intersect -wa -wb -a $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_cnvs_forcount.bed -b $New_Folder/CNVs_at_GWAS_Loci_only_loc_hgsoc.bed | awk 'BEGIN{OFS="\t"; print "Chr", "startpos", "endpos", "Chr0100000", "del_dup", "total_count", "ctrl_count", "hgsoc_case_count", "ChrLoci", "StartLoci", "EndLoci", "Chr:Start-EndLoci"} !seen[$0]++{gsub(/;/, "\t"); print}' >| $New_Folder/CNVs_at_hgsoc_GWAS_Loci_sigCNVs_info.txt


for locus in $(awk '{print $4}' $New_Folder/CNVs_at_GWAS_Loci_only_loc_hgsoc.bed); do
  echo $locus
  echo "sig total hgsoc count"
  awk -v locusid=$locus 'BEGIN{FS="\t"; OFS="\t"; total=0} ($12==locusid){total=total+$6} END{print total}' $New_Folder/CNVs_at_hgsoc_GWAS_Loci_sigCNVs_info.txt
  echo "sig hgsoc case count"
  awk -v locusid=$locus 'BEGIN{FS="\t"; OFS="\t"; total=0} ($12==locusid){total=total+$8} END{print total}' $New_Folder/CNVs_at_hgsoc_GWAS_Loci_sigCNVs_info.txt
  echo "sig hgsoc control count"
  awk -v locusid=$locus 'BEGIN{FS="\t"; OFS="\t"; total=0} ($12==locusid){total=total+$7} END{print total}' $New_Folder/CNVs_at_hgsoc_GWAS_Loci_sigCNVs_info.txt
done


### split by dels and dups

for locus in $(awk '{print $4}' $New_Folder/CNVs_at_GWAS_Loci_only_loc_hgsoc.bed); do
  echo $locus
  echo "sig del count"
  awk -v locusid=$locus 'BEGIN{FS="\t"; OFS="\t"; total=0} ($12==locusid && $5=="del"){total=total+$6} END{print total}' $New_Folder/CNVs_at_hgsoc_GWAS_Loci_sigCNVs_info.txt
  echo "sig dup count"
  awk -v locusid=$locus 'BEGIN{FS="\t"; OFS="\t"; total=0} ($12==locusid && $5=="dup"){total=total+$6} END{print total}' $New_Folder/CNVs_at_hgsoc_GWAS_Loci_sigCNVs_info.txt
done

bedtools intersect -wa -wb -a $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.bed -b $New_Folder/CNVs_at_GWAS_Loci_only_loc_hgsoc.bed | awk 'BEGIN{OFS="\t"; print "Chr", "startpos", "endpos", "ID", "del_dup", "OncID", "ocac_study", "country", "case_control", "probes", "segmean", "length", "Chr0100000", "ChrLoci", "StartLoci", "EndLoci", "Chr:Start-EndLoci"} !seen[$0]++{gsub(/;/, "\t"); print}' >| $New_Folder/CNVs_at_hgsoc_GWAS_Loci_allCNVs_info.txt

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5, $9, $17}' $New_Folder/CNVs_at_hgsoc_GWAS_Loci_allCNVs_info.txt | sort | uniq >| $New_Folder/CNVs_at_hgsoc_GWAS_Loci_allCNVs_info_uniqforcount.txt

for locus in $(awk '{print $4}' $New_Folder/CNVs_at_GWAS_Loci_only_loc_hgsoc.bed); do
  echo $locus
  echo "total hgsoc count"
  awk -v locusid=$locus 'BEGIN{FS="\t"; OFS="\t"} ($6==locusid){print}' $New_Folder/CNVs_at_hgsoc_GWAS_Loci_allCNVs_info_uniqforcount.txt | wc -l
  echo "hgsoc case/control count"
  awk -v locusid=$locus '($6==locusid){print $5}' $New_Folder/CNVs_at_hgsoc_GWAS_Loci_allCNVs_info_uniqforcount.txt | sort | uniq -c
done
