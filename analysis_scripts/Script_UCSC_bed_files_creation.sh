# Script_UCSC_bed_files_creation.sh

New_Folder="/Users/devriesa/Documents/OncoArrayCNVs/New_QC"
CNVs_Folder="/Users/devriesa/Documents/OncoArrayCNVs"

# browser position chr17:43,030,000-43,180,000\n

############################
### Probes (Overall Set) ###
############################

# $12 is OR, $16 is pvalue, $5 is cnv_type, $10 is NumControls, $11 is NumCases
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"AllProbes_OverallSet\" description=\"OCAC Probes from All Cases and Controls Used for Single Probe CNV Testing from Joe Dennis, name is NumControls;NumCases;OR;pvalue;cnv_type, color is red for dels and blue for dups (colors are dark for risk, light for protective, medium colored grey for 0.01<p<0.05, barely colored grey for non-significant (p>=0.05))\" itemRgb=\"on\""; next} ($12>0 && $16<0.01 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $16<0.01 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $16<0.01 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $16<0.01 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "153,153,255,"; next} ($16>=0.01 && $16<0.05 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "153,64,64,"; next} ($16>=0.01 && $16<0.05 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "64,64,153,"; next} ($16>=0.05 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "96,64,64,"; next} ($16>=0.05 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "64,64,96,"; next}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_overall_all_probes_hg38_UCSC.bed

#only significant probes
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"SigProbes_OverallSet\" description=\"Significant (p < 0.05) OCAC Probes from All Cases and Controls Used for Single Probe CNV Testing from Joe Dennis, name is NumControls;NumCases;OR;pvalue;cnv_type, color is red for dels and blue for dups (colors are dark for risk, light for protective)\" itemRgb=\"on\""; next} ($12>0 && $16<0.05 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $16<0.05 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $16<0.05 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $16<0.05 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "153,153,255,"; next} ($16>=0.05){next}' $New_Folder/ocac_overall_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_overall_sig_probes_hg38_UCSC.bed

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"Alt_SigProbes_OverallSet\" description=\"Significant (p < 0.05) OCAC Probes from All Cases and Controls Used for Single Probe CNV Testing from Joe Dennis, name is NumControls;NumCases;OR;pvalue;cnv_type, color is red for dels and blue for dups (colors are dark for risk, light for protective)\" itemRgb=\"on\""; next} ($12>0 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "153,153,255,"; next} ($16>=0.05){next}' $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_overall_sig_probes_hg38_UCSC_alt.bed




############################
### Probes (HGSOC Set) ###
############################


# $12 is OR, $16 is pvalue, $5 is cnv_type, $10 is NumControls, $11 is NumCases
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"AllProbes_HGSOCSet\" description=\"OCAC Probes from HGSOC Cases and Controls Used for Single Probe CNV Testing from Joe Dennis, name is NumControls;NumCases;OR;pvalue;cnv_type, color is red for dels and blue for dups (colors are dark for risk, light for protective, medium colored grey for 0.01<p<0.05, barely colored grey for non-significant (p>=0.05))\" itemRgb=\"on\""; next} ($12>0 && $16<0.01 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $16<0.01 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $16<0.01 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $16<0.01 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "153,153,255,"; next} ($16>=0.01 && $16<0.05 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "153,64,64,"; next} ($16>=0.01 && $16<0.05 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "64,64,153,"; next} ($16>=0.05 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "96,64,64,"; next} ($16>=0.05 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "64,64,96,"; next}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_all_probes_hg38_UCSC.bed

#only significant probes
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"SigProbes_HGSOCSet\" description=\"Significant (p < 0.05) OCAC Probes from HGSOC Cases and Controls Used for Single Probe CNV Testing from Joe Dennis, name is NumControls;NumCases;OR;pvalue;cnv_type, color is red for dels and blue for dups (colors are dark for risk, light for protective)\" itemRgb=\"on\""; next} ($12>0 && $16<0.05 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $16<0.05 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $16<0.05 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $16<0.05 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "153,153,255,"; next} ($16>=0.05){next}' $New_Folder/ocac_hgsoc_cnv_results_v7_all_probes_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_sig_probes_hg38_UCSC.bed

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"Alt_SigProbes_HGSOCSet\" description=\"Significant (p < 0.05) OCAC Probes from HGSOC Cases and Controls Used for Single Probe CNV Testing from Joe Dennis, name is NumControls;NumCases;OR;pvalue;cnv_type, color is red for dels and blue for dups (colors are dark for risk, light for protective)\" itemRgb=\"on\""; next} ($12>0 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $5=="del"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $5=="dup"){print $1, $2, $3, $10";"$11";"$12";"$16";"$5, 0, ".", $2, $3, "153,153,255,"; next} ($16>=0.05){next}' $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_sig_probes_hg38_UCSC_alt.bed





###################################
### Sig CNV Calls (Overall Set) ###
###################################


# $12 is OR, $16 is pvalue, $5 is cnv_type
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"SigCNVs_OverallSet\" description=\"Significant CNVs (p<0.05 probes intersecting) from All Cases and Controls Used for Single Probe CNV Testing by Country from Joe Dennis, name is del_dup;total_count;ctrl_count;oc_case_count;SNP_pos;OR;pvalue, color is red for dels and blue for dups (colors are dark for risk, light for protective)\" itemRgb=\"on\""; next} ($12>0 && $5=="del"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $5=="del"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $5=="dup"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $5=="dup"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16, 0, ".", $2, $3, "153,153,255,"; next}' $New_Folder/ocac_overall_sig05_CNVs_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_overall_sig_CNVs_hg38_UCSC.bed


# $12 is OR, $16 is pvalue, $5 is cnv_type
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"SigCNVs_Unique_OverallSet\" description=\"Significant Unique CNVs (p<0.05 probes intersecting) from All Cases and Controls Used for Single Probe CNV Testing by Country from Joe Dennis, name is del_dup;total_count;ctrl_count;oc_case_count;minSNP_pos;OR;minpvalue;num_sig_probes, color is red for dels and blue for dups (colors are dark for risk, light for protective)\" itemRgb=\"on\""; next} ($12>0 && $5=="del"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16";"$27, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $5=="del"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16";"$27, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $5=="dup"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16";"$27, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $5=="dup"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16";"$27, 0, ".", $2, $3, "153,153,255,"; next}' $New_Folder/ocac_overall_sig05_CNVs_hg38_unique.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_overall_sig_CNVs_hg38_unique_UCSC.bed





#################################
### Sig CNV Calls (HGSOC Set) ###
#################################


# $12 is OR, $16 is pvalue, $5 is cnv_type
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"SigCNVs_HGSOCSet\" description=\"Significant CNVs (p<0.05 probes intersecting) from HGSOC Cases and Controls Used for Single Probe CNV Testing by Country from Joe Dennis, name is del_dup;total_count;ctrl_count;oc_case_count;SNP_pos;OR;pvalue, color is red for dels and blue for dups (colors are dark for risk, light for protective)\" itemRgb=\"on\""; next} ($12>0 && $5=="del"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $5=="del"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $5=="dup"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $5=="dup"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16, 0, ".", $2, $3, "153,153,255,"; next}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_sig_CNVs_hg38_UCSC.bed


# $12 is OR, $16 is pvalue, $5 is cnv_type
awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"SigCNVs_Unique_HGSOCSet\" description=\"Significant Unique CNVs (p<0.05 probes intersecting) from HGSOC Cases and Controls Used for Single Probe CNV Testing by Country from Joe Dennis, name is del_dup;total_count;ctrl_count;oc_case_count;minSNP_pos;OR;minpvalue;num_sig_probes, color is red for dels and blue for dups (colors are dark for risk, light for protective)\" itemRgb=\"on\""; next} ($12>0 && $5=="del"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16";"$27, 0, ".", $2, $3, "255,0,0,"; next} ($12<0 && $5=="del"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16";"$27, 0, ".", $2, $3, "255,153,153,"; next} ($12>0 && $5=="dup"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16";"$27, 0, ".", $2, $3, "0,0,255,"; next} ($12<0 && $5=="dup"){print $18, $19, $20, $22";"$24";"$25";"$26";"$3";"$12";"$16";"$27, 0, ".", $2, $3, "153,153,255,"; next}' $New_Folder/ocac_hgsoc_sig05_CNVs_hg38_unique.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_sig_CNVs_hg38_unique_UCSC.bed




########################
### Combine All UCSC ###
########################

cat $New_Folder/UCSC_hg38_BED_Files/ocac_overall_all_probes_hg38_UCSC.bed $New_Folder/UCSC_hg38_BED_Files/ocac_overall_sig_probes_hg38_UCSC.bed $New_Folder/UCSC_hg38_BED_Files/ocac_overall_sig_CNVs_hg38_UCSC.bed $New_Folder/UCSC_hg38_BED_Files/ocac_overall_sig_CNVs_hg38_unique_UCSC.bed $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_all_probes_hg38_UCSC.bed $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_sig_probes_hg38_UCSC.bed $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_sig_CNVs_hg38_UCSC.bed $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_sig_CNVs_hg38_unique_UCSC.bed >| $New_Folder/UCSC_hg38_BED_Files/ALL_UCSC_tracks.txt


############################################
### Make Colored ChromHMM Consensus BEDs ###
############################################


# 1_Weak_Promoter    255,105,105
# 2_Active_Promoter    255,0,0
# 3_Active_Region    231,84,128
# 4_Active_Enhancer    250,202,0
# 5_Weak_Enhancer    255,252,4
# 6_Insulator    10,190,254
# 7_Transcribed    0,176,80
# 8_Low_signal    245,245,245

### CCOC ###

echo "track name=\"CCOC_consensus_ChromHMM\" description=\"ChromHMM Output for Clear Cell Ovarian Cancer Lines Consensus - Colors: State 1 (Light Red) - Weak Promoter, State 2 (Bright Red) - Active Promoter, State 3 (Purple Pink) - Active Region, State 4 (Orange Yellow) - Active Enhancer, State 5 (Bright Yellow) - Weak Enhancer, State 6 (Light Blue) - Insulator, State 7 (Green) - Transcribed\" visibility=1 itemRgb=\"on\"" >| $New_Folder/UCSC_hg38_BED_Files/CCOC_consensus_ChromHMM.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "1_Weak_Promoter", "0", ".", $2, $3, "255,105,105"}' $New_Folder/consensus_chromHMM_states/state_E1/CCOC_consensus_E1_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/CCOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "2_Active_Promoter", "0", ".", $2, $3, "255,0,0"}' $New_Folder/consensus_chromHMM_states/state_E2/CCOC_consensus_E2_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/CCOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "3_Active_Region", "0", ".", $2, $3, "231,84,128"}' $New_Folder/consensus_chromHMM_states/state_E3/CCOC_consensus_E3_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/CCOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "4_Active_Enhancer", "0", ".", $2, $3, "250,202,0"}' $New_Folder/consensus_chromHMM_states/state_E4/CCOC_consensus_E4_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/CCOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "5_Weak_Enhancer", "0", ".", $2, $3, "255,252,4"}' $New_Folder/consensus_chromHMM_states/state_E5/CCOC_consensus_E5_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/CCOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "6_Insulator", "0", ".", $2, $3, "10,190,254"}' $New_Folder/consensus_chromHMM_states/state_E6/CCOC_consensus_E6_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/CCOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "7_Transcribed", "0", ".", $2, $3, "0,176,80"}' $New_Folder/consensus_chromHMM_states/state_E7/CCOC_consensus_E7_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/CCOC_consensus_ChromHMM.bed

### EnOC ###

echo "track name=\"EnOC_consensus_ChromHMM\" description=\"ChromHMM Output for Clear Cell Ovarian Cancer Lines Consensus - Colors: State 1 (Light Red) - Weak Promoter, State 2 (Bright Red) - Active Promoter, State 3 (Purple Pink) - Active Region, State 4 (Orange Yellow) - Active Enhancer, State 5 (Bright Yellow) - Weak Enhancer, State 6 (Light Blue) - Insulator, State 7 (Green) - Transcribed\" visibility=1 itemRgb=\"on\"" >| $New_Folder/UCSC_hg38_BED_Files/EnOC_consensus_ChromHMM.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "1_Weak_Promoter", "0", ".", $2, $3, "255,105,105"}' $New_Folder/consensus_chromHMM_states/state_E1/EnOC_consensus_E1_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/EnOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "2_Active_Promoter", "0", ".", $2, $3, "255,0,0"}' $New_Folder/consensus_chromHMM_states/state_E2/EnOC_consensus_E2_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/EnOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "3_Active_Region", "0", ".", $2, $3, "231,84,128"}' $New_Folder/consensus_chromHMM_states/state_E3/EnOC_consensus_E3_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/EnOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "4_Active_Enhancer", "0", ".", $2, $3, "250,202,0"}' $New_Folder/consensus_chromHMM_states/state_E4/EnOC_consensus_E4_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/EnOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "5_Weak_Enhancer", "0", ".", $2, $3, "255,252,4"}' $New_Folder/consensus_chromHMM_states/state_E5/EnOC_consensus_E5_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/EnOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "6_Insulator", "0", ".", $2, $3, "10,190,254"}' $New_Folder/consensus_chromHMM_states/state_E6/EnOC_consensus_E6_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/EnOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "7_Transcribed", "0", ".", $2, $3, "0,176,80"}' $New_Folder/consensus_chromHMM_states/state_E7/EnOC_consensus_E7_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/EnOC_consensus_ChromHMM.bed

### FT ###

echo "track name=\"FT_consensus_ChromHMM\" description=\"ChromHMM Output for Clear Cell Ovarian Cancer Lines Consensus - Colors: State 1 (Light Red) - Weak Promoter, State 2 (Bright Red) - Active Promoter, State 3 (Purple Pink) - Active Region, State 4 (Orange Yellow) - Active Enhancer, State 5 (Bright Yellow) - Weak Enhancer, State 6 (Light Blue) - Insulator, State 7 (Green) - Transcribed\" visibility=1 itemRgb=\"on\"" >| $New_Folder/UCSC_hg38_BED_Files/FT_consensus_ChromHMM.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "1_Weak_Promoter", "0", ".", $2, $3, "255,105,105"}' $New_Folder/consensus_chromHMM_states/state_E1/FT_consensus_E1_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/FT_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "2_Active_Promoter", "0", ".", $2, $3, "255,0,0"}' $New_Folder/consensus_chromHMM_states/state_E2/FT_consensus_E2_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/FT_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "3_Active_Region", "0", ".", $2, $3, "231,84,128"}' $New_Folder/consensus_chromHMM_states/state_E3/FT_consensus_E3_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/FT_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "4_Active_Enhancer", "0", ".", $2, $3, "250,202,0"}' $New_Folder/consensus_chromHMM_states/state_E4/FT_consensus_E4_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/FT_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "5_Weak_Enhancer", "0", ".", $2, $3, "255,252,4"}' $New_Folder/consensus_chromHMM_states/state_E5/FT_consensus_E5_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/FT_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "6_Insulator", "0", ".", $2, $3, "10,190,254"}' $New_Folder/consensus_chromHMM_states/state_E6/FT_consensus_E6_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/FT_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "7_Transcribed", "0", ".", $2, $3, "0,176,80"}' $New_Folder/consensus_chromHMM_states/state_E7/FT_consensus_E7_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/FT_consensus_ChromHMM.bed

### HGSOC ###

echo "track name=\"HGSOC_consensus_ChromHMM\" description=\"ChromHMM Output for Clear Cell Ovarian Cancer Lines Consensus - Colors: State 1 (Light Red) - Weak Promoter, State 2 (Bright Red) - Active Promoter, State 3 (Purple Pink) - Active Region, State 4 (Orange Yellow) - Active Enhancer, State 5 (Bright Yellow) - Weak Enhancer, State 6 (Light Blue) - Insulator, State 7 (Green) - Transcribed\" visibility=1 itemRgb=\"on\"" >| $New_Folder/UCSC_hg38_BED_Files/HGSOC_consensus_ChromHMM.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "1_Weak_Promoter", "0", ".", $2, $3, "255,105,105"}' $New_Folder/consensus_chromHMM_states/state_E1/HGSOC_consensus_E1_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/HGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "2_Active_Promoter", "0", ".", $2, $3, "255,0,0"}' $New_Folder/consensus_chromHMM_states/state_E2/HGSOC_consensus_E2_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/HGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "3_Active_Region", "0", ".", $2, $3, "231,84,128"}' $New_Folder/consensus_chromHMM_states/state_E3/HGSOC_consensus_E3_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/HGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "4_Active_Enhancer", "0", ".", $2, $3, "250,202,0"}' $New_Folder/consensus_chromHMM_states/state_E4/HGSOC_consensus_E4_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/HGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "5_Weak_Enhancer", "0", ".", $2, $3, "255,252,4"}' $New_Folder/consensus_chromHMM_states/state_E5/HGSOC_consensus_E5_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/HGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "6_Insulator", "0", ".", $2, $3, "10,190,254"}' $New_Folder/consensus_chromHMM_states/state_E6/HGSOC_consensus_E6_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/HGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "7_Transcribed", "0", ".", $2, $3, "0,176,80"}' $New_Folder/consensus_chromHMM_states/state_E7/HGSOC_consensus_E7_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/HGSOC_consensus_ChromHMM.bed

### IOSE ###

echo "track name=\"IOSE_consensus_ChromHMM\" description=\"ChromHMM Output for Clear Cell Ovarian Cancer Lines Consensus - Colors: State 1 (Light Red) - Weak Promoter, State 2 (Bright Red) - Active Promoter, State 3 (Purple Pink) - Active Region, State 4 (Orange Yellow) - Active Enhancer, State 5 (Bright Yellow) - Weak Enhancer, State 6 (Light Blue) - Insulator, State 7 (Green) - Transcribed\" visibility=1 itemRgb=\"on\"" >| $New_Folder/UCSC_hg38_BED_Files/IOSE_consensus_ChromHMM.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "1_Weak_Promoter", "0", ".", $2, $3, "255,105,105"}' $New_Folder/consensus_chromHMM_states/state_E1/IOSE_consensus_E1_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/IOSE_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "2_Active_Promoter", "0", ".", $2, $3, "255,0,0"}' $New_Folder/consensus_chromHMM_states/state_E2/IOSE_consensus_E2_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/IOSE_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "3_Active_Region", "0", ".", $2, $3, "231,84,128"}' $New_Folder/consensus_chromHMM_states/state_E3/IOSE_consensus_E3_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/IOSE_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "4_Active_Enhancer", "0", ".", $2, $3, "250,202,0"}' $New_Folder/consensus_chromHMM_states/state_E4/IOSE_consensus_E4_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/IOSE_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "5_Weak_Enhancer", "0", ".", $2, $3, "255,252,4"}' $New_Folder/consensus_chromHMM_states/state_E5/IOSE_consensus_E5_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/IOSE_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "6_Insulator", "0", ".", $2, $3, "10,190,254"}' $New_Folder/consensus_chromHMM_states/state_E6/IOSE_consensus_E6_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/IOSE_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "7_Transcribed", "0", ".", $2, $3, "0,176,80"}' $New_Folder/consensus_chromHMM_states/state_E7/IOSE_consensus_E7_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/IOSE_consensus_ChromHMM.bed

### LGSOC ###

echo "track name=\"LGSOC_consensus_ChromHMM\" description=\"ChromHMM Output for Clear Cell Ovarian Cancer Lines Consensus - Colors: State 1 (Light Red) - Weak Promoter, State 2 (Bright Red) - Active Promoter, State 3 (Purple Pink) - Active Region, State 4 (Orange Yellow) - Active Enhancer, State 5 (Bright Yellow) - Weak Enhancer, State 6 (Light Blue) - Insulator, State 7 (Green) - Transcribed\" visibility=1 itemRgb=\"on\"" >| $New_Folder/UCSC_hg38_BED_Files/LGSOC_consensus_ChromHMM.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "1_Weak_Promoter", "0", ".", $2, $3, "255,105,105"}' $New_Folder/consensus_chromHMM_states/state_E1/LGSOC_consensus_E1_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/LGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "2_Active_Promoter", "0", ".", $2, $3, "255,0,0"}' $New_Folder/consensus_chromHMM_states/state_E2/LGSOC_consensus_E2_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/LGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "3_Active_Region", "0", ".", $2, $3, "231,84,128"}' $New_Folder/consensus_chromHMM_states/state_E3/LGSOC_consensus_E3_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/LGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "4_Active_Enhancer", "0", ".", $2, $3, "250,202,0"}' $New_Folder/consensus_chromHMM_states/state_E4/LGSOC_consensus_E4_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/LGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "5_Weak_Enhancer", "0", ".", $2, $3, "255,252,4"}' $New_Folder/consensus_chromHMM_states/state_E5/LGSOC_consensus_E5_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/LGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "6_Insulator", "0", ".", $2, $3, "10,190,254"}' $New_Folder/consensus_chromHMM_states/state_E6/LGSOC_consensus_E6_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/LGSOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "7_Transcribed", "0", ".", $2, $3, "0,176,80"}' $New_Folder/consensus_chromHMM_states/state_E7/LGSOC_consensus_E7_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/LGSOC_consensus_ChromHMM.bed

### MOC ###

echo "track name=\"MOC_consensus_ChromHMM\" description=\"ChromHMM Output for Clear Cell Ovarian Cancer Lines Consensus - Colors: State 1 (Light Red) - Weak Promoter, State 2 (Bright Red) - Active Promoter, State 3 (Purple Pink) - Active Region, State 4 (Orange Yellow) - Active Enhancer, State 5 (Bright Yellow) - Weak Enhancer, State 6 (Light Blue) - Insulator, State 7 (Green) - Transcribed\" visibility=1 itemRgb=\"on\"" >| $New_Folder/UCSC_hg38_BED_Files/MOC_consensus_ChromHMM.bed

awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "1_Weak_Promoter", "0", ".", $2, $3, "255,105,105"}' $New_Folder/consensus_chromHMM_states/state_E1/MOC_consensus_E1_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/MOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "2_Active_Promoter", "0", ".", $2, $3, "255,0,0"}' $New_Folder/consensus_chromHMM_states/state_E2/MOC_consensus_E2_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/MOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "3_Active_Region", "0", ".", $2, $3, "231,84,128"}' $New_Folder/consensus_chromHMM_states/state_E3/MOC_consensus_E3_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/MOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "4_Active_Enhancer", "0", ".", $2, $3, "250,202,0"}' $New_Folder/consensus_chromHMM_states/state_E4/MOC_consensus_E4_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/MOC_consensus_ChromHMM.bed
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "5_Weak_Enhancer", "0", ".", $2, $3, "255,252,4"}' $New_Folder/consensus_chromHMM_states/state_E5/MOC_consensus_E5_hg38.bed >> $New_Folder/UCSC_hg38_BED_Files/MOC_consensus_ChromHMM.bed

### Combine ChromHMM ###

cat $New_Folder/UCSC_hg38_BED_Files/CCOC_consensus_ChromHMM.bed $New_Folder/UCSC_hg38_BED_Files/EnOC_consensus_ChromHMM.bed $New_Folder/UCSC_hg38_BED_Files/FT_consensus_ChromHMM.bed $New_Folder/UCSC_hg38_BED_Files/HGSOC_consensus_ChromHMM.bed $New_Folder/UCSC_hg38_BED_Files/LGSOC_consensus_ChromHMM.bed $New_Folder/UCSC_hg38_BED_Files/MOC_consensus_ChromHMM.bed > $New_Folder/UCSC_hg38_BED_Files/All_consensus_ChromHMM.bed



########################################
### Overall Case Control Unique CNVs ###
########################################

head -n 3 $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt
head -n 3 $New_Folder/ocac_overall_cnvcalls_IDs_hg38_OncIDcounts.txt

awk '(FNR!=1){print $1, $2, $3, $5, $9}' $New_Folder/ocac_overall_cnvcalls_IDs_hg38.txt | sort | uniq -c | awk 'BEGIN{OFS="\t"; print "chr", "start", "end", "del_dup", "case_ctrl", "count"} {print $2, $3, $4, $5, $6, $1}' >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.txt
head -n 3 $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.txt

awk '(FNR!=1){print $6}' $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.txt | sort | uniq -c

awk -v genesymbol="BRCA1" 'BEGIN{OFS="\t"} $4==genesymbol{print}' $CNVs_Folder/Biomart_All_Genes_Hg38.bed >| $New_Folder/BRCA1_geneloc.bed
bedtools intersect -wa -a $New_Folder/ocac_overall_cnvcalls_IDs_hg38.bed -b $New_Folder/BRCA1_geneloc.bed | awk '{gsub(";", "\t"); print}' | wc -l
bedtools intersect -wa -a $New_Folder/ocac_overall_cnvcalls_IDs_hg38.bed -b $New_Folder/BRCA1_geneloc.bed | awk '{gsub(";", "\t"); print}' | awk '($5=="del" && $9=="0")' | wc -l
bedtools intersect -wa -a $New_Folder/ocac_overall_cnvcalls_IDs_hg38.bed -b $New_Folder/BRCA1_geneloc.bed | awk '{gsub(";", "\t"); print}' | awk '($5=="del" && $9=="1")' | wc -l
bedtools intersect -wa -a $New_Folder/ocac_overall_cnvcalls_IDs_hg38.bed -b $New_Folder/BRCA1_geneloc.bed | awk '{gsub(";", "\t"); print}' | awk '($5=="dup" && $9=="0")' | wc -l
bedtools intersect -wa -a $New_Folder/ocac_overall_cnvcalls_IDs_hg38.bed -b $New_Folder/BRCA1_geneloc.bed | awk '{gsub(";", "\t"); print}' | awk '($5=="dup" && $9=="1")' | wc -l

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"Overall_CaseControl_Unique_CNVs\" description=\"CNV Segments from All Cases and Controls colored by deletion/duplication (red/blue) and case/control status (dark/light), name is Count, color is red for dels, blue for dups, dark for case, light for control\" itemRgb=\"on\""; next} ($4=="del" && $5=="1"){print $1, $2, $3, $6, 0, ".", $2, $3, "255,0,0,"; next} ($4=="del" && $5=="0"){print $1, $2, $3, $6, 0, ".", $2, $3, "255,204,204,"; next} ($4=="dup" && $5=="1"){print $1, $2, $3, $6, 0, ".", $2, $3, "0,0,255,"; next} ($4=="dup" && $5=="0"){print $1, $2, $3, $6, 0, ".", $2, $3, "204,204,255,"; next}' $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38_UCSC.bed

head -n 5 $New_Folder/UCSC_hg38_BED_Files/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38_UCSC.bed


######################################
### HGSOC Case Control Unique CNVs ###
######################################

head -n 3 $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt
head -n 3 $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38_OncIDcounts.txt

awk '(FNR!=1){print $1, $2, $3, $5, $9}' $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.txt | sort | uniq -c | awk 'BEGIN{OFS="\t"; print "chr", "start", "end", "del_dup", "case_ctrl", "count"} {print $2, $3, $4, $5, $6, $1}' >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.txt
head -n 3 $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.txt

# BRCA1: chr17	43044294	43170245
awk -v genesymbol="BRCA1" 'BEGIN{OFS="\t"} $4==genesymbol{print}' $CNVs_Folder/Biomart_All_Genes_Hg38.bed >| $New_Folder/BRCA1_geneloc.bed
bedtools intersect -wa -a $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.bed -b $New_Folder/BRCA1_geneloc.bed | awk '{gsub(";", "\t"); print}' | awk '($5=="del" && $9=="0")' | wc -l
bedtools intersect -wa -a $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.bed -b $New_Folder/BRCA1_geneloc.bed | awk '{gsub(";", "\t"); print}' | awk '($5=="del" && $9=="1")' | wc -l
bedtools intersect -wa -a $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.bed -b $New_Folder/BRCA1_geneloc.bed | awk '{gsub(";", "\t"); print}' | awk '($5=="dup" && $9=="0")' | wc -l
bedtools intersect -wa -a $New_Folder/ocac_hgsoc_cnvcalls_IDs_hg38.bed -b $New_Folder/BRCA1_geneloc.bed | awk '{gsub(";", "\t"); print}' | awk '($5=="dup" && $9=="1")' | wc -l

awk '(FNR!=1){print $6}' $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.txt | sort | uniq -c

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"HGSOC_CaseControl_Unique_CNVs\" description=\"CNV Segments from HGSOC Cases and Controls colored by deletion/duplication (red/blue) and case/control status (dark/light), name is Count, color is red for dels, blue for dups, dark for case, light for control\" itemRgb=\"on\""; next} ($4=="del" && $5=="1"){print $1, $2, $3, $6, 0, ".", $2, $3, "255,0,0,"; next} ($4=="del" && $5=="0"){print $1, $2, $3, $6, 0, ".", $2, $3, "255,204,204,"; next} ($4=="dup" && $5=="1"){print $1, $2, $3, $6, 0, ".", $2, $3, "0,0,255,"; next} ($4=="dup" && $5=="0"){print $1, $2, $3, $6, 0, ".", $2, $3, "204,204,255,"; next}' $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38_UCSC.bed

head -n 5 $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38_UCSC.bed


#######################################################
### Overall Case Control Unique Sig (p < 0.05) CNVs ###
#######################################################


awk 'BEGIN{OFS="\t"} (FNR!=1){print $1, $2, $3, $4";"$5";"$6}' $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.bed

awk 'BEGIN{OFS="\t"} (FNR==1){print $1, $2, $3, $4, $5, $6}' $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.header

bedtools intersect -loj -wb -wa -a $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.bed -b $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.bed >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_notab.txt

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_notab.txt >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_nohead.txt


cat $New_Folder/ocac_overall_cnv_results_v7_sig05_probes_hg38.header $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_hg38.header | gsed '$!{:a;N;s/\r//;s/\n/\t/;ta}' >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38.header

cat $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38.header $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_nohead.txt >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_mismatch.txt

awk '(FNR==1||$18=="."){print $0}' $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_nomatch.txt

rm $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_notab.txt
rm $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_nohead.txt

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$5==$21){print $0}' $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38.txt

awk 'BEGIN{OFS="\t"; FS="\t"} (!seen[$18";"$19";"$20";"$21";"$22";"$23]++){print $18, $19, $20, $21, $22, $23}' $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38.txt >| $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_unique_hg38.txt

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"Overall_CaseControl_Unique_CNVs\" description=\"Significant (p < 0.05) CNV Segments from Overall Cases and Controls colored by deletion/duplication (red/blue) and case/control status (dark/light), name is Count, color is red for dels, blue for dups, dark for case, light for control\" itemRgb=\"on\""; next} ($4=="del" && $5=="1"){print $1, $2, $3, $6, 0, ".", $2, $3, "255,0,0,"; next} ($4=="del" && $5=="0"){print $1, $2, $3, $6, 0, ".", $2, $3, "255,204,204,"; next} ($4=="dup" && $5=="1"){print $1, $2, $3, $6, 0, ".", $2, $3, "0,0,255,"; next} ($4=="dup" && $5=="0"){print $1, $2, $3, $6, 0, ".", $2, $3, "204,204,255,"; next}' $New_Folder/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_unique_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_overall_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_UCSC.bed



#######################################################
### HGSOC Case Control Unique Sig (p < 0.05) CNVs ###
#######################################################

awk 'BEGIN{OFS="\t"} (FNR!=1){print $1, $2, $3, $4";"$5";"$6}' $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.bed

awk 'BEGIN{OFS="\t"} (FNR==1){print $1, $2, $3, $4, $5, $6}' $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.header

bedtools intersect -loj -wb -wa -a $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.bed -b $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.bed >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_notab.txt

awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/;/, "\t"); print}' $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_notab.txt >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_nohead.txt


cat $New_Folder/ocac_hgsoc_cnv_results_v7_sig05_probes_hg38.header $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_hg38.header | gsed '$!{:a;N;s/\r//;s/\n/\t/;ta}' >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38.header


cat $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38.header $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_nohead.txt >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_mismatch.txt

awk '(FNR==1||$18=="."){print $0}' $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_nomatch.txt

rm $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_notab.txt
rm $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_nohead.txt


awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1||$5==$21){print $0}' $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_mismatch.txt >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38.txt

awk 'BEGIN{OFS="\t"} (seen[$18";"$19";"$20";"$21";"$22";"$23]++){print $18, $19, $20, $21, $22, $23}' $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38.txt >| $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_unique_hg38.txt

awk 'BEGIN{FS="\t"; OFS="\t"} (FNR==1){print "track name=\"HGSOC_CaseControl_Unique_CNVs\" description=\"Significant (p < 0.05) CNV Segments from HGSOC Cases and Controls colored by deletion/duplication (red/blue) and case/control status (dark/light), name is Count, color is red for dels, blue for dups, dark for case, light for control\" itemRgb=\"on\""; next} ($4=="del" && $5=="1"){print $1, $2, $3, $6, 0, ".", $2, $3, "255,0,0,"; next} ($4=="del" && $5=="0"){print $1, $2, $3, $6, 0, ".", $2, $3, "255,204,204,"; next} ($4=="dup" && $5=="1"){print $1, $2, $3, $6, 0, ".", $2, $3, "0,0,255,"; next} ($4=="dup" && $5=="0"){print $1, $2, $3, $6, 0, ".", $2, $3, "204,204,255,"; next}' $New_Folder/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_unique_hg38.txt >| $New_Folder/UCSC_hg38_BED_Files/ocac_hgsoc_cnvcalls_casecontrol_deldup_counts_sig05_CNVs_hg38_UCSC.bed
