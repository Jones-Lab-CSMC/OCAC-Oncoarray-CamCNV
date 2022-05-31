#!/bin/bash
#
# runFunciVarMarks_NewQC.sh


./FunciVar_Submit_NonCoding_sig05_CNVs_Overall_ChromHMM.R --mark "E1"
./FunciVar_Submit_NonCoding_sig05_CNVs_Overall_ChromHMM.R --mark "E2"
./FunciVar_Submit_NonCoding_sig05_CNVs_Overall_ChromHMM.R --mark "E3"
./FunciVar_Submit_NonCoding_sig05_CNVs_Overall_ChromHMM.R --mark "E4"
./FunciVar_Submit_NonCoding_sig05_CNVs_Overall_ChromHMM.R --mark "E5"
./FunciVar_Submit_NonCoding_sig05_CNVs_Overall_ChromHMM.R --mark "E6"
./FunciVar_Submit_NonCoding_sig05_CNVs_Overall_ChromHMM.R --mark "E7"
# ./FunciVar_Submit_NonCoding_sig05_CNVs_Overall_ChromHMM.R --mark "E8"

./FunciVar_Submit_NonCoding_sig05_CNVs_HGSOC_ChromHMM.R --mark "E1"
./FunciVar_Submit_NonCoding_sig05_CNVs_HGSOC_ChromHMM.R --mark "E2"
./FunciVar_Submit_NonCoding_sig05_CNVs_HGSOC_ChromHMM.R --mark "E3"
./FunciVar_Submit_NonCoding_sig05_CNVs_HGSOC_ChromHMM.R --mark "E4"
./FunciVar_Submit_NonCoding_sig05_CNVs_HGSOC_ChromHMM.R --mark "E5"
./FunciVar_Submit_NonCoding_sig05_CNVs_HGSOC_ChromHMM.R --mark "E6"
./FunciVar_Submit_NonCoding_sig05_CNVs_HGSOC_ChromHMM.R --mark "E7"


wait

./CombineFunciVarEnrichment_NewQC_Final.R
