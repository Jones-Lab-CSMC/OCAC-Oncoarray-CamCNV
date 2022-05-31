# OCAC-Oncoarray-CamCNV

All annotated CNV calls are available in "CNV_Segments_Annotated.txt".

CNV Calling done using the CamCNV pipeline here:
https://github.com/jgd29/CamCNV
Settings used are specified in the supplementary methods.

The by probe association tests were run using the logit program.
http://ccge.medschl.cam.ac.uk/software/mlogit/
There are there no additional options, just:
  -i arg                input file
  -p arg                phenotype file
  -o arg                output file

RAML is available here: http://ccge.medschl.cam.ac.uk/software/raml/
The settings used for RAML were:
raml -r 0.6 -m 0.01 -n 10000
-r arg                r2 at which to merge correlated variants (default 0.9)
-m arg                maximum maf for rare alleles
-n arg                number of permutations (default 10000)

All other main analysis scripts and files for easy running of these scripts are available in "analysis_scripts" and "analysis_files".

The simulated background files for the enrichment analysis are too large for normal Github storage. They can be easily generated with the available uploaded files using "Script_Initial_File_Conversion.sh" if needed.
