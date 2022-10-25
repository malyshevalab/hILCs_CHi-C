# hILCs_CHi-C
Scripts for analysing hILCs CHi-C, including COGS runs and downstream analyses

## rCOGS scripts
These scripts are generic scripts for creating the formatted input files and running rCOGS (the data.table version and multiCOGS). These scripts are kept in the same directory.

### 1. Make_rCOGS_input_files.sh
Script that prepares COGS inputs from given file locations.
### 2. annot_CHi-C_files.R
Run by Make_rCOGS_input_files.sh. Prepares the annotated baitmap and peak matrix.
### 3. run_vep.sh
Runs vep to get coding SNPs for a given set of rsids, produced from running Make_rCOGS_input_files.sh
### 4. run_Mikhails_rCOGS.R
Script for loading input files and running either Standard or multiCOGS. Inlcudes modifications for SuSIE input. Will work with input files produced by scripts above.

## rCOGS runs
Contains the sh scripts used to generate the COGS results in the paper.

## Other analysis scripts

### Alt_isoform_analysis_HRJ.R
Analysis of alternative isoform interactions in hILC PCHiC data.

### COGS_score_contributions_PIR_analysis.Rmd
Comparison between ILC3 COGS genes and other blood type COGS genes that required CHi-C PIRs for prioritisation.

### COGS_clusterProfiler_final.R
Run GO term / GSEA analysis on COGS results.

### compare_COGS_results_scatterplot.R
Make a scatterplot comparing the COGS scores for two condiions; here, standard V multiCOGS

### find_novel_ILC3_genes.R
Look for prior evidence for the prioritised genes based on ILC3 PCHi-C multiCOGS

### disease_ontology.R
Explore the RELI disease enrichment results in ILC3 PIRs and look for overrepresentation of immune diseases

### liftover_javierre_peakm_and_foverlap_Design_files.R
Liftover Javierre peak matrix to hg38 for use in COGS with hg38 GWAS data

### make_manhattan_oneSample.R
Generate Manhattan plot for COGS scores per gene
