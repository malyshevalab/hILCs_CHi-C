# hILCs_CHi-C
Scripts to accompany Malysheva/Ray-Jones et al., 2022: High-resolution promoter interaction analysis in Type 3 Innate Lymphoid Cells implicates Batten Disease gene CLN3 in Crohn’s Disease aetiology 

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
