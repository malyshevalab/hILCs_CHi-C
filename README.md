# hILCs_CHi-C
Scripts to accompany Malysheva/Ray-Jones et al., 2022: High-resolution promoter interaction analysis in Type 3 Innate Lymphoid Cells implicates Batten Disease gene CLN3 in Crohn’s Disease aetiology 

### Alternative_isoform_CHi-C_analysis.R
Analysis of alternative isoform interactions in hILC3 PCHiC data. 

### classicCOGS_vs_impute_vsSuSIE.ipynb
- Comparison between different COGS runs (classic, classic + imputation and SuSIE, i.e. full multiCOGS). 
- Production of heatmaps comparing the features underlying COGS signals.
- Formatting of supplementary tables
- Generating UCSC browser tracks of data

### tocheck still - all below
### Also - 
### comparison with crispr screen?
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
