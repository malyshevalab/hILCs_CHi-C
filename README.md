# hILCs_CHi-C
Scripts to accompany Malysheva/Ray-Jones et al., 2022: High-resolution promoter interaction analysis in Type 3 Innate Lymphoid Cells implicates Batten Disease gene CLN3 in Crohn’s Disease aetiology 

### Alternative_isoform_CHi-C_analysis.R
Analysis of alternative isoform interactions in hILC3 PCHiC data. 

### COGS_comparisons_and_tables.ipynb
- Comparison between different COGS runs (classic, classic + imputation and SuSIE, i.e. full multiCOGS). 
- Production of heatmaps comparing the features underlying COGS signals.
- Formatting of supplementary tables
- Generating UCSC browser tracks of data

### compare_COGS_results_scatterplot.R
Make a scatterplot comparing the COGS scores for two condiions; here, standard V multiCOGS

### make_manhattan_oneSample.R
Generate Manhattan plot for COGS scores per gene

### find_novel_Crohns_genes.R
Look for prior evidence for the prioritised genes in Crohn's Disease multiCOGS in ILC3s and CD4+ T cells

### COGS_clusterProfiler_final.R
Run GO term analysis on COGS/multiCOGS results in Crohn's Disease.

### COGS_clusterProfiler_crossTrait.R
Run GO term analysis on multiCOGS results in additional autoimmune traits.

### multiCOGS_cross_trait_comparison.ipynb
- Make heatmap of multiCOGS scores across autoimmune traits
- Compare and plot GO enrichment across traits
- Make supplementary tables

### disease_ontology.R
Explore the RELI disease enrichment results in ILC3 PIRs and look for overrepresentation of immune diseases

### run_Mageck_pathway.sh
Run the pathway analysis using Mageck to determine enrichment of multiCOGS genes across autoimmune traits (ILC3s) among genes in the IL-22 MNK3 CRISPRi screen (Brown et al., Biorxiv 2025)

### TODO: script for plotting the mageck results.


