### Run gene set enrichment analysis on multiple traits run with multiCOGS/hILCs. Plot a results overview.

## Run in clusterProfiler environment!!!
setwd("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/COGS_results")
library(clusterProfiler)
library(data.table)
library(org.Hs.eg.db)
library(enrichplot)
library(cowplot)
library(forcats)
library(ggplot2)

####################################

cd <- fread("./COGS_out/Version3_revision2/revision_deLange_ILCs_hg38_SuSIE_combinedInteractions_extended_ABC023/COGS_scores_data.table.txt")
astao <- fread("./COGS_out/ASTAO_Ferreira_30929738/ASTAO_Ferreira_30929738_ILC3_SuSIE/COGS_scores_data.table.txt")
uc <- fread("./COGS_out/UC_DeLange_28067908//UC_DeLange_28067908_ILC3_SuSIE/COGS_scores_data.table.txt")
ibd <- fread("./COGS_out/IBD_DeLange_28067908/IBD_DeLange_28067908_ILC3_SuSIE/COGS_scores_data.table.txt")
cel <- fread("./COGS_out/CEL_Dubois_20190752/CEL_Dubois_20190752_ILC3_SuSIE/COGS_scores_data.table.txt")
psc <- fread("./COGS_out/PSC_Ji_27992413/PSC_Ji_27992413_ILC3_SuSIE/COGS_scores_data.table.txt")

#For the plot:
# y axis is pathways
# x axis is traits
# size of dots is geneRatio
# colour is Padj
# label on the top is total no. of genes analysed for each trait
# Maybe just do this for the most sig enriched pathways
# doing in ipynb (multiCOGS_cross_trait_comparison.ipynb)

# IBD core genes - just do GO term analysis
# The following takes as input a table with the required genes as ensg
COGS_enrichGO <- function(cogs) {
  allGenes <- unique(cogs$ensg)
  cg <- unique(cogs[cogs > 0.5, ensg])
  GO <- enrichGO(cg, OrgDb = "org.Hs.eg.db", universe = allGenes, pvalueCutoff=1, qvalueCutoff=1, 
                 readable = TRUE, keyType = "ENSEMBL")
  #print(head(GO, 10))
  #print(dotplot(GO))
  #print(goplot(GO))
  return(GO)
}

### IBD core, with background
ibd_core <- fread("./cross_trait/IBD_core_COGS_genes.txt")
gsea_ibd_core <- COGS_enrichGO(ibd_core)
gsea_ibd_core_dt <- as.data.table(gsea_ibd_core@result)
fwrite_headers(gsea_ibd_core_dt, file = "./cross_trait/ibd_core_gseGO.txt")

pdf(file = "./cross_trait/goPlot_IBD_core.pdf", width = 8, height = 8)
goplot(gsea_ibd_core)
dev.off()
###

gsea_cel <- COGS_enrichGO(cel)
gsea_cel_dt <- as.data.table(gsea_cel@result)
fwrite_headers(gsea_cel_dt, file = "./cross_trait/cel_gseGO.txt")

gsea_psc <- COGS_enrichGO(psc)
gsea_psc_dt <- as.data.table(gsea_psc@result)
fwrite_headers(gsea_psc_dt, file = "./cross_trait/psc_gseGO.txt")

gsea_cd <- COGS_enrichGO(cd)
gsea_cd_dt <- as.data.table(gsea_cd@result)
fwrite_headers(gsea_cd_dt, file = "./cross_trait/cd_gseGO.txt")

gsea_astao <- COGS_enrichGO(astao)
gsea_astao_dt <- as.data.table(gsea_astao@result)
fwrite_headers(gsea_astao_dt, file = "./cross_trait/astao_gseGO.txt")

gsea_uc <- COGS_enrichGO(uc)
gsea_uc_dt <- as.data.table(gsea_uc@result)
fwrite_headers(gsea_uc_dt, file = "./cross_trait/uc_gseGO.txt")

gsea_ibd <- COGS_enrichGO(ibd)
gsea_ibd_dt <- as.data.table(gsea_ibd@result)
fwrite_headers(gsea_ibd_dt, file = "./cross_trait/ibd_gseGO.txt")

