### Run enrichGO and GSEA on our genes.

## Run in clusterProfiler environment
setwd("~/COGS_results")
library(clusterProfiler)
library(data.table)
library(org.Hs.eg.db)
library(enrichplot)
library(cowplot)
library(forcats)
library(ggplot2)
library(dplyr)

############################################################################
####### Find enriched biological functions among our top COGS genes ########
############################################################################

#### Function for running enrichGO from COGs annotated table (not used)
COGS_enrichGO <- function(cogs, score=0.5) {
  sig <- cogs[cogs >= score]
  cg <- unique(sig$ensg)
  GO <- enrichGO(cg, OrgDb = "org.Hs.eg.db", pvalueCutoff=1, qvalueCutoff=1, 
                 readable = TRUE, keyType = "ENSEMBL")
  print(head(GO, 10))
  print(dotplot(GO))
  print(goplot(GO))
  return(GO)
}
####


#### Function for running gene set enrichment analysis

# For Hallmark sets (regular GO) we can just use the gseGO function. For the curated datasets we need to use GSEA function (see below).

cogs_gseGO <- function(cogs) { 
  # make a file with ENSG and score
  # Order ranked genelist
  dt <- copy(cogs)
  ## feature 1: numeric vector
  geneList = dt[,cogs]
  ## feature 2: named vector
  names(geneList) = as.character(dt[,ensg])
  ## feature 3: decreasing orde
  geneList = sort(geneList, decreasing = TRUE)
  GO <- gseGO(geneList, OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, keyType = "ENSEMBL")
  # The results are ordered by Pvalue
  p1 <- gseaplot(GO, geneSetID = GO$ID[1], title = GO$Description[1])
  p2 <- gseaplot(GO, geneSetID = GO$ID[2], title = GO$Description[2])
  p3 <- gseaplot(GO, geneSetID = GO$ID[3], title = GO$Description[3])
  print(cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3]))
  return(GO)
}

####################################

# deLange classic COGS
cad <- fread("./COGS_out/CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_ALL/Annotated_COGS_scores_data.table.txt")
gsea_deLange <- cogs_gseGO(cad)
gsea_deLange_dt <- as.data.table(gsea_deLange@result)
fwrite(gsea_deLange_dt, file = "./COGS_out/CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_ALL/gseGO_hallmarks.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)
pdf(file = "./COGS_out/CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_ALL/gseGO_deLange_hallmarks_gseaPlotALL.pdf", 
    width = 10)
ggplot(gsea_deLange, showCategory=50, aes(NES, Description, fill=p.adjust)) + 
  geom_col() + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"), guide=guide_colorbar(reverse=TRUE)) +
  xlab("Normalized Enrichment Score") + ylab(NULL) +
  ggtitle("GSEA")
dev.off()

# deLange SuSIE
cas <- fread("./COGS_out/CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/Annotated_COGS_scores_data.table.txt")
gsea_deLange_SuSIE <- cogs_gseGO(cas)
gsea_deLange_SuSIE_dt <- as.data.table(gsea_deLange_SuSIE@result)
fwrite(gsea_deLange_SuSIE_dt, file = "./COGS_out/CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/gseGO_deLange_SuSIE_hallmarks.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)
pdf(file = "./COGS_out/CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/gseGO_deLange_SuSIE_hallmarks_gseaPlotALL.pdf", 
    width = 15, height = 20)
ggplot(gsea_deLange_SuSIE, showCategory=200, aes(NES, Description, fill=p.adjust)) + 
  geom_col() + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"), guide=guide_colorbar(reverse=TRUE)) +
  xlab("Normalized Enrichment Score") + ylab(NULL) +
  ggtitle("GSEA")
dev.off()

# Now do the top 50.

pdf(file = "./COGS_out/CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/gseGO_deLange_SuSIE_hallmarks_gseaPlotTOP50.pdf", 
    width = 30, height = 20)
ggplot(gsea_deLange_SuSIE, showCategory=50, aes(NES, Description, fill=p.adjust)) + 
  geom_col() + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"), guide=guide_colorbar(reverse=TRUE)) +
  xlab("Normalized Enrichment Score") + ylab(NULL) + theme(text=element_text(size = 30)) +
  ggtitle("GSEA")
dev.off()

# top 20
pdf(file = "./COGS_out/CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/gseGO_deLange_SuSIE_hallmarks_gseaPlotTOP20.pdf", 
    width = 24, height = 10)
ggplot(gsea_deLange_SuSIE, showCategory=20, aes(NES, Description, fill=p.adjust)) + 
  geom_col() + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"), guide=guide_colorbar(reverse=TRUE)) +
  xlab("Normalized Enrichment Score") + ylab(NULL) + theme(text=element_text(size = 30)) +
  ggtitle("GSEA")
dev.off()

# Compare with and without SuSIE
nosus_small <- gsea_deLange_dt[, .(ID, Description, pvalue)]
names(nosus_small)[3] = "p.value_noSuSIE"
sus_small <- gsea_deLange_SuSIE_dt[, .(ID, Description, pvalue)]
names(sus_small)[3] = "p.value_SuSIE"
both <- merge.data.table(nosus_small, sus_small, on = c("ID", "Description"), all = TRUE)
susie_only <- both[is.na(p.value_noSuSIE)]
fwrite(susie_only, file = "./comparisons/GO_terms_susie_specific_deLange.txt", sep = "\t", 
       quote = F, col.names = T, row.names = F)


### Run GSE GO against gene sets from known CD cell types
setdir <- "~/external_data/IBD_TaMMA/data/human/dge"

sets <- data.table()
files <- list.files(setdir, pattern = ".diffexp.tsv$")
setwd(setdir)

## Get significant DE genes
for(file in files) {
  geneset <- fread(file)
  sig <- geneset[padj < 0.05 & abs(log2FoldChange) >= 2]
  myname <- sub(".diffexp.tsv", "", file)
  sig[, set := myname]
  keep <- sig[, .(Geneid, Gene, set)]
  sets <- rbind(sets, keep, fill = TRUE)
}

fwrite(sets, file = "allDGE_HRJ.txt", sep = "\t", quote = F, row.names = F, col.names = T)

cd_sets <- fread("allDGE_HRJ.txt")
cd_sets <- cd_sets[, .(set, Geneid)]

### We want to remove sets that compare across different tissues, and remove large ones (> 2000 genes)
cd_sets2 <- as.data.table(cd_sets %>% group_by(set) %>% tally())
keep <- cd_sets2[n <= 2000]
keep_same <- keep[set %like% "Blood" | set == "Colon_CD-vs-Colon_Control" | set %like% "Colon_epithelium" |
                set %like% "Colon_fibroblast" | set == "Colon_UC-vs-Colon_CD" | set == "Colon_UC-vs-Colon_Control" | 
                set == "Colonoid_PIBD-vs-Colonoid_Control" | set == "Ileum_CD-vs-Ileum_Control" | set == "Ileum_CD-vs-Ileum_UC" | 
                set == "Ileum_UC-vs-Ileum_Control" | set %like% "Mesenteric_fat" | set == "Rectum_CD-vs-Rectum_Control" | 
                set == "Rectum_UC-vs-Rectum_CD" | set == "Rectum_UC-vs-Rectum_Control" | set %like% "Stool_associated"]
# Now 24 sets
keep_same_sets <- cd_sets[set %in% keep_same$set]
fwrite(keep_same_sets, file = "~/external_data/IBD_TaMMA/data/human/dge/sameTissueDGE_HRJ.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)

#############
### Function for running GSEA using CD sets.

cogs_GSEA <- function(dt, dataset) { 
  # make a file with ENSG and score
  # Order ranked genelist
  ## feature 1: numeric vector
  geneList = dt[,cogs]
  ## feature 2: named vector
  names(geneList) = as.character(dt[,ensg])
  ## feature 3: decreasing order
  geneList = sort(geneList, decreasing = TRUE)
  
  ## Dataset is TERM2GENE set (2 cols with term and gene)
  
  GO <- GSEA(geneList, pvalueCutoff=10, TERM2GENE = dataset, maxGSSize = 2000) # we actually only test 49 sets with this...
  # 500 genes is the default and is appropriate for 10-20K features so 2K is actually high
  # The results are ordered by Pvalue
  print(gseaplot2(GO, geneSetID = GO$ID[1:3], ES_geom = "dot", pvalue_table = TRUE))
  return(GO)
}
##############
setwd("~/COGS_results")

gsea_susie <- cogs_GSEA(cas, keep_same_sets)
gsea_susie_full <- as.data.table(gsea_susie@result)
ggplot(gsea_susie, showCategory=50, aes(NES, Description, fill=p.adjust)) + 
  geom_col() + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"), guide=guide_colorbar(reverse=TRUE)) +
  xlab("Normalized Enrichment Score") + ylab(NULL) +
  ggtitle("GSEA")
fwrite(gsea_susie_full, file = "./COGS_out/CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/cogsGSEA_CD_deLange_SuSIE_results.txt", 
       sep = "\t", quote = F, col.names = T, row.names = F)
# print the gseaplot2 only for the significant ones
pdf(file = "./COGS_out/CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/cogsGSEA_CD_deLange_SuSIE_rankedPlot_sig.pdf")
print(gseaplot2(gsea_susie, geneSetID = gsea_susie$ID[1], ES_geom = "dot", pvalue_table = TRUE))
dev.off()

### Are any significant without SuSIE?
gsea_deLange_CD <- cogs_GSEA(cad, keep_same_sets)
gsea_deLange_CD_full <- as.data.table(gsea_deLange_CD@result)
fwrite(gsea_deLange_CD_full, file = "./COGS_out/CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_ALL/cogsGSEA_CD_results.txt", 
       sep = "\t", quote = F, col.names = T, row.names = F)

# None were significant

####### 