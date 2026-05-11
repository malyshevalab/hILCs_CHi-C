library(data.table)
library(tidyr)
library(dplyr)
setwd("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes")

### Now running for both ILC3s and CD4s (make one combined table.)

delange_ILC3 <- fread("./hILCs/COGS_results/COGS_out/Version3_revision2/revision_deLange_ILCs_hg38_SuSIE_combinedInteractions_extended_ABC023/Annotated_COGS_scores_data.table.txt")
delange_CD4s <- fread("./hILCs/COGS_results/COGS_out/Version3_revision2/revision_deLange_CD4s_hg38_SuSIE_fix_combinedInteractions_extended_ABC023/Annotated_COGS_scores_data.table.txt")

delange_ILC3_sig <- delange_ILC3[cogs > 0.5]
delange_CD4s_sig <- delange_CD4s[cogs > 0.5]
delange_ILC3_sig[, cellType := "ILC3"]
delange_CD4s_sig[, cellType := "CD4"]

ot1 <- fread("./hILCs/previous_evidence/GCST004132-independently-associated-loci.tsv")
ot2 <- fread("./hILCs/previous_evidence/GCST003044-independently-associated-loci.tsv")
ot3 <- fread("./hILCs/previous_evidence/GCST000879-independently-associated-loci.tsv")
ot4 <- fread("./hILCs/previous_evidence/GCST000207-independently-associated-loci.tsv")
ot5 <- fread("./hILCs/previous_evidence/GCST005537-independently-associated-loci.tsv")
all_ot <- rbind(ot1, ot2, ot3, ot4, ot5)

new_ILC3 <- delange_ILC3_sig[!gene %in% all_ot$L2G]
new_CD4s <- delange_CD4s_sig[!gene %in% all_ot$L2G]

ibd <- fread("./hILCs/previous_evidence/Inflammtory Bowel Disease.csv")


new_ILC3_2 <- new_ILC3[!ensg %in% ibd$ensembl_id]
new_ILC3_3 <- new_ILC3_2[!gene %in% ibd$`HGNC gene symbol`]
new_ILC3_3

new_CD4s_2 <- new_CD4s[!ensg %in% ibd$ensembl_id]
new_CD4s_3 <- new_CD4s_2[!gene %in% ibd$`HGNC gene symbol`]
new_CD4s_3

exome <- c("NOD2", "IL23R", "TYK2", "PTPN22", "CARD9", "GPR65", "MST1", "NOD2", "SMAD3", 
                             "PTAFR", "ATG16L1", "SAG", "USP40", "HGFAC", "PDLIM5", "SLC39A8", "PTGER4", "IRGM", 
                             "TAGAP", "DOK2", "SNAPC4", "SEC16A", "RELA", "IL10RA", "LRRK2", "SNX20", "CYLD", 
                             "CCR7", "SDF2L1")



new_ILC3_4 <- new_ILC3_3[!gene %in% exome] # now 71
new_CD4s_4 <- new_CD4s_3[!gene %in% exome] # now 74

rioux <- c("HNF4A", "IFIH1", "SMAD3", "SBNO2", "NFKB1", "NOD2", "ZFP36L1", "IRF1", "GIGYF1", "OTUD3", "AIRE", "PITX1", 
           "KSR1", "DUSP16")

new_ILC3_5 <- new_ILC3_4[!gene %in% rioux]
new_CD4s_5 <- new_CD4s_4[!gene %in% rioux]

### OK then get genes from here https://www.disgenet.org/browser/0/1/0/C0021390/ 
### Either looking at CD, IBD or both. Maybe IBD. I think with IBD we will remove quite a few of the hits.

dis_ibd <- fread("./hILCs/previous_evidence/IBD_C0021390_disease_gda_evidences.txt")
dis_ibd_functional <- dis_ibd[Association_Type != "GeneticVariation"]

new_ILC3_6 <- new_ILC3_5[!gene %in% dis_ibd_functional$Gene] # now there are 63
new_CD4s_6 <- new_CD4s_5[!gene %in% dis_ibd_functional$Gene] # 

dis_cd <- fread("./hILCs/previous_evidence/CD_C0010346_disease_gda_evidences.txt")
dis_cd_functional <- dis_cd[Association_Type != "GeneticVariation"]

new_ILC3_7 <- new_ILC3_6[!gene %in% dis_cd_functional$Gene] # 62
new_CD4s_7 <- new_CD4s_6[!gene %in% dis_cd_functional$Gene] # 68

fwrite(new_ILC3_7, file = "./hILCs/previous_evidence/multiCOGS_our_ILC3_new_genes.txt", sep = "\t", quote =F, row.names = F, col.names = T)
fwrite(new_CD4s_7, file = "./hILCs/previous_evidence/multiCOGS_our_CD4s_new_genes.txt", sep = "\t", quote =F, row.names = F, col.names = T)

#jav <- fread("./hILCs/COGS_results/COGS_out/CD_deLange_Javierre_hg19_allCells/Annotated_COGS_scores_data.table_cell_type_summary.txt")
#new8 <- new7[!ensg %in% jav$ensg]

#### Get the annotation supplementary tables. 
#### One with yes/no and one with raw scores.

delange_ILC3_sig_small <- delange_ILC3_sig[, .(ensg, gene, cogs)]
delange_CD4s_sig_small <- delange_CD4s_sig[, .(ensg, gene, cogs)]
setnames(delange_ILC3_sig_small, "cogs", "ILC3")
setnames(delange_CD4s_sig_small, "cogs", "CD4")

# Make a total join.
both <- merge.data.table(delange_ILC3_sig_small, delange_CD4s_sig_small, by = c("ensg", "gene"), all = T)

options(scipen = 999)

# Make a column of max score
# The following assumes that there are two first cols (ensg and gene name)
myCells <- c("ILC3", "CD4") 
both[, maxScore := do.call(pmax, c(.SD, na.rm = TRUE)),by=c("ensg", "gene"), .SDcols=myCells] # pmax is the parallel maxima.

### Get the name of the column with the max score. In case of ties, paste them.
# R sometimes rounds numbers under the hood. So, the following code handles:
# True zeros (0) vs. small numbers (6e-12) won't be confused
# Tolerates tiny floating-point noise (1.00000000001 vs. 1.0)
# Still works when max is very close to zero (switches to absolute check)
both[, maxCell := apply(.SD, 1, function(row) {
  names(row) <- myCells
  if (all(is.na(row))) return(NA_character_)
  
  max_val <- max(row, na.rm = TRUE)
  
  # Define both relative and absolute tolerances
  rel_tol <- 1e-6
  abs_tol <- 1e-13
  
  is_not_na <- !is.na(row)
  
  # Relative comparison if max_val is large enough, otherwise use absolute
  if (abs(max_val) > abs_tol) {
    is_max <- abs(row - max_val) / abs(max_val) < rel_tol
  } else {
    is_max <- abs(row - max_val) < abs_tol
  }
  
  cells <- names(row)[which(is_not_na & is_max)]
  paste(cells, collapse = ";")
}), .SDcols = myCells]

print(head(unique(both$maxCell)))

### Order by max score.
setorder(both, -maxScore)

setnames(both, c("ensg", "gene"), c("Ensembl Gene ID", "Gene Name"))

### Make the final raw scores table. NO, need to redo with all scores for all genes!!

#myTraits_genes <- c("Ensembl Gene ID", "Gene Name", myCells)

#raw_scores_table <- both[, ..myTraits_genes]
#fwrite_headers(raw_scores_table, file = "./raw_scores_multiCOGS_table.txt")
#print(head(raw_scores_table))

### Make the table of prioritised genes.
#all_annot2 <- both[rowSums(all_annot[, ..myTraits] > 0.5, na.rm = TRUE) >0]

# set score to 0 if NA.
both[is.na(ILC3), ILC3 := 0]
both[is.na(CD4), CD4 := 0]

both[, (myCells) := lapply(.SD, function(x) {
  fifelse(is.na(x), NA_character_, # preserve NA, and make sure column is character type
          fifelse(x > 0.5, "+", "-")) # fifelse is a fast version of ifelse, from data.table
}), .SDcols = myCells]

# get the column order wanted.
myCells_genes_maxScore <- c("Ensembl Gene ID", "Gene Name", "maxScore", "maxCell", myCells)
s9 <- both[, ..myCells_genes_maxScore]
print(head(s9))

#setnames(summaryTable, c("IBD", "UC", "CD", "PSC", "Celiac"), c("Inflammatory Bowel Disease", "Ulcerative Colitis", "Crohn's Disease", "Primary Sclerosing cholangitis", "Celiac Disease"))

########## Now add again the genes with prior evidence.

## Evidence from OpenTargets
all_ot_to_bind <- all_ot[, .(L2G)]
all_ot_to_bind2 <- unique(as.data.table(separate_rows(all_ot_to_bind, L2G, sep = ", ")))
all_ot_to_bind2[, openTargets := "openTargets"]

s9_ot <- all_ot_to_bind2[s9, on = c(L2G = "Gene Name")]
setnames(s9_ot, "L2G", "gene")

## Evidence from IBDDB
ibd <- fread("./hILCs/previous_evidence/Inflammtory Bowel Disease.csv")
ibd2 <- unique(ibd[, .(`HGNC gene symbol`, "IBDDB")])

s9_ot_ibddb <- ibd2[s9_ot, on = c(`HGNC gene symbol` = "gene")]
setnames(s9_ot_ibddb, "HGNC gene symbol", "gene")

## Evidence from rioux et al (Ntunzwenimana 2021)
rioux <- data.table(gene = c("HNF4A", "IFIH1", "SMAD3", "SBNO2", "NFKB1", "NOD2", "ZFP36L1", 
                             "IRF1", "GIGYF1", "OTUD3", "AIRE", "PITX1", "KSR1", "DUSP16"), 
                    rioux = "Ntunzwenimana 2021")

s9_ot_ibddb_rioux <- rioux[s9_ot_ibddb, on = "gene"]

## Evidence from DisGeNET
dis_ibd_functional_bind <- unique(dis_ibd_functional[, .(Gene, Association_Type)])
dis_ibd_functional_bind2 <- as.data.table(dis_ibd_functional_bind %>% 
                                            group_by(Gene) %>% 
                                            summarise(Summary = paste(Association_Type, collapse = " & ")))

dis_ibd_functional_bind2[, Summary := paste("DisGeNET IBD", Summary, sep = ": ")]

dis_cd_functional_bind <- unique(dis_cd_functional[, .(Gene, Association_Type)])
dis_cd_functional_bind2 <- as.data.table(dis_cd_functional_bind %>% 
                                           group_by(Gene) %>% 
                                           summarise(Summary = paste(Association_Type, collapse = " & ")))

dis_cd_functional_bind2[, Summary := paste("DisGeNET CD", Summary, sep = ": ")]

disgenet <- rbind(dis_ibd_functional_bind2, dis_cd_functional_bind2)
disgenet2 <- as.data.table(disgenet %>% 
                             group_by(Gene) %>%
                             summarise(disgenet = paste(Summary, collapse = ", ")))

s9_ot_ibddb_rioux_disgenet <- disgenet2[s9_ot_ibddb_rioux, on = c(Gene = "gene")]

## Evidence from exome study

exome <- data.table(Gene = c("NOD2", "IL23R", "TYK2", "PTPN22", "CARD9", "GPR65", "MST1", "NOD2", "SMAD3", 
           "PTAFR", "ATG16L1", "SAG", "USP40", "HGFAC", "PDLIM5", "SLC39A8", "PTGER4", "IRGM", 
           "TAGAP", "DOK2", "SNAPC4", "SEC16A", "RELA", "IL10RA", "LRRK2", "SNX20", "CYLD", 
           "CCR7", "SDF2L1"), 
           Exome = "Sazonovs 2022 Exome")
s9_ot_ibddb_rioux_disgenet_exome <- exome[s9_ot_ibddb_rioux_disgenet, on = "Gene"]

#### Combine them together

s9_ot_ibddb_rioux_disgenet_exome[is.na(s9_ot_ibddb_rioux_disgenet_exome)] <- ""
s9_ot_ibddb_rioux_disgenet_exome[, `Prior evidence` := paste(openTargets, V2, rioux, disgenet, Exome, sep = ", ")]


s9_ot_ibddb_rioux_disgenet_exome[, c("openTargets", "V2", "rioux", "disgenet", "Exome") := NULL]

library(stringr)
s9_ot_ibddb_rioux_disgenet_exome$test <- str_replace_all(s9_ot_ibddb_rioux_disgenet_exome$`Prior evidence`, ", , , , ", "")
s9_ot_ibddb_rioux_disgenet_exome$test <- str_replace_all(s9_ot_ibddb_rioux_disgenet_exome$test, ", , , ", ", ")
s9_ot_ibddb_rioux_disgenet_exome$test <- str_replace_all(s9_ot_ibddb_rioux_disgenet_exome$test, ", , ", ", ")

s9_ot_ibddb_rioux_disgenet_exome$test <- trimws(s9_ot_ibddb_rioux_disgenet_exome$test) # Remove leading/trailing whitespace
s9_ot_ibddb_rioux_disgenet_exome$test <- sub(",+$", "", s9_ot_ibddb_rioux_disgenet_exome$test) # Remove trailing comma
s9_ot_ibddb_rioux_disgenet_exome$test <- sub("^,+", "", s9_ot_ibddb_rioux_disgenet_exome$test) # Remove initial comma

s9_ot_ibddb_rioux_disgenet_exome$`Prior evidence` <- s9_ot_ibddb_rioux_disgenet_exome$test
s9_ot_ibddb_rioux_disgenet_exome[, test := NULL]
setnames(s9_ot_ibddb_rioux_disgenet_exome, c("Gene", "maxScore", "maxCell"), c("Gene Name", "Max Gene Score", "Cell type(s) with Max Gene Score"))

######### write ###########
fwrite(s9_ot_ibddb_rioux_disgenet_exome, file = "./hILCs/COGS_results/COGS_out/Version3_revision2/ILC3_and_CD4_suppl_table_priorit_genes_with_prior_evidence.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)
######### ##### ###########
#setnames(summaryTable, c("maxScore", "maxTrait"), c("Max Gene Score", "Trait(s) with Max Gene Score"))

mynew <- s9_ot_ibddb_rioux_disgenet_exome[`Prior evidence` == ""]
myold <- s9_ot_ibddb_rioux_disgenet_exome[`Prior evidence` != ""]

######### Do a quick comparison with "ILC3 PIR specific" genes.
pirs <- fread("hILCs/COGS_results/comparisons/ILC_Javierre_PIR_analysis.txt")

pirs_new <- pirs[ensg %in% mynew$`Ensembl Gene ID` & `Prioritised based on Javierre PIRs` == "no" & BLOOD_Promoter == "no"]
pirs_old <- pirs[ensg %in% myold$`Ensembl Gene ID` & `Prioritised based on Javierre PIRs` == "no" & BLOOD_Promoter == "no"]






