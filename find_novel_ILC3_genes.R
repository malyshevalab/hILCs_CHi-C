library(data.table)
library(tidyr)
library(dplyr)
setwd("~")

#### Get the annotation table.
s9 <- fread("./COGS_results/COGS_out/CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/Annotated_COGS_scores_data.table.txt")

## Evidence from OpenTargets
ot1 <- fread("./previous_evidence/GCST004132-independently-associated-loci.tsv")
ot2 <- fread("./previous_evidence/GCST003044-independently-associated-loci.tsv")
ot3 <- fread("./previous_evidence/GCST000879-independently-associated-loci.tsv")
ot4 <- fread("./previous_evidence/GCST000207-independently-associated-loci.tsv")
ot5 <- fread("./previous_evidence/GCST005537-independently-associated-loci.tsv")
all_ot <- rbind(ot1, ot2, ot3, ot4, ot5)

all_ot_to_bind <- all_ot[, .(L2G)]
all_ot_to_bind2 <- unique(as.data.table(separate_rows(all_ot_to_bind, L2G, sep = ", ")))
all_ot_to_bind2[, openTargets := "openTargets"]

s9_ot <- all_ot_to_bind2[s9, on = c(L2G = "gene")]
setnames(s9_ot, "L2G", "gene")

## Evidence from IBDDB
ibd <- fread("./previous_evidence/Inflammtory Bowel Disease.csv")
ibd2 <- unique(ibd[, .(`HGNC gene symbol`, "IBDDB")])

s9_ot_ibddb <- ibd2[s9_ot, on = c(`HGNC gene symbol` = "gene")]
setnames(s9_ot_ibddb, "HGNC gene symbol", "gene")

## Evidence from rioux et al (Ntunzwenimana 2021)
rioux <- data.table(gene = c("HNF4A", "IFIH1", "SMAD3", "SBNO2", "NFKB1", "NOD2", "ZFP36L1", 
                             "IRF1", "GIGYF1", "OTUD3", "AIRE", "PITX1", "KSR1", "DUSP16"), 
                    rioux = "Ntunzwenimana 2021")

s9_ot_ibddb_rioux <- rioux[s9_ot_ibddb, on = "gene"]

## Evidence from DisGeNET
## from here https://www.disgenet.org/browser/0/1/0/C0021390/ 
## Either looking at CD, IBD or both. 

dis_ibd <- fread("./previous_evidence/IBD_C0021390_disease_gda_evidences.txt")
dis_ibd_functional <- dis_ibd[Association_Type != "GeneticVariation"]
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

s9_ot_ibddb_rioux_disgenet_exome_sig <- s9_ot_ibddb_rioux_disgenet_exome[cogs >= 0.5]
s9_ot_ibddb_rioux_disgenet_exome_sig[is.na(s9_ot_ibddb_rioux_disgenet_exome_sig)] <- ""
s9_ot_ibddb_rioux_disgenet_exome_sig[, `Prior evidence` := paste(openTargets, V2, rioux, disgenet, Exome, sep = ", ")]


s9_ot_ibddb_rioux_disgenet_exome_sig[, c("openTargets", "V2", "rioux", "disgenet", "Exome") := NULL]


setnames(s9_ot_ibddb_rioux_disgenet_exome_sig, c("Gene", "cogs"), c("gene", "ILC3_cogs"))

s9_ot_ibddb_rioux_disgenet_exome_sig$test <- str_replace_all(s9_ot_ibddb_rioux_disgenet_exome_sig$`Prior evidence`, ", , , , ", "")
s9_ot_ibddb_rioux_disgenet_exome_sig$test <- str_replace_all(s9_ot_ibddb_rioux_disgenet_exome_sig$test, ", , , ", ", ")
s9_ot_ibddb_rioux_disgenet_exome_sig$test <- str_replace_all(s9_ot_ibddb_rioux_disgenet_exome_sig$test, ", , ", ", ")

s9_ot_ibddb_rioux_disgenet_exome_sig$test <- trimws(s9_ot_ibddb_rioux_disgenet_exome_sig$test) # Remove leading/trailing whitespace
s9_ot_ibddb_rioux_disgenet_exome_sig$test <- sub(",+$", "", s9_ot_ibddb_rioux_disgenet_exome_sig$test) # Remove trailing comma
s9_ot_ibddb_rioux_disgenet_exome_sig$test <- sub("^,+", "", s9_ot_ibddb_rioux_disgenet_exome_sig$test) # Remove initial comma

s9_ot_ibddb_rioux_disgenet_exome_sig$`Prior evidence` <- s9_ot_ibddb_rioux_disgenet_exome_sig$test
s9_ot_ibddb_rioux_disgenet_exome_sig[, test := NULL]

######### write ###########
fwrite(s9_ot_ibddb_rioux_disgenet_exome_sig, file = "./COGS_results/COGS_out/CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/Annotated_COGS_scores_data.table_with_prior_info_SuSIE_combinedInteractions.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)
######### ##### ###########





