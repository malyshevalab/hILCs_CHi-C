library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
setwd("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/disease_annotation")
### Run in DT_DPLYR env
### Updating for revision, 14th April 2025

#reli <- fread("~/HRJ_monocytes/hILCs/RELI/20221019_ILC3_wABC_publication_RELI_results_summary.txt")
reli <- as.data.table(read_xlsx("~/HRJ_monocytes/hILCs/RELI/Disease_RELI_results_min_10.xlsx"))
reli <- reli[ATAC_library == "ILC3_chicago_fres_bin_5kb_abc_023_fres_extended_peakm_Jan25.txt_PIR_intersect_RE_hg19"]
#reli <- reli[ATAC_library == "CD4_chicago_fres_5kb_abc_023_fres_extended_peakm_Jan25.txt_PIR_intersect_RE_hg19"]

reli_sig <- reli[Padj < 0.05]

### Need to:
### Sort by enrichment
### Re-label according to Leah's suggestions (but now using phenotype, not disease, as many are traits)
setnames(reli_sig, c("Overlap", "Total"), c("N_SNP_Overlap", "N_SNP_Total"))
reli_sig_write <- reli_sig[, .(Phenotype, N_SNP_Overlap, N_SNP_Total, Ratio, Mean, Std_Dev, `Z_score`, Enrichment, P_val, Padj)]
setorder(reli_sig_write, `P_val`)

fwrite(reli_sig_write, file = "./reli_sig_mergedPM_2025_ILC.txt", sep = "\t", quote = F, row.names = F, col.names = T)
#fwrite(reli_sig_write, file = "./reli_sig_mergedPM_2025_CD4.txt", sep = "\t", quote = F, row.names = F, col.names = T)

## Connect with accession no.s
access <- fread("./GWAS_1.0.2a_EUR_PMID_Accession.txt")
access2 <- as.data.table(access %>% separate_rows(`Study accession`, sep = ", "))
access2[, disease := str_replace_all(`Phenotype name`, "_", " ")] # see if this connects
access2[, disease := str_replace_all(disease, "Sj\x9agren", "Sjögren")] # special case name!
reli[, Phenotype := str_replace_all(Phenotype, "Sj\x9agren", "Sjögren")] # special case name!

reli_access <- access2[reli, on = c(`Phenotype name` = "Phenotype"), allow.cartesian = TRUE]
reli_access[is.na(PMIDs)] # none


# Make a table showing accession numbers
reli_access_small <- unique(reli_access[, c("disease", "Phenotype name", "PMIDs", "Study accession")])
reli_access_small_group <- as.data.table(reli_access_small %>% group_by(disease) %>% summarise(all_PMIDs = paste(PMIDs, collapse = ","), 
                                                                                 all_accession = paste(`Study accession`, collapse = ",")))

setnames(reli_access_small_group, c("disease", "all_PMIDs", "all_accession"), c("Phenotype", "PMIDs", "Study Accession Numbers"))
fwrite(reli_access_small_group, file = "./reli_access_small_group_phenos_2025_ILC.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)
# We are actually looking at 495 traits, if we filter for > 10 loci. 

old <- fread("./reli_access_small_group_phenos.txt")
# 763 was original number of traits.
#library(ontologyIndex)

#ontology <- get_ontology("/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/disease_annotation/HumanDO.obo", 
#                         extract_tags="everything")
#ontology2 <- as.data.table(ontology)
#ontology3 <- ontology2[, .(id, name, is_a)]
#autoimmune <- ontology3[is_a %like% "DOID:417" | is_a == "DOID:2914" | is_a == "DOID:0060031" | is_a == "DOID:438" | is_a == "DOID:0050589"] 

#reli[, disease := tolower(disease)]
#ontology3[, name := tolower(name)]
#ontology3[, name2 := str_replace(name, "'", " ")]
#on_reli <- ontology3[reli, on = c(name2 = "disease")]
#on_reli_match <- on_reli[!is.na(is_a)]

#is_a <- fread("hILCs/disease_annotation/is_a", sep = "!", fill = TRUE, header = FALSE)
#is_a[, c("is_a", "DOID") := tstrsplit(V1, split = ": ")]
#is_a[, c("V1", "is_a") := NULL]
#setnames(is_a, "V2", "diseaseType")


#on_reli_match_w_type <- is_a[on_reli_match, on = c(DOID = "is_a")]
#unique(on_reli_match_w_type$diseaseType)

####### Now using the term EFO0005140 as immune diseases. We are signigicantly enriched for autoimmune diseases.
####### From the paper here https://www.sciencedirect.com/science/article/pii/S0092867421004293?via%3Dihub#sec4
####### "Dynamic landscape of immune cell-specific gene regulation in immune-mediated diseases"

####### Using "immune system disease" EFO:0000540

immune <- fread("./efotraits_EFO_0005140-associations-2022-09-20.csv")
immune[, type := "immune"]
# can also use this file: "gwas-association-downloaded_2022-09-20-EFO_1000870-withChildTraits.tsv" Both were downloaded from GWAS catalog under EFO_0005140

# Now connect based on study accession number.

reli_im <- immune[reli_access, on = "Study accession", allow.cartesian = TRUE]
reli_im[is.na(type), type := "other"]
look <- unique(reli_im[, .(disease, `Phenotype name`, Overlap, Total, Enrichment, Padj, type)])

look2 <- as.data.table(look %>% group_by(disease, `Phenotype name`, Overlap, Total, Enrichment, Padj) %>% summarise(type = paste(type, collapse = ",")))

# No. of total immune diseases in RELI analysis
immune <- look2[type %like% "immune"] 
no.immune.total <- length(unique(immune$disease))
other <- length(unique(look2$disease)) - no.immune.total
total <- length(unique(look2$disease))

# No. of significant signals
sig <- look2[Padj < 0.05] 
no.sig <- length(unique(sig$disease))

sig_immune <- sig[type %like% "immune"]
no.sig.immune <- length(unique(sig_immune$disease))
sig_other <- length(unique(sig$disease)) - no.sig.immune 

phyper(no.sig.immune-1, no.immune.total, total-no.immune.total, no.sig,lower.tail= FALSE)
# P = 1.076882e-05 1-tailed hypergeometric test
