library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
setwd("~/disease_annotation")
reli <- fread("~/RELI/20221019_ILC3_wABC_publication_RELI_results_summary.txt")

reli_sig <- reli[Enrichment >= 2 & BH_corrected_pval < 0.01 & Total >= 10]

### Need to:
### Sort by enrichment
### Re-label for supplementary table
setnames(reli_sig, c("disease", "Overlap", "Total", "STD"), c("Disease", "N_SNP_Overlap", "N_SNP_Total", "Std_Dev"))
reli_sig_write <- reli_sig[, .(Disease, `P-val`, BH_corrected_pval, N_SNP_Overlap, N_SNP_Total, Ratio, Mean, Std_Dev, `Z-score`)]
setorder(reli_sig_write, `P-val`)

fwrite(reli_sig_write, file = "./reli_sig_mergedPM.txt", sep = "\t", quote = F, row.names = F, col.names = T)

## Connect with accession no.s
access <- fread("./GWAS_1.0.2a_EUR_PMID_Accession.txt")
access2 <- as.data.table(access %>% separate_rows(`Study accession`, sep = ", "))
access2[, disease := str_replace_all(`Phenotype name`, "_", " ")] 
access2[, disease := str_replace_all(disease, "Sj\x9agren", "Sjögren")] # special case name!
reli[, disease := str_replace_all(disease, "Sj\x9agren", "Sjögren")] # special case name!

reli_access <- access2[reli, on = "disease", allow.cartesian = TRUE]
reli_access[is.na(PMIDs)]

# Make a table showing accession numbers
reli_access_small <- unique(reli_access[, c("disease", "Phenotype name", "PMIDs", "Study accession")])
reli_access_small_group <- as.data.table(reli_access_small %>% group_by(disease) %>% summarise(all_PMIDs = paste(PMIDs, collapse = ","), 
                                                                                 all_accession = paste(`Study accession`, collapse = ",")))

setnames(reli_access_small_group, c("disease", "all_PMIDs", "all_accession"), c("Disease", "PMIDs", "Study Accession Numbers"))
fwrite(reli_access_small_group, file = "./reli_access_small_group_phenos.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)

####### Are we significantly enriched for autoimmune diseases?
####### Using the term EFO0005140 as immune diseases. 
####### From the paper here https://www.sciencedirect.com/science/article/pii/S0092867421004293?via%3Dihub#sec4
####### "Dynamic landscape of immune cell-specific gene regulation in immune-mediated diseases"

immune <- fread("./efotraits_EFO_0005140-associations-2022-09-20.csv") # downloaded from GWAS catalog under EFO_0005140
immune[, type := "immune"]

# Now connect based on study accession number.

reli_im <- immune[reli_access, on = "Study accession", allow.cartesian = TRUE]
reli_im[is.na(type), type := "other"]
look <- unique(reli_im[, .(disease, `Phenotype name`, Overlap, Total, Enrichment, BH_corrected_pval, type)])

look2 <- as.data.table(look %>% group_by(disease, `Phenotype name`, Overlap, Total, Enrichment, BH_corrected_pval) %>% summarise(type = paste(type, collapse = ",")))

# No. of total immune diseases in RELI analysis
# Actually needs to be no. of total that had > 10 snps in the analysis, because these are the ones that we consider.

considered <- look2[Total >= 10]

immune <- considered[type %like% "immune"] 
no.immune.total <- length(unique(immune$disease))
other <- length(unique(considered$disease)) - no.immune.total
total <- length(unique(considered$disease))

# No. of significant signals
sig <- considered[Enrichment >= 2 & BH_corrected_pval < 0.01] 
no.sig <- length(unique(sig$disease))

sig_immune <- sig[type %like% "immune"]
no.sig.immune <- length(unique(sig_immune$disease))
sig_other <- length(unique(sig$disease)) - no.sig.immune #### up to here

# 1-tailed hypergeometric test
phyper(no.sig.immune-1, no.immune.total, total-no.immune.total, no.sig,lower.tail= FALSE)
