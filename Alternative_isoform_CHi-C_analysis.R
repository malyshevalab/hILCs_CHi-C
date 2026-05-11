library(Chicago)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(RColorBrewer)

######################################################
bin_dir <- "~/hILCs_all_merged_bin_Step2/data"   #####
des_dir <- "~/Design/Human_hg38_bin5K_sol_baits" #####
######################################################

#########################################################################
######## General properties of per-gene interactions in the PCHiC data
#########################################################################

bin <- readRDS(paste(bin_dir, "hILCs_all_merged_bin_Step2.Rds", sep = "/"))
baits <- fread(paste(des_dir, "human_DpnII_5K_sol_baits.baitmap", sep = "/"))

names(baits) = c("Chr", "baitStart", "baitEnd", "baitID", "Gene")

### Split baits into individual genes
baits_names <- baits[, list(Gene = unlist(strsplit(Gene, "/"))), by=c("Chr", "baitStart", "baitEnd", "baitID")]

### Restrict to a nominal score 3.
int <- bin@x[score >= 3] 

### Intersect interactions with baitIDs (bait-to-bait are listed twice)
setkey(baits_names, baitID)
bint <- int[baits_names, on = "baitID", nomatch = NULL, allow.cartesian = TRUE]
bintG <- bint[Gene != "off_target"]

### Save bintG object
fwrite(bintG, file = "ILC3_all_Chicago_score3_with_baits_per_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### How many TSS per gene?
TSS <- as.data.table(bintG %>% distinct(Gene, baitID) %>% group_by(Gene) %>% tally)

### Get the median, and the number of genes that had more than one TSS in the data.
median(TSS$n)
length(unique(TSS$Gene))
more_than_one <- TSS[n > 1]
length(unique(more_than_one$Gene))

# Get the mode.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

mynumbers <- more_than_one$n
check <- getmode(mynumbers)
print(check)

#### How many interactions per TSS or per gene? (need to include TSS with no interactions, though!)
sbintG <- bintG[score >= 5 & N >= 5]

intTSS <- as.data.table(sbintG %>% distinct(Gene, baitID, otherEndID) %>% group_by(Gene, baitID) %>% tally)
noint1 <- bintG[!baitID %in% sbintG$baitID]
noint1[, n := 0]
nointTSS <- unique(noint1[, .(Gene, baitID, n)])
all_intTSS <- rbind(intTSS, nointTSS)

# check the average no. interactions per TSS
check <- unique(intTSS[, .(baitID, n)])
median(check$n)

pdf(file = "Interactions_per_BaitID_in_CHiC_data.pdf")
hist(all_intTSS$n, breaks = 100, main = "Number of significant interactions per baitID - including zeroes")
dev.off()

all_intTSS_byGene <- as.data.table(all_intTSS %>% group_by(Gene) %>% dplyr::summarize(total_interactions = sum(n)))

pdf(file = "Interactions_per_Gene_in_CHiC_data.pdf")
hist(all_intTSS_byGene$total_interactions, breaks = 100, main = "Number of significant interactions per Gene - including zeroes")
dev.off()

###########################################################################
######## Looking at alternative interactions between different promoters.
###########################################################################

### Table of significant interactions
### Count the no. of baited promoters
### Group by gene and otherEndID.
### Compare this no. with the total no. of baited promoters.

# 1. Only include baited promoters that have at least one interaction. At score 5 and N>=5
# For these promoters, make sure that the associated gene has more than one TSS. Join to genes at score5.

genes_baits_sig <- copy(sbintG)
genes_baits_sig[, gene_bait := paste(Gene, baitID, sep = "_")] # make an ID for gene_bait
bintG_multiTSS1 <- as.data.table(genes_baits_sig %>% distinct(Gene, baitID) %>% group_by(Gene) %>% tally())
multi <- bintG_multiTSS1[n>1]

allScores_genes_baits <- genes_baits_sig[Gene %in% multi$Gene]
length(unique(allScores_genes_baits$Gene))
# Now 1661 genes.

# To plot no. TSS we actually need to use int without restricting to score3
TSSCopy <- copy(bintG_multiTSS1)
TSSCopy[, n := as.factor(n)]
pdf(file = "Baited_interacting_TSS_per_gene_in_CHiC_data_N5.pdf")
p <- ggplot(TSSCopy, aes(x = n))
p + geom_bar() + ylab("Number of genes") + xlab("Number of targeted promoters") + theme(text = element_text(size = 20))
dev.off()

# 2. Add the col for the no. of baited promoters that had interactions.
allScores_gbp <- bintG_multiTSS1[allScores_genes_baits, on = "Gene"]
setnames(allScores_gbp, "n", "No.promoters")
# in this file, scores are 5.
# make a file with scores > 3. Will query this dataset to look for interactions with other promoters.
score3 <- bintG[score >=3] # Here we will not restrict on N.
# query this file for "shared" interactions.
score5_dt <- copy(allScores_gbp)
score3_dt <- copy(score3)
fwrite(score5_dt, file = "./score5_dt", sep = "\t", quote = F, row.names = F, col.names = T)
fwrite(score3_dt, file = "./score3_dt", sep = "\t", quote = F, row.names = F, col.names = T)
rm(bin)


# 3. Count how many promoters (baits) each PIR interacts with for the same gene.
# in the score3 table, we can only have promoters and other ends who are also in score 5 table (but interacting in different combinations). 
score3_keep2 <- score3_dt[baitID %in% score5_dt$baitID]
# We can also have other ends + / -1
score5_dt[, otherEndID_plusone := otherEndID + 1]
score5_dt[, otherEndID_minusone := otherEndID - 1]
score3_keep1 <- score3_keep2[otherEndID %in% score5_dt$otherEndID | otherEndID %in% score5_dt$otherEndID_plusone | otherEndID %in% score5_dt$otherEndID_minusone]
# Gene must also be in score5
score3_keep <- score3_keep1[Gene %in% score5_dt$Gene]
score3_keep[, otherEndID_plusone := otherEndID + 1]
score3_keep[, otherEndID_minusone := otherEndID - 1]
# Have to ask, for each PIR, does it (or one plus/minus one) interact with which promoters?

score5_small <- unique(score5_dt[, .(Gene, otherEndID)])
names(score5_small)[2] = "score5_otherEndID"
score3_small <- unique(score3_keep[, .(Gene, baitID, otherEndID, otherEndID_plusone, otherEndID_minusone, score)])

# Look at each PIR in score5 data.table. (called them score5_otherEndID)
# Combine with score3_small, based on Gene
# Only keep the ones where DpnID or DpnID_plusone or DpnID_minusone == score5_DpnID (replace all below pseudo_otherEndID)
# Then count promoters.

interaction1 <- score3_small[score5_small, on = c("Gene"), allow.cartesian = TRUE]
interaction2 <- interaction1[score5_otherEndID == otherEndID | score5_otherEndID == otherEndID_plusone | score5_otherEndID == otherEndID_minusone]

dt_small <- interaction2[, .(Gene, baitID, score5_otherEndID)]

res <- as.data.table(dt_small %>% distinct %>% group_by(Gene, score5_otherEndID) %>% distinct %>% tally())
setnames(res, "n", "No.interactingPromoters")

dt2 <- res[score5_dt, on = c("Gene", score5_otherEndID = "otherEndID"), nomatch = NULL]
# Now ask how the total no. promoters and the interacting promoters compare.
dt2[No.interactingPromoters == 1, PIR_Category := "Distinct"]
dt2[No.interactingPromoters < No.promoters & No.interactingPromoters > 1, PIR_Category := "Partially_shared"]
dt2[No.interactingPromoters == No.promoters, PIR_Category := "Fully_shared"]
dt2[No.interactingPromoters > No.promoters, PIR_Category := "Problem"]

categories <- copy(dt2)


fwrite(categories, file = "./ILC3_designated_otherEnds_isoforms_withDists_threeCategories_binned_corrected_N5.txt", 
       sep = "\t", quote = F, row.names = F, col.names = T)

###############################################################
################## Look at individual genes, make plots
###############################################################

### Look at INPP4B
inp <- categories[Gene == "INPP4B"]
inpp <- unique(inp[, . (Gene, baitID, score5_otherEndID, otherEndID_plusone, otherEndID_minusone, PIR_Category)])
setorder(inpp, baitID, score5_otherEndID)

### Plot 1. The proportion of all PIRs that are in each category
cats <- unique(categories[, .(score5_otherEndID, PIR_Category)])
cats2 <- cats[, .N, by=.(PIR_Category)]
display.brewer.all()
colourcount = length(unique(cats$PIR_Category))
getPalette = colorRampPalette(brewer.pal(6, "Set2"))
cats2$PIR_category <- factor(cats2$PIR_Category, levels = c("Fully_shared", "Partially_shared", "Distinct"))
bp<- ggplot(cats2, aes(x="", y=N, fill=PIR_category)) + scale_fill_manual(values=getPalette(colourcount))
pie <- bp + geom_col() +coord_polar("y", start=0)
pdf(file = "./PIR_categories_pie_binned_N5.pdf")
pie + theme(text = element_text(size = 15))
dev.off()


### Plot 2. A comparison of the no.s of PIRs in each category, per gene.
# Get the category back when the count is zero.
to_plot1 <- unique(categories[, .(Gene, score5_otherEndID, PIR_Category)]) 
to_plot <- as.data.table(to_plot1 %>% group_by(Gene, PIR_Category) %>% tally())

mygenes <- unique(to_plot$Gene)
mygenes2 <- rep(mygenes, 3)
mygenes3 <- sort(mygenes2)
mylength <- length(unique(mygenes3))

mycats <- rep(c("Distinct", "Partially_shared", "Fully_shared"), mylength)

dt_cats <- data.table(Gene = mygenes3, PIR_Category = mycats)

to_plot_allCats <- to_plot[dt_cats, on = c("Gene", "PIR_Category")]
to_plot_allCats[is.na(n), n := 0] 
to_plot_allCats$PIR_category <- factor(to_plot_allCats$PIR_Category, levels = c("Fully_shared", "Partially_shared", "Distinct"))

l <- ggplot(to_plot_allCats, aes(x = PIR_category, y = n, fill = PIR_category)) +
  scale_fill_manual(values=getPalette(colourcount))
pdf(file = "./no_of_PIRs_per_cat_per_gene_binned_N5.pdf")
l + geom_boxplot() + ylab("Number of PIRs per gene") + xlab("PIR category") +
  theme(text = element_text(size = 20), legend.position = "none")
dev.off()


######### Ask the following, from the perspective of PIRs:######### 
### 1. What is the mean distance of the PIR to all promoters of the given gene, irrespective of interactions being present?
# Note, only include alternative promoters that have at least one PIR
### 2. Plot the distribution of these distances for shared/non-shared PIRs
# What is a shared PIR? Could either - 
# Interact with >1 PIR for a given gene
# Interact with all PIRs for a given gene

### Add the DpnII locations back to our results.
rmap <- fread(paste0(des_dir, "/human_DpnII_5K_sol_baits.rmap"))
names(rmap) = c("Chr", "DpnStart", "DpnEnd", "DpnID")
rmap[, DpnMid := round((DpnStart + DpnEnd)/2)]

all_Results <- unique(categories[, .(Gene, baitID, score5_otherEndID, No.interactingPromoters, No.promoters, score, 
                                     PIR_Category)])
ar1 <- rmap[all_Results, on = c(DpnID = "baitID")]
setnames(ar1, c("DpnStart", "DpnEnd", "DpnID", "Chr", "DpnMid"), c("bait_DpnStart", "bait_DpnEnd", "bait_DpnID", "baitChr", "bait_DpnMid"))
ar2 <- rmap[ar1, on = c(DpnID = "score5_otherEndID")]
setnames(ar2, c("DpnStart", "DpnEnd", "DpnID", "Chr", "DpnMid"), c("OE_DpnStart", "OE_DpnEnd", "OE_DpnID", "OEChr", "OE_DpnMid"))
ar2[OEChr != baitChr] # remove trans! (there is now only one)
ar <- ar2[baitChr == OEChr]

### Need to calculate the distances to all promoters, regardless of interaction. So, make a table with all combinations.
### Essentially, make a new table of gene-TSS, and join with this one, making new rows. then calc average.
ba <- unique(ar[, .(Gene, bait_DpnID, bait_DpnMid)])
oe <- unique(ar[, .(Gene, OE_DpnID, OE_DpnMid)])
ba_oe <- ba[oe, on = "Gene", nomatch = NULL, allow.cartesian = TRUE]
### In this table, OE_DpnID is the "PIR"
ba_oe[, PIR_TSSdist := abs(OE_DpnMid - bait_DpnMid)]
hist(ba_oe$PIR_TSSdist)

ba_oe_dis <- as.data.table(ba_oe %>% group_by(OE_DpnID) %>% mutate(av_PIR_TSSdist = round(mean(PIR_TSSdist))))

ba_oe_dis_2merge <- unique(ba_oe_dis[, .(Gene, OE_DpnID, av_PIR_TSSdist)])

ar_avDist1 <- ba_oe_dis_2merge[ar, on = c("Gene", "OE_DpnID")]
ar_avDist <- ar_avDist1[!is.na(av_PIR_TSSdist)] 

fwrite(ar_avDist, file = "./ILC3_designated_otherEnds_isoforms_withDists_threeCategories_binned_N5.txt", sep = "\t", quote = F, 
       row.names = F, col.names = T)


### Plot the distribution of these distances for shared/non-shared PIRs, as they have already been designated.
getPalette = colorRampPalette(brewer.pal(6, "Set2"))
ar_avDist <- fread("./ILC3_designated_otherEnds_isoforms_withDists_threeCategories_binned_N5.txt")
colourcount = length(unique(ar_avDist$PIR_Category))
to_plot <- unique(ar_avDist[, .(OE_DpnID, av_PIR_TSSdist, PIR_Category)])
to_plot$PIR_category <- factor(to_plot$PIR_Category, levels = c("Fully_shared", "Partially_shared", "Distinct"))

pdf(file = "./av_PIR_TSS_dist_in_threeCategories_binned_N5.pdf")
p <- ggplot(to_plot, aes(x = PIR_category, y = av_PIR_TSSdist, fill = PIR_category)) +
  scale_fill_manual(values=getPalette(colourcount))
p + geom_boxplot() + xlab("PIR Category") + ylab("Average PIR-promoter distance") +
  theme(text = element_text(size = 20))
dev.off()

# Look at the log
to_plot[, log_av_PIR_TSSdist := log(av_PIR_TSSdist)]

pdf(file = "./log_av_PIR_TSS_dist_in_categories_binned_N5.pdf")
q <- ggplot(to_plot, aes(x = PIR_category, y = log_av_PIR_TSSdist, fill = PIR_category)) + 
  scale_fill_manual(values=getPalette(colourcount))
q + geom_boxplot() + xlab("PIR Category") + ylab("Log average PIR-promoter distance") +
  theme(text = element_text(size = 20))
dev.off()

fully_shared <- to_plot[PIR_category == "Fully_shared", av_PIR_TSSdist]
partially_shared <- to_plot[PIR_category == "Partially_shared", av_PIR_TSSdist]
distinct <-  to_plot[PIR_category == "Distinct", av_PIR_TSSdist]

print(median(fully_shared)) 
print(median(partially_shared)) 
print(median(distinct)) 

kruskal.test(av_PIR_TSSdist ~ PIR_category, data = to_plot)
pairwise.wilcox.test(to_plot$av_PIR_TSSdist, to_plot$PIR_category, p.adj="bonferroni", exact = F)
# groups are non identical
# each value represents a different "PIR-Gene" pairing.

##########################################################
################ Looking at distance cutoffs
##########################################################

# At a cutoff of 1 E06

mydt <- copy(to_plot)

######### Plotting function 1
get_plot <- function(dt, maxdist) {
  to_plot <- dt[av_PIR_TSSdist <= maxdist]
fully_shared <- to_plot[PIR_category == "Fully_shared", av_PIR_TSSdist]
partially_shared <- to_plot[PIR_category == "Partially_shared", av_PIR_TSSdist]
distinct <-  to_plot[PIR_category == "Distinct", av_PIR_TSSdist]

print(ggplot() +
  geom_density(aes(fully_shared, color = "Fully Shared")) +
  geom_density(aes(partially_shared, color = "Partially Shared")) +
  geom_density(aes(distinct, color = "Distinct")) +
    xlab("Average PIR-TSS distance, bp") +
    ylab("PIR density per category") + 
    theme(text = element_text(size = 15)) + 
  scale_color_manual(name = "PIR Category", values=c(`Fully Shared` = "#33CC99", `Partially Shared` = "#CC66FF", `Distinct` = "#FFCC33")))
}
######### 

pdf(file = "./PIR_cat_density_plot_1m.pdf")
get_plot(mydt, 1000000)
dev.off()

get_plot(mydt, 500000)

png(file = "./PIR_cat_density_plot_150K.png")
get_plot(mydt, 150000)
dev.off()
png(file = "./PIR_cat_density_plot_100K.png")
get_plot(mydt, 100000)
dev.off()
png(file = "./PIR_cat_density_plot_50K.png")
get_plot(mydt, 50000) # at 15kb, there are more "distinct" interactions than other categories.
dev.off()
get_plot(mydt, 25000) 

# At which distance are the distance effects eliminated?

######### plotting function 2
get_median_dist <- function(dt, distance) {
  to_plot <- dt[av_PIR_TSSdist <= distance]
  fully_shared <- to_plot[PIR_category == "Fully_shared", av_PIR_TSSdist]
  partially_shared <- to_plot[PIR_category == "Partially_shared", av_PIR_TSSdist]
  distinct <-  to_plot[PIR_category == "Distinct", av_PIR_TSSdist]
  
  print(paste0("The median distance for fully shared PIRs under ", distance, "bp is ", median(fully_shared)))
  print(paste0("The median distance for partially shared PIRs under ", distance, "bp is ", median(partially_shared))) 
  print(paste0("The median distance for distinct PIRs under ", distance, "bp is ", median(distinct)))
  
  p <- ggplot(to_plot, aes(x = PIR_category, y = av_PIR_TSSdist, fill = PIR_category)) +
    scale_fill_manual(values=getPalette(colourcount))
  print(p + geom_boxplot() + xlab("PIR Category") + ylab("Average PIR-promoter distance") +
    theme(text = element_text(size = 20)))
  
  print(kruskal.test(av_PIR_TSSdist ~ PIR_category, data = to_plot))
  print(pairwise.wilcox.test(to_plot$av_PIR_TSSdist, to_plot$PIR_category, p.adj="bonferroni", exact = F))
}
######### 

get_median_dist(mydt, 1000000)
get_median_dist(mydt, 150000)
get_median_dist(mydt, 100000)
get_median_dist(mydt, 50000)
get_median_dist(mydt, 25000)
get_median_dist(mydt, 10000)
get_median_dist(mydt, 9700)

png(file = "./PIRs_dist_numbers_per_cat_50kb.png")
under50 <- mydt[av_PIR_TSSdist <= 50000]
p <- ggplot(under50, aes(x = PIR_category, y = av_PIR_TSSdist, fill = PIR_category)) +
  scale_fill_manual(values=getPalette(colourcount))
print(p + geom_boxplot() + xlab("PIR Category") + ylab("Average PIR-promoter distance") +
        theme(text = element_text(size = 15)))
dev.off()

# what is an eg of a gene with PIRs in 50kb dist range
in_50 <- ar_avDist[av_PIR_TSSdist <= 50000]
# remove genes in other dist ranges
not_50 <- ar_avDist[av_PIR_TSSdist > 50000]
in_50_only <- in_50[!Gene %in% not_50$Gene]

### Get a pie chart for these.
for_pie <- in_50_only[, .(OE_DpnID, PIR_Category)]
for_pie_unique <- unique(for_pie)
for_pie2 <- for_pie_unique[, .N, by=.(PIR_Category)]


colourcount = length(unique(for_pie2$PIR_Category))
getPalette = colorRampPalette(brewer.pal(6, "Set2"))
for_pie2$PIR_category <- factor(for_pie2$PIR_Category, levels = c("Fully_shared", "Partially_shared", "Distinct"))
bp<- ggplot(for_pie2, aes(x="", y=N, fill=PIR_category)) + scale_fill_manual(values=getPalette(colourcount))
pie <- bp + geom_col() +coord_polar("y", start=0)
pdf(file = "./PIR_categories_pie_binned_N5_50kb.pdf")
pie + theme(text = element_text(size = 15))
dev.off()

########### 

