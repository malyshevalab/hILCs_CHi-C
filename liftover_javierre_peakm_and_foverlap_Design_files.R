### Liftover Javierre data to GRCh38. Using previously prepared design files in hg38.

peakm <- fread("~/Javierre/PCHiC_peak_matrix_cutoff5.txt")

### Make bait and OE. Each one has an ID for the interaction pair. 
peakm[, intID := paste(baitID, oeID, sep = "_")]
baits <- peakm[, .(baitChr, baitStart, baitEnd, intID)]
peakm[which(duplicated(peakm$intID))] # none. 

fwrite(baits, file = "baits_peakm5.bed", sep = "\t", quote = F, row.names = F, col.names = F)

oe <- peakm[, .(oeChr, oeStart, oeEnd, intID)]
fwrite(oe, file = "oe_peakm5.bed", sep = "\t", quote = F, row.names = F, col.names = F)

### Now doing liftover on the command line
#~/bin/liftOver -bedPlus=3 baits_peakm5.bed ~/bin/liftover/hg19ToHg38_nochr.over.chain baits_peakm5_hg38.bed baits_peakm5_liftover_unmapped.bed
#~/bin/liftOver -bedPlus=3 oe_peakm5.bed ~/bin/liftover/hg19ToHg38_nochr.over.chain oe_peakm5_hg38.bed oe_peakm5_liftover_unmapped.bed

baits_hg38 <- fread("./baits_peakm5_hg38.bed")
names(baits_hg38) = c("baitChr38", "baitStart38", "baitEnd38", "intID")
oe_hg38 <- fread("./oe_peakm5_hg38.bed")
names(oe_hg38) = c("oeChr38", "oeStart38", "oeEnd38", "intID")

peakm_no_coords <- copy(peakm)
peakm_no_coords[, c("baitChr", "baitStart", "baitEnd", "oeChr", "oeStart", "oeEnd") := NULL]
baits_peakm_hg38 <- baits_hg38[peakm_no_coords, on = "intID"]
baits_oe_peakm_hg38 <- oe_hg38[baits_peakm_hg38, on = "intID"]
baits_oe_peakm_hg38[which(duplicated(baits_oe_peakm_hg38$intID))] 

baits_oe_peakm_hg38[is.na(baitChr38)]
baits_oe_peakm_hg38[is.na(oeChr38)]

baits_oe_peakm_hg38_mapped <- baits_oe_peakm_hg38[!is.na(baitChr38) & !is.na(oeChr38)]

final <- length(unique(baits_oe_peakm_hg38_mapped$intID))
started <- length(unique(peakm$intID))

started - final # we lose 666 interactions out of ~700K

baits_oe_peakm_hg38_mapped_write <- baits_oe_peakm_hg38_mapped[, .(baitChr38, baitStart38, baitEnd38, baitID, baitName, oeChr38, 
                                                                   oeStart38, oeEnd38, oeID, oeName, dist, Mon, Mac0, Mac1, Mac2, 
                                                                   Neu, MK, EP, Ery, FoeT, nCD4, tCD4, aCD4, naCD4, nCD8, tCD8, nB, 
                                                                   tB, clusterID, clusterPostProb)]
setnames(baits_oe_peakm_hg38_mapped_write, c("baitChr38", "baitStart38", "baitEnd38", "oeChr38", "oeStart38", "oeEnd38"), 
         c("baitChr", "baitStart", "baitEnd", "oeChr", "oeStart", "oeEnd"))


fwrite(baits_oe_peakm_hg38_mapped_write, file = "./PCHiC_peak_matrix_cutoff5_hg38Liftover.txt", sep = "\t", quote = F, row.names = F, col.names = T)


#### Intersect with previously lifted over rmap and baitmap.
jav <- "PCHiC_peak_matrix_cutoff5_hg38Liftover.txt"
rmap <- "Human_HindIII_hg38.rmap"
baitmap <- "Human_HindIII_hg38.baitmap"

jav_pm <- fread(jav)
rmap_dt <- fread(rmap)

# get the correct coords and bait IDs and OE IDs from this rmap.
# check they match the baitmap
baitmap_dt <- fread(baitmap)
baitmap_dt

setkey(jav_pm, baitChr, baitStart, baitEnd)
jav_pm_fo1 <- foverlaps(rmap_dt, jav_pm, by.x = c("V1", "V2", "V3"), nomatch = NULL)

# Look for duplicates.
jav_pm_fo1_check <- unique(jav_pm_fo1[, .(V1, baitStart, baitEnd, baitID, V2, V3, V4)])
jav_pm_fo1_check[which(duplicated(jav_pm_fo1_check[, .(V2)]))]
jav_pm_fo1_check[V2 == 145709804] # overlap by 2bp
jav_pm_fo1_check[V2 == 145921318] # overlap by 2bp
jav_pm_fo1_check[V2 == 47382881] # overlap by 2bp

jav_pm_fo1_check[which(duplicated(jav_pm_fo1_check[, .(baitStart)]))]
jav_pm_fo1_check[baitStart == 145705060]

# Require a minimum overlap of 20 bp and allow them to be split in new.
# Need to use a different package - genomicRanges
library(GenomicRanges)

jav_pm[, int_id := paste(baitID, oeID, sep = "_")]
jav_pm_small <- jav_pm[, .(baitChr, baitStart, baitEnd, baitID, oeChr, oeStart, oeEnd, oeID, int_id)]

names(rmap_dt) = c("Chr", "hg38Start", "hg38End", "hg38ID")

mybaits <- makeGRangesFromDataFrame(jav_pm_small, keep.extra.columns = T, ignore.strand = T, 
                                    seqnames.field = "baitChr", start.field = "baitStart", 
                                    end.field = "baitEnd")

newfrags <- makeGRangesFromDataFrame(rmap_dt, keep.extra.columns = T, ignore.strand = T, 
                                     seqnames.field = "Chr", start.field = "hg38Start", end.field = "hg38End")

newbaits <- as.data.table(mergeByOverlaps(mybaits, newfrags, minoverlap = 20, ignore.strand = T))

newbaits_small <- newbaits[, .(newfrags.seqnames, newfrags.start, newfrags.end, newfrags.hg38ID, 
                               mybaits.oeChr, mybaits.oeStart, mybaits.oeEnd, mybaits.oeID, int_id)]

names(newbaits_small) = c("baitChr", "baitStart", "baitEnd", "baitID", "hg19_oeChr", "hg19_oeStart", 
                          "hg19_oeEnd", "hg19_oeID", "int_id")


myoe <- makeGRangesFromDataFrame(newbaits_small, keep.extra.columns = T, ignore.strand = T, 
                                    seqnames.field = "hg19_oeChr", start.field = "hg19_oeStart", 
                                    end.field = "hg19_oeEnd")

newoe <- as.data.table(mergeByOverlaps(myoe, newfrags, minoverlap = 20, ignore.strand = T))

newoe_small <- newoe[, .(myoe.baitChr, myoe.baitStart, myoe.baitEnd, myoe.baitID, 
                         newfrags.seqnames, newfrags.start, newfrags.end, newfrags.hg38ID, int_id)]

names(newoe_small) = c("baitChr", "baitStart", "baitEnd", "baitID", "oeChr", 
                       "oeStart", "oeEnd", "oeID", "int_id")

### Now join back to the original peakmatrix, based on int_id (the old bait_oe)

jav_pm2 <- copy(jav_pm)
jav_pm2[, c("baitChr", "baitStart", "baitEnd", "baitID", "oeChr", "oeStart", "oeEnd", "oeID") := NULL]

final_pm <- newoe_small[jav_pm2, on = "int_id", nomatch = NULL]

final_pm_write <- final_pm[, .(baitChr, baitStart, baitEnd, baitID, baitName, oeChr, oeStart, oeEnd, oeID, oeName, 
                               dist, Mon, Mac0, Mac1, Mac2, Neu, MK, EP, Ery, FoeT, nCD4, tCD4, aCD4, naCD4, nCD8, 
                               tCD8, nB, tB, clusterID, clusterPostProb)]

fwrite(final_pm_write, file = "PCHiC_peak_matrix_cutoff5_hg38Liftover_foverlapped.txt", sep = "\t", 
       quote = F, row.names = F, col.names = T)

###