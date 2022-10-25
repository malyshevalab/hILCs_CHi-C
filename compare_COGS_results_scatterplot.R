######## Compare the COGS results across multiple datasets

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

mydir <- "~/COGS_results/"
results <- paste0(mydir, "COGS_out")
outdir = paste0(mydir, "comparisons")
dir.create(outdir)
list.files(results)

##### To run with pre-annotated results
make_comparison_plot_annotated <- function(res1, name1, res2, name2, outdir, outname, imageType) {
  results1 <- fread(paste(results, res1, sep = "/"))
  results2 <- fread(paste(results, res2, sep = "/"))
  results1[, chr37 := as.character(chr37)]
  results2[, chr37 := as.character(chr37)]
  results1[, chr38 := as.character(chr38)]
  results2[, chr38 := as.character(chr38)]
  
  setkey(results1, ensg, chr37, chr38, gene, minTss37, maxTss37, minTss38, maxTss38)
  
  both <- merge.data.table(results1, results2, on = c("ensg", "chr37", "chr38", "gene", "minTss37", "maxTss37", "minTss38", "maxTss38"), all = TRUE)
  
  both[cogs.x >= 0.5 & cogs.y < 0.5, ':=' (label = gene, 
                                         group = name1)]
  both[cogs.x < 0.5 & cogs.y >= 0.5, ':=' (label = gene, 
                                         group = name2)]
  both[cogs.x >= 0.5 & cogs.y >= 0.5, ':=' (label = gene, 
                                          group = "Standard and MultiCOGS")]

  # Fix NAs 
  both[cogs.x >= 0.5 & is.na(cogs.y), ':=' (group = name1,
                                           label = gene, 
                                           cogs.y = 0)]
  
  both[is.na(cogs.x) & cogs.y >= 0.5, ':=' (group = name2, 
                                            label = gene, 
                                            cogs.x = 0)]

  
  cogs1 <- paste0("cogs_", name1)
  cogs2 <- paste0("cogs_", name2)
  
  setnames(both, "cogs.x", cogs1)
  setnames(both, "cogs.y", cogs2)
  
  # Make a scatter plot.
  if(imageType == "pdf") {
    pdf(file = paste0(outdir, outname, ".pdf"), width = 12, height = 12)
    p <- ggplot(both, aes(x = get(cogs1), y = get(cogs2), label = label, colour = group))
    
    print(p + geom_point() +
            geom_text_repel(aes(label = label),
                            #box.padding   = 0.2, 
                            #point.padding = 0.2,
                            segment.color = 'grey50',
                            size = 4, 
                            max.overlaps = 20) + 
            theme(text = element_text(size = 20), panel.background = element_rect(fill = "white", colour = "grey50"), 
                  legend.position = "bottom") +
            xlab(name1) +
            ylab(name2)) 
    dev.off()
    } else { 
      if(imageType == "jpeg") {
        jpeg(file = paste0(outdir, outname, ".jpg"))
        p <- ggplot(both, aes(x = get(cogs1), y = get(cogs2), label = label, colour = group))
        print(p + geom_point() +
                geom_text_repel(aes(label = label),
                                #box.padding   = 0.2, 
                                #point.padding = 0.2,
                                segment.color = 'grey50',
                                size = 3, 
                                max.overlaps = 20) + 
                xlab(name1) +
                ylab(name2) +
                theme_classic())
        dev.off() 
        } else { if(imageType == "png") {
          png(file = paste0(outdir, outname, ".png"))
          p <- ggplot(both, aes(x = get(cogs1), y = get(cogs2), label = label, colour = group))
          print(p + geom_point() +
                  geom_text_repel(aes(label = label),
                                  #box.padding   = 0.2, 
                                  #point.padding = 0.2,
                                  segment.color = 'grey50',
                                  size = 3, 
                                  max.overlaps = 20) + 
                  xlab(name1) +
                  ylab(name2) +
                  theme_classic())
          dev.off() 
        } else { print("Please specify pdf, jpeg or png")
        }
        }
    }
  return(both)
}
#####



deLange_susie <- make_comparison_plot_annotated(res1 = "CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_ALL/Annotated_COGS_scores_data.table.txt",
                                                name1 = "Standard COGS", 
                                                res2 = "CD_deLange_ILCs_hg38_SuSIE_combinedInteractions_ALL/Annotated_COGS_scores_data.table.txt", 
                                                name2 = "MultiCOGS", 
                                                outdir = paste0(mydir, "comparisons/"), 
                                                outname = "CD_deLange_Classic_vs_SuSIE_Scatterplot", 
                                                imageType = "pdf")

