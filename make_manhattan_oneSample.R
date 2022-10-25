#!/usr/bin/env Rscript

suppressMessages(library(argparser))
args = commandArgs(trailingOnly=T)

p <- arg_parser("Making rCOGS Manhattan", name="Rscript make_manhattan_multiset.R")
p <- add_argument(p, arg="--rCOGS1",
                  help="Full path to Annotated rCOGS results")
p = add_argument(p, arg="--outdir",
                 help="Full path to output directory", default=".")
p = add_argument(p, arg="--verbose",
                 help = "Flag specifying whether to print process steps.", flag=TRUE)

opts = parse_args(p, args)

rCOGS_results1 = opts[["rCOGS1"]]
outdir = opts[["outdir"]]
verb = opts[["verbose"]]

library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggrepel)


### 1. Read the COGS results
gene_res1 <- fread(rCOGS_results1)

### 2. Make the manhattan plot.
###### 
ifelse(!dir.exists(file.path(outdir, "Manhattan")), dir.create(file.path(outdir, "Manhattan")), FALSE)
######
### Note, this first bit gets rid of all chrs other than 1-22
gene_res1[, chr := as.numeric(chr38)]

## First get the cumulative position of each gene.
don <- as.data.table(gene_res1 %>% 
                       
                       # Compute chromosome size
                       group_by(chr) %>% 
                       summarise(chr_len=max(minTss38)) %>% 
                       
                       # Calculate cumulative position of each chromosome
                       mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
                       select(-chr_len) %>%
                       
                       # Add this info to the initial dataset
                       left_join(gene_res1, ., by=c("chr"="chr")) %>%
                       
                       # Add a cumulative position of each SNP
                       arrange(chr, minTss38) %>%
                       mutate( BPcum=minTss38+tot))

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of gene in bp, but just show the chromosome name instead.
axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Then make the plot using ggplot2
myname = basename(rCOGS_results1)
out_prefix_man = paste(outdir, "Manhattan", myname, sep = "/")

pdf(file = paste(out_prefix_man, "_Manhattan_score5.pdf", sep = ""), width = 10, height = 7)
# Make a list of gene names that we want to label
don[, label := ifelse(cogs >= 0.5, TRUE, FALSE)]

# Run ggplot
ggplot(don, aes(x=BPcum, y=cogs, label=gene)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=1, size=0.8) +
  scale_color_manual(values = rep(c("blue", "skyblue"), 22 )) +
  
  # custom axes:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) + ##here
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  ylim(0, 1) +
  xlab("Chromosome") +
  ylab("COGS score") +
  ggtitle(myname) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    text = element_text(size = 15)) +
  
  # Add cutoff line
  geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
  
  # Add labels for genes with score of > 0.3
  geom_text_repel(aes(label=ifelse(label == TRUE, as.character(gene),'')),
                  angle = 45, size = 3.5, max.time = 120.0, color = "black", force = 2, lineheight = 2, segment.size = 0.1, 
                  segment.alpha = 0.5, box.padding = 0.5, max.overlaps = Inf)

dev.off()
