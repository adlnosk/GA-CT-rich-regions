
########################################################################################
# Plot coverage curves, groupped by lengths, use samtools depth

########################################
## GACT_plot_coverages.R
# used in rules: plot_assembly_curves, plot_fail_curves 
# pass arguments

options(scipen=999)

intable = snakemake@input
plotfile = snakemake@output

path=snakemake@params[["path"]]
setwd(path)


# load libraries
library(dplyr)
library(ggplot2)

df3 <- read.table(intable, h=T)

pdf(plotfile, width=13.8, height=7.31)
ggplot(df3, aes(x=pos2, y=DP_norm, group=wID, color=wID)) + geom_line() + theme_minimal() + facet_wrap(~length) + theme(legend.position="none") + scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")
ggplot(df3, aes(x=pos2, y=DP, group=wID, color=wID)) + geom_line() + theme_minimal() + facet_wrap(~length) + theme(legend.position="none") + scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")

dev.off()


