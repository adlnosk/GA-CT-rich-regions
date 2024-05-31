
########################################################################################
# Plot coverage curves, groupped by lengths, use samtools depth

########################################
## GACT_plot_coverages.R
# used in rules: plot_assembly_curves, plot_fail_curves 
# pass arguments

spec = snakemake@params[["species"]]
flank = snakemake@params[["plotting_flank"]]

fa_file = snakemake@input[["ass_regions"]]
fq_file = snakemake@input[["reads_regions"]]
sam = snakemake@input[["aligned_reads"]]
ref = snakemake@input[["ref"]]

plotfile=snakemake@output[["plot"]]
plot_ends=snakemake@output[["plot_ends"]]

table=snakemake@output[["table"]]

outdepth=snakemake@params[["depth_prefix"]]

# load libraries
library(dplyr)
library(ggplot2)

# FASTA table
fa <- read.table(fa_file)
colnames(fa)<-c("contigID", "window_start", "window_end", "length", "sequence")
fa$window_ID <- paste0(fa$contigID,"_", fa$window_start,"_", fa$window_end,"_", fa$sequence)

print(table(fa$length))

ofa <- fa[fa$length>0,]

map <- ofa %>% group_by(contigID) %>% mutate(num = 1:n())
map$groupID <- paste0(map$contigID,"_",map$num)

ofa$st <- ofa$window_start - flank
ofa$en <- ofa$window_end + flank
ofa <- ofa[,c("window_ID","length","contigID","st","en")]
ofa[,"st"] <- ifelse(ofa[,"st"]<0,0,ofa[,"st"])
ofa$wID <- paste0(ofa$contigID,"_", ofa$st,"_", ofa$en)
out <-  ofa[,c("contigID","st","en", "window_ID")]

write.table(out, paste0(outdepth,".bed"), col.names=F, row.names=F, quote=F)

# -------------------------------------------------------
# run samtools depth, wait for output
system(paste0("module load bioinfo/samtools/1.19; samtools depth -G UNMAP,SECONDARY -b ", outdepth, ".bed --reference ", ref ," ", sam, " -o ", outdepth, ".depth" ), intern = TRUE)
#--------------------------------------------------------

df <- read.table(paste0(outdepth, ".depth"), h=F)
colnames(df) <- c("contigID", "pos", "DP")
df$a <- paste0(df$contigID, "_", df$pos)
ofa$a <- paste0(ofa$contigID, "_", ofa$st)

# create dictionary to assign windowID to samtools depth output
contig_positions <- list()
for (i in 1:nrow(ofa)) {
  # Generate a sequence of positions from start to end
  positions <- paste0(ofa$window_ID[i], ".", ofa$contigID[i],"_", seq(ofa$st[i], ofa$en[i]))
  contig_positions[[i]] <- positions
  }
require(reshape2)
dd <- melt(contig_positions)
library(stringr)
convertor <- as.data.frame(str_split_fixed(dd[,1], fixed("."), 2)[, 1])
convertor[,2] <- str_split_fixed(dd[,1], fixed("."), 2)[, 2]
colnames(convertor) <- c("window_ID", "a")
dfcon <- df %>% left_join(convertor, by="a")
df_ofa <- dfcon %>% left_join(ofa, by="window_ID")


# SET PLOTTING POSITIONS
# regions at contig starts (shorter than 20 000, starting with 0, descending DP)
df_half_down <- df_ofa %>% group_by(wID) %>% mutate(pos = row_number()) %>% filter(max(pos) < 15000 ) %>% filter(first(DP) > last(DP) ) %>% mutate(pos2  = row_number() )
table(df_half_down[,c("wID","length")])
hdlist <- unique(df_half_down[,"window_ID"])

pdf(plot_ends,width=13.8, height=7.31)

# pdf(paste0("/home/apoublan/public_html/reports/results/GA_CT_regions/",spec,"_facet_sizes_fa0_10flanks_down.pdf"),width=13.8, height=7.31)
ggplot(df_half_down, aes(x=pos2, y=DP, group=wID, color=wID)) + geom_line() + xlim(0,20000) + theme_minimal() + facet_wrap(~length) + theme(legend.position="none")
#dev.off()

# regions at contig ends (shorter than 20 000, ascending DP)
df_half_up <- df_ofa %>% group_by(wID) %>% mutate(pos = row_number()) %>% filter(max(pos) < 15000 ) %>% filter(first(DP) < last(DP) ) %>% mutate(pos2  = row_number() )
table(df_half_up[,c("wID","length")])
hulist <- unique(df_half_up[,"window_ID"])

#pdf(paste0("/home/apoublan/public_html/reports/results/GA_CT_regions/",spec,"_facet_sizes_fa0_10flanks_up.pdf"),width=13.8, height=7.31)
ggplot(df_half_up, aes(x=pos2, y=DP, group=wID, color=wID)) + geom_line() + xlim(0,20000) + theme_minimal() + facet_wrap(~length) + theme(legend.position="none")
dev.off()

# regions of whole lengths (20 000), or shorter (starting in the middle - ascending DP)
df_ofa2 <- df_ofa[! df_ofa$window_ID %in% hdlist$window_ID,]
df_rest <- df_ofa2 %>% group_by(wID) %>% mutate(pos = row_number()) %>% mutate(pos2  = row_number()+((20000+length)-max(pos)))
df2 <- rbind(df_half_down, df_rest)

# SET NORMALISED DEPTH (per wID group)
# merge both
df3 <- df2 %>% group_by(wID) %>% mutate(DP_norm = (DP / max(DP)) * 100)


mu <- mean(df2$DP)
sd <- sd(df2$DP)
df2$zscore <- (df2$DP-mu)/sd

# pdf("facet_sizes_fa200_10flanks.pdf", width=13.8, height=7.31)
#pdf(paste0("/home/apoublan/public_html/reports/results/GA_CT_regions/",spec,"_facet_sizes_fa0_10flanks.pdf"),width=13.8, height=7.31)

pdf(plotfile, width=13.8, height=7.31)
ggplot(df3, aes(x=pos2, y=DP_norm, group=wID, color=wID)) + geom_line() + theme_minimal() + facet_wrap(~length) + theme(legend.position="none") + scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")
ggplot(df3, aes(x=pos2, y=DP, group=wID, color=wID)) + geom_line() + theme_minimal() + facet_wrap(~length) + theme(legend.position="none") + scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")

dev.off()

write.table(df3, table, col.names=T, row.names=F, quote=F)

# system(paste0("convert /home/apoublan/public_html/reports/results/GA_CT_regions/",spec,"_facet_sizes_fa0_10flanks.pdf /home/apoublan/public_html/reports/results/GA_CT_regions/",spec,"_facet_sizes_fa0_10flanks.png"), intern = TRUE)



