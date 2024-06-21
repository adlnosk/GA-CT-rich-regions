
########################################################################################
# Plot coverage curves, groupped by lengths, use samtools depth

########################################
## GACT_plot_coverages.R
# used in rules: plot_assembly_curves, plot_fail_curves 
# pass arguments

options(scipen=999)

fa_file = snakemake@input[["ass_regions"]]
fq_file = snakemake@input[["reads_regions"]]
sam = snakemake@input[["aligned_reads"]]
ref = snakemake@params[["ref"]]

plot_ends=snakemake@output[["plot_ends"]] # PDF with contig ENDS
table_plot=snakemake@output[["table_plot_all"]] # table with coordinates to plot all
table_ends=snakemake@output[["table_ends"]] # table with main info for contig ENDS
table_contigs=snakemake@output[["table_contigs"]]

spec = snakemake@params[["species"]]
outdepth=snakemake@params[["depth_prefix"]]
flank=snakemake@params[["plotting_flank"]]

path=snakemake@params[["path"]]
setwd(paste0(path,spec))


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

ofa$keep <- ofa$en - ofa$st
out <-  ofa[which(ofa$keep>0),c("contigID","st","en", "window_ID")]

write.table(out, paste0(outdepth,".bed"), col.names=F, row.names=F, quote=F)

# -------------------------------------------------------
# run samtools depth, wait for output
system(paste0("module load bioinfo/samtools/1.19; samtools depth -G UNMAP,SECONDARY -b ", outdepth, ".bed --reference ", ref ," ", sam, " -o ", outdepth, ".depth" ), intern = TRUE)
system(paste0("mkdir -p PLOTS"), intern = TRUE)
# OPTION TO FILTER DEPTHS (every second line): 
# system(paste0("awk -i inplace 'FNR%2' ", outdepth, ".depth"), intern = TRUE)
#--------------------------------------------------------

### clean 
rm(fa)

df <- read.delim(paste0(outdepth, ".depth"), h=F)
colnames(df) <- c("contigID", "pos", "DP")
df$a <- paste0(df$contigID, "_", df$pos)
ofa$a <- paste0(ofa$contigID, "_", ofa$st)

# create dictionary to assign windowID to samtools depth output

#contig_positions <- list()
#for (i in 1:nrow(ofa)) {
#  # Generate a sequence of positions from start to end
#  positions <- paste0(ofa$window_ID[i], ".", ofa$contigID[i],"_", seq(ofa$st[i], ofa$en[i]))
#  contig_positions[[i]] <- positions
#  }
# faster than the loop:
contig_positions <- vector("list", nrow(ofa))
contig_positions <- unlist(lapply(seq_len(nrow(ofa)), function(i) paste0(ofa$window_ID[i], ".", ofa$contigID[i], "_", seq(ofa$st[i], ofa$en[i], by=1))))


require(reshape2)
dd <- melt(contig_positions)
library(stringr)
convertor <- as.data.frame(str_split_fixed(dd[,1], fixed("."), 2)[, 1])
convertor[,2] <- str_split_fixed(dd[,1], fixed("."), 2)[, 2]
colnames(convertor) <- c("window_ID", "a")
dfcon <- df %>% left_join(convertor, by="a")
df_ofa <- dfcon %>% left_join(ofa, by="window_ID")


# SET PLOTTING POSITIONS
# regions found at contig starts (shorter than 70% of plotting size, starting with 0, descending DP)
size=(2*flank)*0.7
df_half_down <- df_ofa %>% group_by(wID) %>% mutate(pos = row_number()) %>% filter(max(pos) < size ) %>% filter(first(DP) > last(DP) ) %>% mutate(pos2  = row_number() )
table(df_half_down[,c("wID","length")])
hdlist <- unique(df_half_down[,"window_ID"])

######### PLOTTING CONTIG ENDS ############
print("Exporting PDF with contig ends - right and left.")
pdf(plot_ends,width=13.8, height=7.31, bg="white")

if(nrow(df_half_down)>0){
# ENDS downward
ggplot(df_half_down, aes(x=pos2, y=DP, group=wID, color=wID)) + geom_line() + xlim(0,(2*flank)) + theme_minimal() + facet_wrap(~length) + theme(legend.position="none")
}

# ENDS upward = regions at contig ends (shorter than 20 000, ascending DP)
df_half_up <- df_ofa %>% group_by(wID) %>% mutate(pos = row_number()) %>% filter(max(pos) < size ) %>% filter(first(DP) < last(DP) ) %>% mutate(pos2  = row_number() )
table(df_half_up[,c("wID","length")])
hulist <- unique(df_half_up[,"window_ID"])

if(nrow(df_half_up)>0){
ggplot(df_half_up, aes(x=pos2, y=DP, group=wID, color=wID)) + geom_line() + xlim(0,(2*flank)) + theme_minimal() + facet_wrap(~length) + theme(legend.position="none")
}

dev.off()


############## WRITING ENDS TABLE ###########
hdlist$side <- "R"
hulist$side <- "L"
ends_out <- rbind(hdlist, hulist)
print(nrow(ends_out))

library(tidyr)
ends_out2 <- ends_out %>% separate(window_ID, into = c("contigID", "start", "end", "sequence"), sep = "_")

print("Writing table with contig ends.")
write.table(ends_out2, table_ends, col.names=T, row.names=F, quote=F)

contig_pairs <- df %>%
  group_by(readID) %>%
  summarise(contig_pair = paste(sort(contigID), collapse = "_")) %>%
  ungroup()

# Count the number of readIDs for each contigID pair
contig_pair_counts <- contig_pairs %>%
  group_by(contig_pair) %>%
  summarise(count = n()) %>%
  ungroup()

# Split the contig_pair back into separate columns if needed
contig_pair_counts <- contig_pair_counts %>%
  separate(contig_pair, into = c("contigID1", "contigID2"), sep = "_")

# Print the resulting data frame
print(contig_pair_counts)
print("Writing table with contig pairs and number of supporting reads.")
write.table(contig_pair_counts, table_contigs, col.names=T, row.names=F, quote=F)



#############################################
# prepare for plotting - write.out the file
# regions of whole lengths - ONLY (20 000)
df_ofa2 <- df_ofa[! df_ofa$window_ID %in% ends_out$window_ID,]
df_rest <- df_ofa2 %>% group_by(wID) %>% mutate(pos = row_number()) %>% mutate(pos2  = row_number()+(((2*flank)+length)-max(pos)))
#df2 <- rbind(df_half_down, df_rest)
df2 <- df_rest

# SET NORMALISED DEPTH (per wID group)
# merge both
df3 <- df2 %>% group_by(wID) %>% mutate(DP_norm = (DP / max(DP)) * 100)

mu <- mean(df2$DP)
sd <- sd(df2$DP)
df2$zscore <- (df2$DP-mu)/sd

print("Writing table for plotting all.")
print(head(df3))
write.table(df3, table_plot, col.names=T, row.names=F, quote=F)




