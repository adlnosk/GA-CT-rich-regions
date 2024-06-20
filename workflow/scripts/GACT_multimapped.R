
###########################################################################

# SEARCH FOR MULTIMAPPED READS wit clipped CIGAR, mapped on the edges- EXPORT TABLE WITH READ_ID, LENGTH, CONTIG_ID, CONTIG_POSITIONS (st and en), distance with S or H from CIGAR

## select only well-mapped reads to only 2 contig buts (not more, otherwise cannot be trusted)

# use reads multimapped: (primary and secondary => without unmapped and supplementary = 800 + 4) with soft/hard clip

#spec = "Medicago_sativa"

# pass arguments

options(scipen=999)

spec = snakemake@params[["species"]]

fa_file = snakemake@input[["ass_regions"]]
fq_file = snakemake@input[["reads_regions"]]
ref = snakemake@input[["assembly"]]


locs_out = snakemake@output[["locs"]]
targets_out = snakemake@output[["targets"]]
multimapped_txt = snakemake@input[["multimapped_reads"]]

path=snakemake@params[["path"]]
setwd(path)

# load libraries
library(dplyr)
library(ggplot2)

# -------------------------------------------------------
# run samtools view to get the multimapped reads
# system(paste0("module load bioinfo/samtools/1.19; samtools view -F 0x104 ", sam , " | awk '$6 ~ /H|S/{print $0}' | awk '{print $1, $2, $3, $4, $6, length ($10)-1 }'  > ", multimapped_txt), intern = TRUE)
#--------------------------------------------------------

fq <- read.table(fq_file)
colnames(fq)<-c("readID", "window_start", "window_end", "length", "sequence")
a <- fq[,c("readID", "window_start", "window_end", "sequence")]

b <- read.table(multimapped_txt, h=F)
colnames(b) <- c("readID", "FLAG", "contigID", "contig_start", "CIGAR", "read_len")

c <- merge(a,b, by.x="readID") # multimapped reads with detected region

# number of soft clipped bases from CIGAR - left and right
#install.packages("gsubfn")
library(gsubfn)
c$leftclip <- rep(0,nrow(c))
c$rightclip <- rep(0,nrow(c))
c$leftclip <- apply(c,1,function(x){sum(as.numeric(strapply(as.character(x["CIGAR"]),"^(\\d+)[HS]", simplify=c)))})
c$rightclip <- apply(c,1,function(x){sum(as.numeric(strapply(as.character(x["CIGAR"]),"(\\d+)[HS]$",simplify=c)))})

c$start <- c$window_start+c$contig_start+c$leftclip # actual coordinates accounting for leftclip
c$end <- c$window_end+c$contig_start+c$leftclip # actual coordinates accounting for leftclip

# prepared for selection of reads with multiple records (or filter just before the end)
t <- table(c$readID)
tdup <- as.data.frame(t[t>0])
cdup <- c[is.element(c$readID, tdup$Var1),]
cdup$st <- cdup$window_start+cdup$contig_start - 10000 #for plotting
cdup$en <- cdup$window_end+cdup$contig_start + 10000
#cdup$window_ID <- paste0(cdup$readID, "_", cdup$contigID, "_", cdup$window_start,"_", cdup$window_end,"_", cdup$sequence)
cdup$window_length <- cdup$window_end - cdup$window_start

cout <- cdup[order(cdup[,"readID"]),c("readID", "window_length", "window_start", "window_end", "contigID", "start", "end", "st", "en", "leftclip", "rightclip", "contig_start" )]

# which reads are longer than beginning? - need lengths of contigs
fai <- read.table(paste0(ref, ".fai"), h=F)
colnames(fai) <- c("V1", "contig_max")

outfai <- merge(cout, fai[,c(1,2)], by.x="contigID", by.y="V1")

outfai[,"st"] <- ifelse(outfai[,"st"] < 0, 0, outfai[,"st"])


# select reads with either contig_start or dist_end < 50,000
# edges <- outfai[which((outfai$en + outfai$rightclip) > outfai$contig_max | (outfai$contig_start-outfai$leftclip) < 0),]

# reads with window start < 0 and end > max
edges <- outfai[which(outfai$end > outfai$contig_max | (outfai$start - outfai$leftclip) < outfai$contig_start ),]
target <- edges[order(edges$readID),]
targets <- target %>% group_by(readID) %>% distinct(readID, contigID, .keep_all = TRUE) %>% filter(n() == 2)


# clipped reads that are > once (are duplicated) and on the edge
# targets <- cdup[is.element(cdup$readID, edges$readID),]

print("Writing targets.")
write.table(targets[,c("readID", "contigID", "start", "end", "contig_start", "leftclip", "rightclip")], targets_out,  col.names=T, row.names=F, quote=F)  

# locs for IGV
#out <- as.data.frame(targets[,c("contigID","readID", "st","en")])
#out[,"st2"] <- formatC(out[,"st"], format="d", big.mark=",") #add comas into numbers
#out[,"en2"] <- formatC(out[,"en"], format="d", big.mark=",")
#out[,"readID2"] <- gsub("/", "_", out$readID)
#locs <- paste0(out[,"contigID"],":",out[,"st2"],"-",out[,"en2"])
#locs <- paste0("goto ", out[,"contigID"],":",out[,"st2"],"-",out[,"en2"],"\n","sort position\n", "collapse\n", "snapshot ", out[,"readID2"], "_", out[,"contigID"], ".png\n")
#print("Writing locations for IGV (locs).")
#write.table(locs, locs_out, col.names=F, row.names=F, quote=F)

print("END OF SCRIPT")

