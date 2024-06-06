
###########################################################################

# SEARCH FOR MULTIMAPPED READS wit clipped CIGAR, mapped on the edges- EXPORT TABLE WITH READ_ID, LENGTH, CONTIG_ID, CONTIG_POSITIONS (st and en), distance with S or H from CIGAR

## select only well-mapped reads to only 2 contig buts (not more, otherwise cannot be trusted)

# use reads multimapped: (primary and secondary => without unmapped and supplementary = 800 + 4) with soft/hard clip

#spec = "Medicago_sativa"

# pass arguments

spec = snakemake@params[["species"]]

fa_file = snakemake@input[["ass_regions"]]
fq_file = snakemake@input[["reads_regions"]]
sam = snakemake@input[["aligned_reads"]]
ref = snakemake@input[["ref"]]

locs_out = snakemake@output[["locs"]]
targets_out = snakemake@output[["targets"]]
multimapped_txt = snakemake@output[["multimapped_reads"]]


path=snakemake@params[["path"]]
setwd(path)


# load libraries
library(dplyr)
library(ggplot2)


# -------------------------------------------------------
# run samtools view to get the multimapped reads
system(paste0("module load bioinfo/samtools/1.19; samtools view -F 0x104 ", sam , " | awk '$6 ~ /H|S/{print $0}' | awk '{print $1, $2, $3, $4, $6, length ($10)-1 }'  > ", multimapped_txt), intern = TRUE)
#--------------------------------------------------------

# a <- read.table("fq_pos_list_all")
# colnames(a) <- c("readID", "window_start", "window_end", "sequence")
# b<- read.table("multimapped_fq_fa_pos.txt", h=F)
# colnames(b) <- c("readID", "FLAG", "contigID", "contig_start", "CIGAR", "read_len")
# c <- merge(a,b, by.x="readID") # multimapped with detected region

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
fai <- read.table(paste0("/work/project/briefwp3/Adela/", spec, "/assembly/", spec, ".asm.bp.p_ctg.fa.fai"), h=F)
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

write.table(targets[,c("readID", "contigID", "start", "end", "contig_start", "leftclip", "rightclip")], targets_out,  col.names=T, row.names=F, quote=F)  

# locs for IGV
out <- as.data.frame(targets[,c("contigID","readID", "st","en")])
out[,"st2"] <- formatC(out[,"st"], format="d", big.mark=",") #add comas into numbers
out[,"en2"] <- formatC(out[,"en"], format="d", big.mark=",")
locs <- paste0(out[,"contigID"],":",out[,"st2"],"-",out[,"en2"])
locs <- paste0("goto ", out[,"contigID"],":",out[,"st2"],"-",out[,"en2"],"\n","sort position\n", "collapse\n", "snapshot ", out[,"readID"], "_", out[,"contigID"], ".png\n")
head(locs)
dim(locs)
write.table(locs, locs_out, col.names=F, row.names=F, quote=F)

print("END OF SCRIPT")

############################################

# IGV batch script - insert in the beginning of `locs`:


#mkdir -p IGV_multimapped_edges

#cd IGV_multimapped_edges
#mv ../locs temp
#(echo -e new\\n\
#genome /work/project/briefwp3/Adela/${sp}/assembly/${sp}.asm.bp.p_ctg.fa\\n\
#snapshotDirectory /work/project/briefwp3/Adela/${sp}/assembly/GA_CT/IGV_shots\\n\
#load /work/project/briefwp3/Adela/${sp}/assembly/map_raw/${sp}_mapped_sorted.bam\\n ; cat temp ; echo exit) > locs
#rm temp
#
#rm *.png
#module load devel/java/17.0.6 bioinfo/IGV/2.16.1; igv.sh -b locs
#
#convert *png ${sp}_multi_edges.pdf
#mv ${sp}_multi_edges.pdf $HOME/public_html/reports/results/GA_CT_regions/




