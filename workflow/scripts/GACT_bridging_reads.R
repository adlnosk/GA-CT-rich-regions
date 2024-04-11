##################### reads vs contigs ###############
# spec = "Dactylis_glomerata"

# spec = "Medicago_sativa"

# pass arguments
args < - commandArgs(trailingOnly=TRUE)
print(args[1])

spec = args[1]


setwd(paste0("/work/project/briefwp3/Adela/", spec, "/assembly/GA_CT"))

# load libraries
library(dplyr)
library(ggplot2)

# get mapped reads (PRIMARY AND Supplementary) and fasta regions
# calculate number of reads in the fasta regions that are BRIDGING - aligning around the region


sam = paste0("/work/project/briefwp3/Adela/", spec,
             "/assembly/map_raw/", spec, "_mapped_sorted.bam")
print(sam)

# 0x100 = secondary, 0x4=unmapped
system(paste0("module load bioinfo/samtools/1.19; samtools view -F 0x104 ", sam,
       " | awk '{print $1, $3, $4, length($10)}' > fq_fa_fqstart.txt"), intern=TRUE)

# samtools view -F 0x804 ../map_raw/Medicago_sativa_mapped.sam | awk '{print $1, $3, $4, length($10)}' > fq_fa_fqstart.txt

# mapped reads
b < - read.table("fq_fa_fqstart.txt")
# FASTA table
fa < - read.table("fa_enriched_regions_GACT_slide.txt_IDs")
colnames(fa) < -c("readID", "contigID", "window_start",
                  "window_end", "length", "sequence")
fa$window_ID < - paste0(fa$contigID, "_", fa$window_start, "_", fa$window_end, "_", fa$sequence)
mapfa < - merge(b, fa, by.x="V2", by.y="contigID")
# reads regions
fq < - read.table("fq_enriched_regions_GACT_slide.txt_IDs")
fq$V6 < - sub('.', '', fq$V2)  # remove first character from the name to match
fq2 < - fq[, c(6, 3, 4)]
colnames(fq2) < - c("V1", "V2", "V3")
c < - merge(fq2, mapfa, by.x="V1", by.y="V1")
c$fq_length < - c$V3.x - c$V2.x
c$fq_wst < - c$V2.x + c$V3.y
c$fq_wen < - c$V3.x + c$V3.y

main < - c[, c("V2.y", "V1", "window_start", "window_end", "fq_wst", "fq_wen", "length", "fq_length", "V4", "V3.y")]
colnames(main) < - c("contigID", "readID", "fa_wst", "fa_wen", "fq_wst",
                     "fq_wen", "fa_length", "fq_length", "read_length", "contig_start")
# fq region = fa region:
main2 < - main[which(main$fq_wst <= main$fa_wst & main$fq_wen > main$fa_wen),]

plotdf < - read.table("plot_depths_fa200_10flanks.txt", h=T)

# loop over several flanking sequence sizes: 10 kb (= whole plotted region), 5 kb, 100 bp, 50 bp, 10 bp

for (n in c(5000, 100, 50, 10)) {

    print("###########################################################################")
    print(n)
    main3 < - main[which((main$contig_start) <= (main$fa_wst - n) & (main$contig_start + main$read_length) > (main$fa_wen + n)),]
    main3$wID < - paste0(main3$contigID, "_", (main3$fa_wst-10000), "_", (main3$fa_wen+10000))
    main4 < - main3 % > % group_by(contigID, fa_wst) % > % group_by(readID) % > % filter(row_number() == 1)
    main5 < - main4 % > % group_by(contigID, fa_wst) % > % mutate(nreads=n())
    print(nrow(main5))
    print(length(table(main5$contigID)))
    print(table(main5$contigID))

    # plot by group (wanna see if the non-drop in coverage is related to the number of reads bridging

    plmain < - merge(plotdf, main5, by="wID")


    pdf(paste0(spec, "_facet_bridging_", n, "bp.pdf"),
        width=13.8, height=7.31, onefile=TRUE)

    p1 < - ggplot(plmain, aes(x=pos2, y=DP_norm, group=wID, color=wID)) + geom_line() + theme_minimal() + facet_wrap(~wID) + theme(legend.position="none") + scale_x_discrete(labels=NULL, breaks=NULL) + labs(x="")

    p2 < - ggplot(plmain, aes(x=pos2, y=DP, group=wID, color=wID)) + geom_line() + theme_minimal() + facet_wrap(~wID) + theme(legend.position="none") + scale_x_discrete(labels=NULL, breaks=NULL) + labs(x="")
    print(p1)
    print(p2)

    dev.off()

    system(paste0("convert ", spec, "_facet_bridging_", n, "bp.pdf ",
           spec, "_facet_bridging_", n, "bp.png"), intern=TRUE)
}
