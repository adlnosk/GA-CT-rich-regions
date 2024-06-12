# script to get reads that are bridging the GA/CT regions (aligning on both sides)

fa_file = snakemake@input[["ass_regions"]]
fq_file = snakemake@input[["reads_regions"]]
multimapped_txt = snakemake@input[["multimapped_reads"]]

plot_out = snakemake@output[["plot"]]

path=snakemake@params[["path"]]
setwd(path)

# load libraries
library(dplyr)
library(ggplot2)

#fa_file="/work/project/briefwp3/Adela/GA_CT/results/Medicago_sativa_GACT_assembly.txt"
#fq_file="/work/project/briefwp3/Adela/GA_CT/results/Medicago_sativa_GACT_raw_reads.txt"
#sam="/work/project/briefwp3/Adela/GA_CT/results/mapped_reads/Diplocyclos_palmidus_fail_sorted.bam"

# mapped reads
fa <- read.table(fa_file)
colnames(fa)<-c("contigID", "fa_wst", "fa_wen", "fa_length", "sequence")
fa$window_ID < - paste0(fa$contigID, "_", fa$fa_wst, "_", fa$fa_wen, "_", fa$sequence)

b <- read.table(multimapped_txt, h=F)
colnames(b) <- c("readID", "FLAG", "contigID", "contig_start", "CIGAR", "read_length")

mapfa <- merge(fa,b, by.x="contigID") # multimapped reads with detected region

fq <- read.table(fq_file)
colnames(fq)<-c("readID", "fq_window_start", "fq_window_end", "fq_length", "sequence")

fq2 < - fq[, c(1,2,3,4)] #readID, start, end
c < - merge(fq2, mapfa, by="readID")

c$fq_wst < - c$fq_window_start + c$contig_start
c$fq_wen < - c$fq_window_end + c$contig_start


main < - c[,c("contigID", "readID", "fa_wst", "fa_wen", "fq_wst",
                     "fq_wen", "fa_length", "fq_length", "read_length", "contig_start")]
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


#    pdf(paste0(spec, "_facet_bridging_", n, "bp.pdf"),
#        width=13.8, height=7.31, onefile=TRUE)

    pdf(plot_out,width=13.8, height=7.31, onefile=TRUE)

    p1 < - ggplot(plmain, aes(x=pos2, y=DP_norm, group=wID, color=wID)) + geom_line() + theme_minimal() + facet_wrap(~wID) + theme(legend.position="none") + scale_x_discrete(labels=NULL, breaks=NULL) + labs(x="")

    p2 < - ggplot(plmain, aes(x=pos2, y=DP, group=wID, color=wID)) + geom_line() + theme_minimal() + facet_wrap(~wID) + theme(legend.position="none") + scale_x_discrete(labels=NULL, breaks=NULL) + labs(x="")
    print(p1)
    print(p2)

    dev.off()

#    system(paste0("convert ", spec, "_facet_bridging_", n, "bp.pdf ",
#           spec, "_facet_bridging_", n, "bp.png"), intern=TRUE)
}
