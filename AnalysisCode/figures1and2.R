# Run classification algorithm separately (Run mBE_pipeline.sh as per documentation) 
# */classified.bed obtained from classification algorithm 

###########################################################################
# Collect all signals used in classification into one single table

# ENCODE_cCREs_hg38-1.tab downloaded from UCSC Genome Browser's Table Browser https://genome.ucsc.edu/cgi-bin/hgTables
encode <- read.table("ENCODE_cCREs_hg38-1.tab", skip=1, sep="\t") %>%
    select(V4, V11) %>%
    set_colnames(c("name", "ccre"))
myclasses <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")

get_signals <- function(){
    rawdatas <- lapply(paste0(basedir, "/", signal_files), read.table, skip=3)
    rawdatas <- lapply(rawdatas, function(x){x %>% mutate(across(everything(), ~replace_na(.x, 0)))})
    datas <- rawdatas
    mydata <- bind_cols(lapply(datas, function(x){x %>% rowSums(.)}))
    colnames(mydata) <- signals
    mydata <- mydata %>% mutate(across(everything(), function(x) {log(((x*10000)/sum(x))+0.01)}))
    mydata <- mydata %>% mutate(across(everything(), function(x) {(x-mean(x))/sd(x)}))
    lojed <- read.table(loj_file) %>%
        select(V4,V10) %>%
        set_colnames(c("peakname", "name")) %>%
        left_join(encode) %>%
        group_by(peakname) %>% dplyr::slice(1) %>%
        replace_na(list(ccre = "notccre"))
    currpeaks <- read.table(currpeaks_file) %>%
        select(V1:V4) %>%
        set_colnames(c("chr", "start", "end", "peakname")) %>% 
        left_join(lojed)
    origpeaks <- read.table(origpeaks_file) %>%
        select(V1:V4) %>%
        set_colnames(c("chr", "start", "end", "peakname"))
    mydata <- mydata %>% bind_cols(currpeaks)
    plotdata <- mydata %>% filter(ccre != "notccre")
    classified <- read.table(paste0(basedir, "/", "classified.bed")) %>% set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass"))
    plotdata <- plotdata %>% bind_cols(classified %>% select(newclass)) %>% mutate(newclass=factor(newclass, levels = myclasses)) %>% select(-c(chr, start, end)) %>% left_join(origpeaks)
    return(plotdata %>% select(-c(peakname, name, chr, start, end, ccre)) %>% pivot_longer(!newclass, names_to = "mark", values_to = "z"))
    # return(plotdata %>% select(-c(peakname, name, chr, start, end, ccre)))
}

basedir <- "/HMEC/"
signal_files <- c("H3K27ac", "H3K4me1", "H3K4me3_signal", "H2Az", "H3K27me3", "mH2A1", "mH2A2", "CTCF", "H4K12ac")
signals <- c("H3K27ac", "H3K4me1", "H3K4me3", "H2Az", "H3K27me3", "mH2A1", "mH2A2", "CTCF", "H4K12ac")
loj_file <- "HMEC_loj_cCREs.bed"
currpeaks_file <- "HMEC-ATAC-K4m1-outBL-resized.bed"
origpeaks_file <- "HMEC-ATAC-K4m1-outBL.bed"
HMEC <- get_signals()

basedir <- "/NHM/"
signal_files <- c("mH2A2", "H3K27ac", "H3K4me1", "H3K4me3_correct_signal", "H2Az", "H3K27me3", "mH2A1", "H4K12ac")
signals <- c("mH2A2", "H3K27ac", "H3K4me1", "H3K4me3", "H2Az", "H3K27me3", "mH2A1", "H4K12ac")
loj_file <- "NHM_loj_cCREs.bed"
currpeaks_file <- "NHM-ATAC-K4m1-outBL-resized.bed"
origpeaks_file <- "NHM-ATAC-K4m1-outBL.bed"
NHM <- get_signals()

basedir <- "/MCF7/"
signal_files <- c("mH2A2", "H3K27ac", "H3K4me1", "H3K4me3_signal", "H3K27me3", "mH2A1", "CTCF", "H4K12acVeh", "H4K12acE2")
signals <- c("mH2A2", "H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "mH2A1", "CTCF", "H4K12acVeh", "H4K12acE2")
loj_file <- "MCF7_loj_cCREs.bed"
currpeaks_file <- "MCF7-ATAC-K4m1-outBL-resized.bed"
origpeaks_file <- "MCF7-ATAC-K4m1-outBL.bed"
MCF7 <- get_signals()

basedir <- "/231L/"
signal_files <- c("H3K27ac", "H3K4me1", "H3K27me3", "ATAC")
signals <- c("H3K27ac", "H3K4me1", "H3K27me3", "ATAC")
loj_file <- "231L_loj_cCREs.bed"
currpeaks_file <- "231L-ATAC-K4m1-outBL-resized.bed"
origpeaks_file <- "231L-ATAC-K4m1-outBL.bed"
v231L <- get_signals()

basedir <- "/HepG2/"
signal_files <- c("H3K27ac", "H3K4me1", "H3K4me3_signal", "H3K27me3", "mH2A1", "mH2A2", "CTCF")
signals <- c("H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "mH2A1", "mH2A2", "CTCF")
loj_file <- "HepG2_loj_cCREs.bed"
currpeaks_file <- "HepG2-ATAC-K4m1-outBL-resized.bed"
origpeaks_file <- "HepG2-ATAC-K4m1-outBL.bed"
HepG2 <- get_signals()

all_signals <- bind_rows(list(HMEC = HMEC, NHM = NHM, MCF7 = MCF7, `231L` = v231L, HepG2 = HepG2), .id = "CellLine")
write.table(all_signals, file = "all_signals.tab", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#########################################  Fig 1d and Supp Fig 1d  ###########################################################

all_signals <- read.table(file = "all_signals.tab")
colnames(all_signals) <- c("CellLine", "newclass", "mark", "z")
plotdata <- all_signals %>% filter(!(mark %in% c("ATAC", "H4K12acVeh", "H4K12acE2", "H4K12ac"))) %>% 
    mutate(mark = factor(mark, levels = c("H3K27ac", "H3K4me1", "H2Az", "H3K4me3", "H3K27me3", "mH2A1", "mH2A2", "CTCF"))) %>%
    mutate(CellLine = factor(CellLine, levels = c("HMEC", "NHM", "MCF7", "231L", "HepG2")))

# Mann-Whitney test with Bonferroni correction to test for significance of enrichment
combos <- expand.grid(c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE"), c("H3K4me1", "H3K27ac", "H2Az", "H3K4me3", "H3K27me3", "CTCF", "mH2A1", "mH2A2"), c("HMEC", "NHM", "MCF7", "HepG2"))
combos_split <- split(combos, seq(nrow(combos)))
get_tests <- function(combo){
    print(paste0(as.character(combo$Var1), " ", as.character(combo$Var2), " ", as.character(combo$Var3)))
    class_val <- as.character(combo$Var1)
    mark_val <- as.character(combo$Var2)
    cl_val <- as.character(combo$Var3)
    if (!(mark_val %in% (curr %>% filter(CellLine == cl_val))$mark)) {return (data.frame())}
    temp <- curr %>% filter(CellLine == cl_val & mark == mark_val)
    x <- temp %>% filter(newclass == class_val) %>% select(z) %>% unlist() %>% unname()
    y <- temp %>% filter(newclass != class_val) %>% select(z) %>% unlist() %>% unname()
    p <- wilcox.test(x = x, y = y)$p.value
    if (cl_val == "HMEC"){ padj <- p * 40 }else{ padj <- p * 35 } # Bonferroni correction
    notsigbool <- padj >= 0.05
    return(data.frame(CellLine = cl_val, newclass = class_val, mark = mark_val, notsig = notsigbool))
}
curr <- plotdata
notsig <- bind_rows(lapply(combos_split, get_tests))
notsig <- notsig %>% filter(notsig == TRUE) %>% mutate(notsig = "*")

plotdata1 <- plotdata %>% group_by(CellLine, newclass, mark) %>% summarize(median_z = median(z)) %>% left_join(notsig) %>% filter(CellLine != "231L") %>%
    mutate(newclass=replace(newclass, newclass=="H3K4me3", "APL")) %>%
    mutate(newclass=replace(newclass, newclass=="ATAConly", "ATAC-only")) %>%
    mutate(newclass = factor(newclass, levels = rev(c("Active", "APL", "ATAC-only", "Inactive", "mBE")))) %>%
    mutate(CellLine = factor(CellLine, levels = c("HMEC", "NHM", "MCF7", "HepG2"))) %>%
    mutate(mark = factor(mark, levels = c("H3K4me1", "H3K27ac", "H2Az", "H3K4me3", "H3K27me3", "CTCF", "mH2A1", "mH2A2")))

pdf(file='fig1d.pdf', width=8, height=3)
p1 <- ggplot(plotdata1, aes(x = mark, y = newclass)) +
    facet_grid(cols = vars(CellLine), scales = "free_x") +
    geom_tile(aes(fill = median_z)) +
    geom_point(data = plotdata1 %>% filter(notsig == "*"), shape=1, size=1) +
    # coord_equal() +
    scale_fill_gradient2(mid = "lightgray", low = "blue", high = "red", limits = c(-2.05, 2.05)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), strip.background = element_blank())
p1
dev.off()

plotdata1 <- plotdata %>% filter(CellLine != "231L") %>%
    mutate(newclass=replace(newclass, newclass=="H3K4me3", "APL")) %>%
    mutate(newclass=replace(newclass, newclass=="ATAConly", "ATAC-only")) %>%
    mutate(newclass = factor(newclass, levels = c("Active", "APL", "ATAC-only", "Inactive", "mBE"))) %>%
    mutate(CellLine = factor(CellLine, levels = c("HMEC", "NHM", "MCF7", "HepG2"))) %>%
    mutate(mark = factor(mark, levels = c("H3K4me1", "H3K27ac", "H2Az", "H3K4me3", "H3K27me3", "CTCF", "mH2A1", "mH2A2")))
pdf(file='fig1d-supp.pdf', width=12, height=10)
ggplot(plotdata1, aes(x = mark, y = z)) +
    facet_grid(cols = vars(newclass), rows = vars(CellLine)) +
    geom_boxplot(outlier.alpha = 0.1, outlier.shape = NA) +
    coord_cartesian(ylim = c(-4, 4)) +
    # theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

############ Fig 1e and Supp Fig 2a and 2b ###########################

# *.tab.gz files from computeMatrix command used in mBE_pipeline.sh for HMEC
# DeepTools plotHeatmap is used for plotting heatmaps
# plotHeatmap -m ATAC.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out ATAC.pdf --samplesLabel ATAC
# plotHeatmap -m H3K4me1.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H3K4me1.pdf --samplesLabel H3K4me1
# plotHeatmap -m H3K4me3.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H3K4me3.pdf --samplesLabel H3K4me3
# plotHeatmap -m H3K27ac.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H3K27ac.pdf --samplesLabel H3K27ac
# plotHeatmap -m H2Az.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H2Az.pdf --samplesLabel H2Az
# plotHeatmap -m H3K27me3.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H3K27me3.pdf --samplesLabel H3K27me3
# plotHeatmap -m mH2A1.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out mH2A1.pdf --samplesLabel mH2A1
# plotHeatmap -m mH2A2.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out mH2A2.pdf --samplesLabel mH2A2
# plotHeatmap -m CTCF.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out CTCF.pdf --samplesLabel CTCF
# plotHeatmap -m Conserv.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out Conserv.pdf --samplesLabel Conserv
# plotHeatmap -m DNase.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out DNase.pdf --samplesLabel DNase
# plotHeatmap -m DNAM.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out DNAM.pdf --samplesLabel '5mC'
# plotHeatmap -m DNAM5hmC.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out DNAM5hmC.pdf --samplesLabel '5hmC'
# plotHeatmap -m H3K36me3.tab.gz --colorList 'white,blue' --sortRegions keep -out H3K36me3.pdf --samplesLabel H3K36me3
# plotHeatmap -m H2BK12ac.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H2BK12ac.pdf --samplesLabel H2BK12ac
# plotHeatmap -m H2BK120ac.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H2BK120ac.pdf --samplesLabel H2BK120ac

############ Fig 1f and Supp Fig 1c (Homer annotation of classified CRE) #################

HMEC <- read.table("HMEC_classes.bed") %>% set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass"))
write.table(HMEC %>% select(chr, start, end, peakname) %>% mutate(strand = "."), file = "HMEC.classified.homer.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
MCF7 <- read.table("MCF7_classes.bed") %>% set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass"))
write.table(MCF7 %>% select(chr, start, end, peakname) %>% mutate(strand = "."), file = "MCF7.classified.homer.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
NHM <- read.table("NHM_classes.bed") %>% set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass"))
write.table(NHM %>% select(chr, start, end, peakname) %>% mutate(strand = "."), file = "NHM.classified.homer.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
HepG2 <- read.table("HepG2_classes.bed") %>% set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass"))
write.table(HepG2 %>% select(chr, start, end, peakname) %>% mutate(strand = "."), file = "HepG2.classified.homer.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
DF <- read.table("DF_classes.bed") %>% set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass"))
write.table(DF %>% select(chr, start, end, peakname) %>% mutate(strand = "."), file = "DF.classified.homer.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# cd /mBE/figures
# annotatePeaks.pl HMEC.classified.homer.bed hg19 -annStats HMEC.annStats > HMEC.ann 2> 2 &
# annotatePeaks.pl MCF7.classified.homer.bed hg19 -annStats MCF7.annStats > MCF7.ann 2> 2 &
# annotatePeaks.pl NHM.classified.homer.bed hg19 -annStats NHM.annStats > NHM.ann 2> 2 &
# annotatePeaks.pl HepG2.classified.homer.bed hg19 -annStats HepG2.annStats > HepG2.ann 2> 2 &
# annotatePeaks.pl DF.classified.homer.bed mm9 -annStats DF.annStats > DF.ann 2> 2 &

annlevels <- c("promoter-TSS", "TTS", "5' UTR", "exon", "intron", "3' UTR", "Intergenic", "non-coding")

ann <- mapply(
    function(annFile, classified){
        read.table(annFile, sep = "\t", skip = 1, quote = "", comment.char = "") %>% 
            select(V1, V8) %>% 
            set_colnames(c("peakname", "annotation")) %>%
            separate(annotation, c("annotation", "throw"), sep = " \\(") %>%
            select(-throw) %>%
            mutate(annotation = factor(annotation, levels = annlevels)) %>%
            left_join(classified %>% select(peakname, newclass))
    },
    c("HMEC.ann", "MCF7.ann", "NHM.ann", "HepG2.ann", "DF.ann"),
    list(HMEC, MCF7, NHM, HepG2, DF), SIMPLIFY = FALSE
)
names(ann) <- c("HMEC", "MCF7", "NHM", "HepG2", "DF")
ann <- bind_rows(ann, .id = "cellline")
ann <- ann %>% 
    mutate(cellline = factor(cellline, levels = c("HMEC", "MCF7", "NHM", "HepG2", "DF"))) %>%
    mutate(newclass = factor(newclass, levels = c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE"))) %>%
    mutate(annotation = factor(annotation, levels = annlevels))
# anntable <- data.frame(table(ann %>% select(cellline, annotation, newclass))) %>% mutate(cellline = factor(cellline, levels = c("HMEC", "MCF7", "NHM", "HepG2", "DF")))

library(ggpie)

plots <- lapply(c("HMEC", "MCF7", "NHM", "HepG2", "DF"), function(x){if (x != "DF"){ggpie(ann %>% filter(cellline == x), annotation, newclass, nrow=1, label.size=0, border.color="white") + theme(legend.position = "none") + ggtitle(x)}else{ggpie(ann %>% filter(cellline == x), annotation, newclass, nrow=1, label.size=0, border.color="white") + theme(legend.position = "bottom") + ggtitle(x)}})

pdf(file='homer_pies.pdf', width=12, height=12)
wrap_plots(plots, ncol = 1)
dev.off()

###################### Fig 2a (eRNA estimation from RNAseq data) #################

# wget https://www.encodeproject.org/files/ENCFF284MYQ/@@download/ENCFF284MYQ.bam -O HMEC_polyA-rna.bam
# wget https://www.encodeproject.org/files/ENCFF499SYJ/@@download/ENCFF499SYJ.bam -O HMEC_total-rna.bam
# wget https://www.encodeproject.org/files/ENCFF238UQC/@@download/ENCFF238UQC.bam -O HMEC_small-rna.bam

# bedtools multicov -bams HMEC_polyA-rna.bam -bed HMEC-ATAC-K4m1-outBL.bed > HMEC_polyA-rna.encode.bed &
# bedtools multicov -bams HMEC_total-rna.bam -bed HMEC-ATAC-K4m1-outBL.bed > HMEC_total-rna.encode.bed &
# bedtools multicov -bams HMEC_small-rna.bam -bed HMEC-ATAC-K4m1-outBL.bed > HMEC_small-rna.encode.bed &

myclasses <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")
color_key <- c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")

origpeaks <- read.table("HMEC-ATAC-K4m1-outBL.bed") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname"))

files <- c("HMEC_polyA-rna.encode.bed", "HMEC_total-rna.encode.bed", "HMEC_small-rna.encode.bed")
samplenames <- c("enc.polyA", "enc.total", "enc.small")

exprs <- origpeaks %>%
    bind_cols(lapply(files, function(x){read.table(x) %>% select(V5)})) %>%
    set_colnames(c(colnames(origpeaks), samplenames)) %>%
    mutate(length = end-start+1)

ccre_exprs <- read.table("HMEC_classes.bed") %>% 
    set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass")) %>% 
    select(peakname, newclass) %>%
    left_join(exprs)

d0 <- DGEList(ccre_exprs %>% select(enc.polyA, enc.total, enc.small))
# d0 <- calcNormFactors(d0)
# cutoff <- 1
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# d <- d0[-drop,] 
# filtered_ccre_exprs <- ccre_exprs[-drop,] 
d <- d0
filtered_ccre_exprs <- ccre_exprs

logRPKM <- data.frame(rpkm(d, filtered_ccre_exprs$length, log=TRUE), newclass = filtered_ccre_exprs$newclass) %>%
    pivot_longer(!newclass, names_to="dataset", values_to="logRPKM") %>%
    arrange(dataset) %>% 
    mutate(newclass = factor(newclass, levels = myclasses)) %>%
    mutate(dataset = factor(dataset, levels = c("enc.polyA", "enc.total", "enc.small")))

p1 <- ggplot(logRPKM %>% filter(dataset == "enc.polyA") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.1, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(-6), 2^(15))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()

pdf(file='HMEC_expr_polyARNA_ENCFF284MYQ.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

p1 <- ggplot(logRPKM %>% filter(dataset == "enc.total") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.1, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(-6), 2^(15))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()

pdf(file='HMEC_expr_totalRNA_ENCFF499SYJ.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

p1 <- ggplot(logRPKM %>% filter(dataset == "enc.small") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.1, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(-6), 2^(5))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()

pdf(file='HMEC_expr_smallRNA_ENCFF238UQC.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()


#####

# fasterq-dump SRR5228548
# ref=${references}/STAR/hg19
# STAR --runMode alignReads --runThreadN 16 --genomeDir ${ref} --readFilesIn SRR5228548.1_1.fastq SRR5228548.1_2.fastq --outFileNamePrefix SRR5228548_ --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
# samtools index SRR5228548_Aligned.sortedByCoord.out.bam
# bedtools multicov -bams SRR5228548_Aligned.sortedByCoord.out.bam -bed NHM-ATAC-K4m1-outBL.bed > SRR5228548.enh.bed &


myclasses <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")
color_key <- c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")

origpeaks <- read.table("NHM-ATAC-K4m1-outBL.bed") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname"))

files <- c("SRR5228548.enh.bed")
samplenames <- c("SRR5228548")

exprs <- origpeaks %>%
    bind_cols(lapply(files, function(x){read.table(x) %>% select(V5)})) %>%
    set_colnames(c(colnames(origpeaks), samplenames)) %>%
    mutate(length = end-start+1)

ccre_exprs <- read.table("NHM_classes.bed") %>% 
    set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass")) %>% 
    select(peakname, newclass) %>%
    left_join(exprs)

d0 <- DGEList(ccre_exprs %>% select(SRR5228548))
d <- d0
filtered_ccre_exprs <- ccre_exprs

logRPKM <- data.frame(rpkm(d, filtered_ccre_exprs$length, log=TRUE), newclass = filtered_ccre_exprs$newclass) %>%
    pivot_longer(!newclass, names_to="dataset", values_to="logRPKM") %>%
    arrange(dataset) %>% 
    mutate(newclass = factor(newclass, levels = myclasses))

p1 <- ggplot(logRPKM %>% filter(dataset == "SRR5228548") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(1), 2^(16))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()
    # ggtitle("SRR5228548")

pdf(file='NHM_expr_totalRNA_SRR5228548.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

##########

# wget https://www.encodeproject.org/files/ENCFF712EIK/@@download/ENCFF712EIK.bam -O MCF7_polyA-rna.bam
# wget https://www.encodeproject.org/files/ENCFF599PRI/@@download/ENCFF599PRI.bam -O MCF7_small-rna.bam
# wget https://www.encodeproject.org/files/ENCFF522BMH/@@download/ENCFF522BMH.fastq.gz -O MCF7_micro-rna.fastq.gz

# ref=${references}/STAR/hg19
# STAR --runMode alignReads --runThreadN 16 --genomeDir ${ref} --readFilesIn MCF7_micro-rna.fastq.gz --outFileNamePrefix MCF7_micro-rna_ --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat
# samtools index MCF7_micro-rna_Aligned.sortedByCoord.out.bam

# bedtools multicov -bams MCF7_polyA-rna.bam -bed MCF7-ATAC-K4m1-outBL.bed > MCF7_polyA-rna.encode.bed &
# bedtools multicov -bams MCF7_small-rna.bam -bed MCF7-ATAC-K4m1-outBL.bed > MCF7_small-rna.encode.bed &
# bedtools multicov -bams MCF7_micro-rna_Aligned.sortedByCoord.out.bam -bed MCF7-ATAC-K4m1-outBL.bed > MCF7_micro-rna.encode.bed &

myclasses <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")
color_key <- c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")

origpeaks <- read.table("MCF7-ATAC-K4m1-outBL.bed") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname"))

files <- c("MCF7_polyA-rna.encode.bed", "MCF7_small-rna.encode.bed", "MCF7_micro-rna.encode.bed")
samplenames <- c("enc.polyA", "enc.small", "enc.micro")

exprs <- origpeaks %>%
    bind_cols(lapply(files, function(x){read.table(x) %>% select(V5)})) %>%
    set_colnames(c(colnames(origpeaks), samplenames)) %>%
    mutate(length = end-start+1)

ccre_exprs <- read.table("MCF7_classes.bed") %>% 
    set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass")) %>% 
    select(peakname, newclass) %>%
    left_join(exprs)

d0 <- DGEList(ccre_exprs %>% select(enc.polyA, enc.small, enc.micro))
d <- d0
filtered_ccre_exprs <- ccre_exprs

logRPKM <- data.frame(rpkm(d, filtered_ccre_exprs$length, log=TRUE), newclass = filtered_ccre_exprs$newclass) %>%
    pivot_longer(!newclass, names_to="dataset", values_to="logRPKM") %>%
    arrange(dataset) %>% 
    mutate(newclass = factor(newclass, levels = myclasses)) %>%
    mutate(dataset = factor(dataset, levels = c("enc.polyA", "enc.small", "enc.micro")))

p1 <- ggplot(logRPKM %>% filter(dataset == "enc.small") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.1, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(-4), 2^(0))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()
    # ggtitle("ENCFF599PRI")

pdf(file='MCF7_expr_smallRNA_ENCFF599PRI.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

p1 <- ggplot(logRPKM %>% filter(dataset == "enc.polyA") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.1, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(-5), 2^(15))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()
    # ggtitle("ENCFF599PRI")

pdf(file='MCF7_expr_polyARNA_ENCFF712EIK.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

p1 <- ggplot(logRPKM %>% filter(dataset == "enc.micro") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.1, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(-6), 2^(-1))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()
    # ggtitle("ENCFF599PRI")

pdf(file='MCF7_expr_microRNA_ENCFF522BMH.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

#######

# wget -nv https://www.encodeproject.org/files/ENCFF642RLH/@@download/ENCFF642RLH.bam -O HepG2_small-rna.bam &
# wget -nv https://www.encodeproject.org/files/ENCFF770XVY/@@download/ENCFF770XVY.bam -O HepG2_total-rna.bam &
# wget -nv https://www.encodeproject.org/files/ENCFF989ZMQ/@@download/ENCFF989ZMQ.bam -O HepG2_polyA-rna.bam &

# bedtools multicov -bams HepG2_polyA-rna.bam -bed HepG2-ATAC-K4m1-outBL.bed > HepG2_polyA-rna.encode.bed &
# bedtools multicov -bams HepG2_total-rna.bam -bed HepG2-ATAC-K4m1-outBL.bed > HepG2_total-rna.encode.bed &
# bedtools multicov -bams HepG2_small-rna.bam -bed HepG2-ATAC-K4m1-outBL.bed > HepG2_small-rna.encode.bed &

myclasses <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")
color_key <- c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")

origpeaks <- read.table("HepG2-ATAC-K4m1-outBL.bed") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname"))

files <- c("HepG2_polyA-rna.encode.bed", "HepG2_total-rna.encode.bed", "HepG2_small-rna.encode.bed")
samplenames <- c("enc.polyA", "enc.total", "enc.small")

exprs <- origpeaks %>%
    bind_cols(lapply(files, function(x){read.table(x) %>% select(V5)})) %>%
    set_colnames(c(colnames(origpeaks), samplenames)) %>%
    mutate(length = end-start+1)

ccre_exprs <- read.table("HepG2_classes.bed") %>% 
    set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass")) %>% 
    select(peakname, newclass) %>%
    left_join(exprs)

d0 <- DGEList(ccre_exprs %>% select(enc.polyA, enc.total, enc.small))
d <- d0
filtered_ccre_exprs <- ccre_exprs

logRPKM <- data.frame(rpkm(d, filtered_ccre_exprs$length, log=TRUE), newclass = filtered_ccre_exprs$newclass) %>%
    pivot_longer(!newclass, names_to="dataset", values_to="logRPKM") %>%
    arrange(dataset) %>% 
    mutate(newclass = factor(newclass, levels = myclasses)) %>%
    mutate(dataset = factor(dataset, levels = c("enc.polyA", "enc.total", "enc.small")))

p1 <- ggplot(logRPKM %>% filter(dataset == "enc.total") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.1, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(-6), 2^(15))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()

pdf(file='HepG2_expr_totalRNA_ENCFF770XVY.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

p1 <- ggplot(logRPKM %>% filter(dataset == "enc.polyA") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(-6), 2^(15))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()

pdf(file='HepG2_expr_PolyARNA_ENCFF989ZMQ.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

p1 <- ggplot(logRPKM %>% filter(dataset == "enc.small") %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^(-6), 2^(15))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()

pdf(file='HepG2_expr_smallRNA_ENCFF642RLH.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

###################### Salmon ####################################

# salmon quant -i /references/gtfs/Hsapiens_index -l A -1 /mBE/Encode2021/HMEC/total-rna_R1.fastq.gz -2 /mBE/Encode2021/HMEC/total-rna_R2.fastq.gz -p 8 --validateMappings -o /mBE/HMEC1/salmon_quant/totalRNA 1> 1 2> 2 &
# salmon quant -i /references/gtfs/Hsapiens_index -l A -1 /mBE/Encode2021/MCF7/polyA-rna_R1.fastq.gz -2 /mBE/Encode2021/MCF7/polyA-rna_R2.fastq.gz -p 8 --validateMappings -o /mBE/MCF7/salmon_quant/polyARNA 1> 1 2> 2 &
# salmon quant -i /references/gtfs/Hsapiens_index -l A -1 /mBE/Encode2021/MCF7/polyA-rna_rep2_R1.fastq.gz -2 /mBE/Encode2021/MCF7/polyA-rna_rep2_R2.fastq.gz -p 8 --validateMappings -o /mBE/MCF7/salmon_quant/polyARNA_rep2 1> 1 2> 2 &
# salmon quant -i /references/gtfs/Hsapiens_index -l A -1 /mBE/HepG2/total-rna_R1.fastq.gz -2 /mBE/HepG2/total-rna_R2.fastq.gz -p 8 --validateMappings -o /mBE/HepG2/salmon_quant/totalRNA 1> 1 2> 2 &
# salmon quant -i /references/gtfs/Hsapiens_index -l A -1 /mBE/HepG2/total-rna_rep2_R1.fastq.gz -2 /mBE/HepG2/total-rna_rep2_R2.fastq.gz -p 8 --validateMappings -o /mBE/HepG2/salmon_quant/totalRNA_rep2 1> 1 2> 2 &

# salmon quant -i /references/gtfs/Hsapiens_index -l A -1 /RNAseq-data/NHM/SRR5228548_1.fastq -2 /RNAseq-data/NHM/SRR5228548_2.fastq -p 8 --validateMappings -o /mBE/NHM/salmon_quant/totalRNA 1> 1 2> 2 &
# salmon quant -i /references/gtfs/Hsapiens_index -l A -1 /RNAseq-data/NHM/SRR5228549_1.fastq -2 /RNAseq-data/NHM/SRR5228549_2.fastq -p 8 --validateMappings -o /mBE/NHM/salmon_quant/totalRNA_rep2 1> 1 2> 2 &

# macroH2A1.1 = XM_005272132.2
# https://www.ncbi.nlm.nih.gov/nuccore/XM_005272132.2
# http://useast.ensembl.org/Homo_sapiens/Transcript/Summary?g=ENSG00000113648;r=5:135334410-135399211;t=ENST00000312469

# macroH2A1.2 = NM_004893.3
# https://www.ncbi.nlm.nih.gov/nuccore/NM_004893
# http://useast.ensembl.org/Homo_sapiens/Transcript/Summary?g=ENSG00000113648;r=5:135334383-135399887;t=ENST00000304332

setwd("/mBE/")

files <- c(
    "/mBE/HMEC1/salmon_quant/totalRNA/quant.sf", 
    "/mBE/NHM/salmon_quant/totalRNA/quant.sf",
    "/mBE/NHM/salmon_quant/totalRNA_rep2/quant.sf",
    "/mBE/MCF7/salmon_quant/polyARNA/quant.sf",
    "/mBE/MCF7/salmon_quant/polyARNA_rep2/quant.sf",
    "/mBE/HepG2/salmon_quant/totalRNA/quant.sf",
    "/mBE/HepG2/salmon_quant/totalRNA_rep2/quant.sf"
)

samplenames <- c("HMEC_Rep1", "NHM_Rep1", "NHM_Rep2", "MCF7_Rep1", "MCF7_Rep2", "HepG2_Rep1", "HepG2_Rep2")
celllines <- c("HMEC", "NHM", "MCF7", "HepG2")

exprs <- read.table(files[[1]], skip = 1) %>% select(V1) %>%
    bind_cols(lapply(files, function(x){read.table(x, skip = 1) %>% select(V4)})) %>%
    set_colnames(c("TXID", samplenames))

TXs <- c("ENST00000312469.4", "ENST00000304332.4", "ENST00000373255.4")
# TXs <- c("ENST00000312469.4", "ENST00000304332.4")
variants <- c("mH2A1.1", "mH2A1.2", "mH2A2")
mapping <- data.frame(TXID = TXs, variant = variants)

exprs1 <- exprs %>% 
    pivot_longer(!TXID, names_to = "sample", values_to = "TPM")
exprs2 <- exprs1 %>%
    separate(col = "sample", sep = "_", into = c("cellline", "replicate")) %>%
    filter(TXID %in% TXs) %>% 
    left_join(mapping) %>% 
    mutate(variant = factor(variant, levels = variants), cellline = factor(cellline, levels = celllines))
mu <- exprs2 %>% group_by(cellline, variant) %>% summarise(avgTPM = mean(TPM))

pdf(file='macroVariants_RNAseq.pdf', width=8, height=6)
ggplot() +
    geom_col(data = mu, aes(x = variant, y = avgTPM, fill = variant)) +
    geom_point(data = exprs2, aes(x = variant, y = TPM), color = "black") +
    facet_grid(cols = vars(cellline)) +
    # scale_y_continuous(trans='log10') +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


###################### Fig 2b (TCGA RNAseq data) #################

# 15808_enhancer_BRCA_normal_113.tsv file obtained from authors of Chen et al. Cell 2018. (A Pan-Cancer Analysis of Enhancer Expression in Nearly 9000 Patient Samples)
# enhancer_15808_BRCA_1095_RPKM.txt

# Get enhancers that intersect with classified CRE from HMEC
# tail -n +2 15808_enhancer_BRCA_normal_113.tsv | cut -f1 | sed 's/[:|-]/\t/g' > 15808_enhancers.bed
# bedtools intersect -loj -a HMEC-ATAC-K4m1-outBL.bed -b 15808_enhancers.bed > HMEC.TCGA.bed

myclasses <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")
color_key <- c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")

normal <- read.table("15808_enhancer_BRCA_normal_113.tsv", skip = 1)
normal.headers <- read.table("15808_enhancer_BRCA_normal_113.tsv", nrows = 1)

classmap <- read.table("HMEC_classes.bed") %>% 
    set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass")) %>% 
    select(peakname, newclass) %>% mutate(newclass = factor(newclass, levels = myclasses))
normal <- normal %>% 
    separate(V1, c("chr", "se"), sep = ":") %>% 
    separate(se, c("start", "end"), sep = "-", convert = TRUE) %>% 
    set_colnames(c("chr", "start", "end", normal.headers))

HMEC.TCGA <- read.table("HMEC.TCGA.bed") %>% filter(V6 != -1) %>% select(V5:V7,V4) %>% set_colnames(c("chr", "start", "end", "peakname"))
HMEC.TCGA.expr <- HMEC.TCGA %>% left_join(normal) %>% 
    pivot_longer(!c(chr, start, end, peakname), names_to="sample", values_to="RPKM") %>% 
    left_join(classmap)

logRPKM <- HMEC.TCGA.expr %>% 
    select(peakname, RPKM) %>% 
    group_by(peakname) %>% 
    summarize(meanRPKM = mean(RPKM)) %>% 
    left_join(classmap) %>% drop_na() %>% 
    mutate(logRPKM = log2(meanRPKM))
ntable <- logRPKM %>% group_by(newclass) %>% summarise(n = n())

p1 <- ggplot(logRPKM %>% mutate(logRPKM = 2^(logRPKM))) +
    geom_boxplot(mapping = aes(x=newclass, y=logRPKM, fill=newclass), outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # stat_summary(mapping = aes(x=newclass, y=logRPKM), fun.data = mean_sd, geom = "errorbar", lwd = 0.5, width = 0.25, color = "black") + #se
    # ggrastr::rasterise(geom_jitter(size = 0.3, stroke = 0, shape = 16, mapping = aes(x=newclass, y=logRPKM, fill=newclass), show.legend = FALSE), dpi = 300) +
    geom_text(data = ntable, mapping = aes(x = newclass, label = paste0(n)), y = -16, color = "black") +
    scale_fill_manual(values = color_key) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    ylab("RPKM") +
    scale_x_discrete(labels=c("Active" = "Active", "H3K4me3" = "APL", "ATAConly" = "ATAC-only", "Inactive" = "Inactive", "mBE" = "mBE")) +
    theme_cowplot()
    # ggtitle("TCGA Normal")

pdf(file='HMEC_expr_TCGA.pdf', width=5, height=5)
plot(p1 + theme(legend.position = "none", axis.title.x=element_blank()))
dev.off()

####################### Super Enhancer prediction of HMEC CRE using LILY ########################

# wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E119-H3K27ac.narrowPeak.gz
# wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E119-H3K27ac.broadPeak.gz
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E119-H3K27ac.fc.signal.bigwig
# mv E119-H3K27ac.narrowPeak E119-H3K27ac_peaks.narrowPeak
# mv E119-H3K27ac.broadPeak E119-H3K27ac_regions.bed
# mv E119-H3K27ac.fc.signal.bigwig E119-H3K27ac.bw
# cat /softwares/LILY/scripts/runLILY.R | R --slave --args E119-H3K27ac /outdir/ 12500 2500 /softwares/rose2/annotation/hg19_refseq.ucsc
# # (Note: I did some debugging/tweaks in runLILY.R to work around a bug in calculate_cutoff)
# # (Note: Disable .Rprofile before running LILY)
# grep SE E119-H3K27ac.scores.bed > HMEC.lily.se.bed
# bedtools intersect -loj -a HMEC-ATAC-K4m1-outBL.bed -b HMEC.lily.se.bed > HMEC.SE.bed

# wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E059-H3K27ac.narrowPeak.gz
# wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E059-H3K27ac.broadPeak.gz
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E059-H3K27ac.fc.signal.bigwig
# mv E059-H3K27ac.narrowPeak E059-H3K27ac_peaks.narrowPeak
# mv E059-H3K27ac.broadPeak E059-H3K27ac_regions.bed
# mv E059-H3K27ac.fc.signal.bigwig E059-H3K27ac.bw
# cat /softwares/LILY/scripts/runLILY.R | R --slave --args E059-H3K27ac /outdir/ 12500 2500 /softwares/rose2/annotation/hg19_refseq.ucsc
# grep SE E059-H3K27ac.scores.bed > NHM.lily.se.bed

# wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E118-H3K27ac.narrowPeak.gz
# wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E118-H3K27ac.broadPeak.gz
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E118-H3K27ac.fc.signal.bigwig
# mv E118-H3K27ac.narrowPeak E118-H3K27ac_peaks.narrowPeak
# mv E118-H3K27ac.broadPeak E118-H3K27ac_regions.bed
# mv E118-H3K27ac.fc.signal.bigwig E118-H3K27ac.bw
# cat /softwares/LILY/scripts/runLILY.R | R --slave --args E118-H3K27ac /outdir/ 12500 2500 /softwares/rose2/annotation/hg19_refseq.ucsc
# grep SE E118-H3K27ac.scores.bed > HepG2.lily.se.bed

# wget https://www.encodeproject.org/files/ENCFF855RCK/@@download/ENCFF855RCK.bam -O H3K27ac.bam -nv &
# wget https://www.encodeproject.org/files/ENCFF587UZE/@@download/ENCFF587UZE.bam -O H3K27ac_input.bam -nv &
# HMCan H3K27ac.bam H3K27ac_input.bam HMCan_config.txt MCF7-H3K27ac
# samtools faidx /references/fastas/hg19.fa
# wigToBigWig -clip MCF7-H3K27ac.wig /references/fastas/hg19.fa.fai MCF7-H3K27ac.bw
# grep SE MCF7-H3K27ac.scores.bed > MCF7.lily.se.bed

# intervene pairwise -i pairwise_1/* --filenames --compute frac --htype color --output pairwise_1/ --figsize 10 10

###################### Fig 2c (Expression of genes in HMEC - genes that are associated with classified CRE through GeneHancer link) #################

# Download RPKM quantification at Protein Coding regions, in Roadmap Reference Epigenomes
# wget https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.RPKM.pc.gz
# gunzip 57epigenomes.RPKM.pc.gz

# Process and clean gene names/ids in GeneHancer files
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz
# gunzip gencode.v19.annotation.gff3.gz
# grep -P "\tgene\t" gencode.v19.annotation.gff3 > gencode.v19.gene.gff3
# python process_gencodegff3.py &
# python clean_genehancer.py > GeneHancer.clean.bed &
# bedtools intersect -loj -a HMEC-ATAC-K4m1-outBL.bed -b GeneHancer.clean.bed > HMEC.gh.bed
# bedtools intersect -loj -a HMEC.lily.se.bed -b GeneHancer.clean.bed > HMEC.SE.gh.bed

myclasses <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")
classmap <- read.table("HMEC_classes.bed") %>% 
    set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass")) %>% 
    select(peakname, newclass) %>% mutate(newclass = factor(newclass, levels = myclasses))
expr_tbl <- read.table("57epigenomes.RPKM.pc", sep = "\t", skip = 1) %>% 
    select(V1, V53) %>% # V53 is E119 (HMEC)
    set_colnames(c("gene", "RPKM"))
genehancers <- read.table("HMEC.gh.bed") %>% 
    set_colnames(c("chr", "start", "end", "peakname", "ghchr", "ghstart", "ghend", "ghid", "gene", "score")) %>%
    filter(ghchr != ".") %>%
    left_join(classmap) %>%
    drop_na()
se <- read.table("HMEC.SE.gh.bed") %>% 
    select(V1:V3, V7:V12) %>%
    set_colnames(c("chr", "start", "end", "ghchr", "ghstart", "ghend", "ghid", "gene", "score")) %>% 
    filter(ghchr != ".")
aliases <- read.table("aliases", sep = "\t", comment.char = '', quote = '', na.strings = '', skip = 1) %>% 
    set_colnames(c("gene", "previous", "alias", "ensembl")) %>%
    drop_na(ensembl)

gene_tbl <- genehancers %>% 
    group_by(newclass, gene) %>% 
    dplyr::summarise(n = all()) %>% 
    filter(newclass != "H3K4me3") %>% # Removing class H3K4me3
    pivot_wider(names_from = newclass, values_from = n, values_fill = FALSE) %>%
    relocate(gene, Active, ATAConly, Inactive, mBE) %>%
    group_by(Active, ATAConly, Inactive, mBE)

gene_tbl %>% group_keys()
categories <- c("mBEonly", "InactiveOnly", "Combi", "ATAConly", "Combi", "Combi", "Combi", "ActiveCombi", "ActiveCombi", "ActiveCombi", "ActiveCombi", "ActiveCombi", "ActiveCombi", "ActiveCombi", "ActiveCombi")

category_levels=c("SE", "ActiveCombi", "Combi", "ATAConly", "InactiveOnly", "mBEonly")

gene_tbl <- data.frame(gene_tbl, category=recode_factor(gene_tbl %>% group_indices(), !!!categories)) %>%
    mutate(category = factor(category, levels=category_levels)) %>% select(gene, category)
gene_tbl <- bind_rows(gene_tbl, 
    se %>% select(gene) %>% unique() %>% mutate(category="SE"))

proc_tbl <- bind_rows(
        gene_tbl %>% filter(!str_detect(gene, "^ENSG")) %>% left_join(aliases) %>% mutate(gene = ensembl) %>% select(gene, category),
        gene_tbl %>% filter(str_detect(gene, "^ENSG"))) %>% 
    drop_na(gene) %>%
    left_join(expr_tbl) %>% 
    drop_na(RPKM) %>%
    mutate(category = factor(category, levels=category_levels))

my_comparisons <- list( c("SE", "ActiveCombi"), c("ActiveCombi", "InactiveOnly"), c("ActiveCombi", "mBEonly") )
gh_colors <- c(ActiveCombi = "#0f9448", ATAConly = "#E0AC69", mBEonly = "#2b598b", InactiveOnly = "#f15a2b", Combi = "#8B8680", SE = "#A5CF76")

ntable <- proc_tbl %>% group_by(category) %>% summarise(n = n())
pdf(file='HMEC_expr_GeneHancer.pdf', width=6, height=5)
ggplot(proc_tbl, aes(x = reorder(category, RPKM, FUN = median), y = RPKM, fill = category)) + 
    geom_boxplot(outlier.size = 0.7, outlier.stroke = 0, outlier.shape = 16, color = "black", lwd = 0.3, fatten = 2) +
    # ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16), dpi = 300) +
    geom_text(data = ntable, mapping = aes(x = category, label = paste0(n)), y = -11, color = "black") +
    scale_fill_manual(values = gh_colors) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))) +
    coord_cartesian(ylim = c(2^-11, 2^14)) +
    ylab("RPKM") +
    theme_cowplot() +
    theme(legend.position = "none", axis.title.x=element_blank())
dev.off()

#################### Fig 2d #########################

# mkdir FINAL_CLASSES/Active
# grep "Active" NHM/cCREs/classified.bed > FINAL_CLASSES/Active/A.bed
# grep "Active" HMEC1/cCREs/classified.bed > FINAL_CLASSES/Active/B.bed
# grep "Active" MCF7/cCREs/classified.bed > FINAL_CLASSES/Active/C.bed
# grep "Active" HepG2/cCREs/classified.bed > FINAL_CLASSES/Active/D.bed
# mkdir FINAL_CLASSES/Inactive
# grep "Inactive" NHM/cCREs/classified.bed > FINAL_CLASSES/Inactive/A.bed
# grep "Inactive" HMEC1/cCREs/classified.bed > FINAL_CLASSES/Inactive/B.bed
# grep "Inactive" MCF7/cCREs/classified.bed > FINAL_CLASSES/Inactive/C.bed
# grep "Inactive" HepG2/cCREs/classified.bed > FINAL_CLASSES/Inactive/D.bed
# mkdir FINAL_CLASSES/ATAConly
# grep "ATAConly" NHM/cCREs/classified.bed > FINAL_CLASSES/ATAConly/A.bed
# grep "ATAConly" HMEC1/cCREs/classified.bed > FINAL_CLASSES/ATAConly/B.bed
# grep "ATAConly" MCF7/cCREs/classified.bed > FINAL_CLASSES/ATAConly/C.bed
# grep "ATAConly" HepG2/cCREs/classified.bed > FINAL_CLASSES/ATAConly/D.bed

# intervene venn -i FINAL_CLASSES/Active/*.bed --names=NHM,HMEC,MCF7,HepG2 --output FINAL_CLASSES/Active/

##################### Fig 2e,f ######################

# bedtools intersect -a HMEC_mBE.bed -b MCF7_mBE.bed -u > HMEC_MCF7_mBE.bed
# bedtools intersect -a HMEC_APL.bed -b MCF7_APL.bed -u > HMEC_MCF7_APL.bed
# bedtools intersect -a HMEC_Active.bed -b MCF7_Active.bed -u > HMEC_MCF7_Active.bed
# bedtools intersect -a HMEC_Inactive.bed -b MCF7_Inactive.bed -u > HMEC_MCF7_Inactive.bed
# bedtools intersect -a HMEC_ATAConly.bed -b MCF7_ATAConly.bed -u > HMEC_MCF7_ATAConly.bed
# Run the resulting bed files in Cistrome-GO: http://go.cistrome.org/ and download results for KEGG pathways
# 1651531462_CistromeGO_go_kegg_result.txt mBE
# 1651531881_CistromeGO_go_kegg_result.txt APL
# 1651595932_CistromeGO_go_kegg_result.txt Active
# 1651596435_CistromeGO_go_kegg_result.txt Inactive
# 1651596597_CistromeGO_go_kegg_result.txt ATAConly

setwd("/mBE/")
files <- c("1651531462_CistromeGO_go_kegg_result.txt", "1651531881_CistromeGO_go_kegg_result.txt", "1651595932_CistromeGO_go_kegg_result.txt", "1651596435_CistromeGO_go_kegg_result.txt", "1651596597_CistromeGO_go_kegg_result.txt")
# files <- c("HMEC_mBE.txt", "HMEC_APL.txt", "HMEC_Active.txt", "HMEC_Inactive.txt", "HMEC_ATAConly.txt")
# files <- c("MCF7_mBE.txt", "MCF7_APL.txt", "MCF7_Active.txt", "MCF7_Inactive.txt", "MCF7_ATAConly.txt")
cistrome <- lapply(files, function(f){
    read.table(f, skip = 9, sep = "\t") %>% 
    select(V1, V5, V6) %>% set_colnames(c("term", "Padj", "count")) %>% 
    arrange(Padj) %>% head(n = 5) %>% 
    mutate(neg_log10_pval = -log10(Padj)) %>%
    separate(col = term, sep = "\\(", into = c("term", NA))
})
names(cistrome) <- c("mBE", "APL" ,"Active", "Inactive", "ATAConly")
cistrome <- bind_rows(cistrome, .id = "class")



setwd("/mBE/")

mBE <- read.table("1651531462_CistromeGO_go_kegg_result.txt", skip = 9, sep = "\t") %>% 
    select(V1, V5) %>% set_colnames(c("term", "Padj")) %>% 
    arrange(Padj) %>% head(n = 5) %>% 
    mutate(neg_log10_pval = -log10(Padj)) %>%
    separate(col = term, sep = "\\(", into = c("term", NA)) %>%
    rownames_to_column() %>% 
    mutate(rowname = factor(rowname, levels = rev(as.integer(rowname))))
APL <- read.table("1651531881_CistromeGO_go_kegg_result.txt", skip = 9, sep = "\t") %>% 
    select(V1, V5) %>% set_colnames(c("term", "Padj")) %>% 
    arrange(Padj) %>% head(n = 5) %>% 
    mutate(neg_log10_pval = -log10(Padj)) %>%
    separate(col = term, sep = "\\(", into = c("term", NA)) %>%
    rownames_to_column() %>% 
    mutate(rowname = factor(rowname, levels = rev(as.integer(rowname))))

p1 <- ggplot(mBE, aes(y = rowname, x = neg_log10_pval, fill = Padj < 0.05)) +
    geom_col(color = "black") + 
    scale_y_discrete(breaks=rownames(mBE), labels=mBE$term) +
    scale_fill_manual(values = c("grey", "#D55E00")) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    scale_x_continuous(expand = expansion(mult = c(0, .1))) +
    theme_cowplot() +
    theme(axis.title=element_blank(), axis.text=element_text(size = 9), axis.ticks = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
p2 <- ggplot(APL, aes(y = rowname, x = neg_log10_pval, fill = Padj < 0.05)) +
    geom_col(color = "black") + 
    scale_y_discrete(breaks=rownames(APL), labels=APL$term) +
    scale_fill_manual(values = c("grey", "#D55E00")) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    scale_x_continuous(expand = expansion(mult = c(0, .1))) +
    theme_cowplot() +
    theme(axis.title=element_blank(), axis.text=element_text(size = 9), axis.ticks = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

pdf(file='CistromeGO_KEGG.pdf', width=8, height=2)
p1 | p2
dev.off()


#################### Supp Fig 1a ####################

# Run classification algorithm separately (Run mBE_pipeline.sh as per documentation) 
# *_signals.txt obtained from classification algorithm 

# python clustering_stats.py HMEC_signals.txt > HMEC_cstats.txt &
# python clustering_stats.py NHM_signals.txt > NHM_cstats.txt &
# python clustering_stats.py MCF7_signals.txt > MCF7_cstats.txt &
# python clustering_stats.py 231L_signals.txt > 231L_cstats.txt &
# python clustering_stats.py HepG2_signals.txt > HepG2_cstats.txt &

cstat_files <- c("HMEC_cstats.txt", "NHM_cstats.txt", "MCF7_cstats.txt", "231L_cstats.txt", "HepG2_cstats.txt")
cell_lines <- c("HMEC", "NHM", "MCF7", "231L", "HepG2")

cstats <- lapply(cstat_files, read.table)
names(cstats) <- cell_lines
cstats <- bind_rows(cstats, .id = "cell.line") %>% 
    set_colnames(c("cell.line", "k", "avg_silhouette_score", "CH_score")) %>%
    mutate(k = factor(k, levels = as.character(2:11))) %>%
    mutate(cell.line = factor(cell.line, levels = cell_lines))

pdf(file='cstats.pdf', width=12, height=6)
p1 <- ggplot(cstats, aes(x = k, y = avg_silhouette_score, color = cell.line, group = cell.line)) +
    geom_point() +
    geom_line() +
    scale_x_discrete(breaks = as.character(2:11)) +
    theme_cowplot()
p2 <- ggplot(cstats, aes(x = k, y = CH_score, color = cell.line, group = cell.line)) +
    geom_point() +
    geom_line() +
    theme_cowplot()
plot(p1 | p2)
dev.off()

################### Supp Fig 1b #######################

# plotdata.rds from HMEC run of mBE_pipeline.sh classification algorithm
plotdata <- readRDS("plotdata.rds")
color_key <- c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")
p10 <- ggplot(plotdata %>% filter(ccre == "PLS")) +
    ggrastr::rasterise(geom_point(aes(x=mH2A2, y=H3K27ac, color=newclass), size = 0.5, stroke = 0, shape = 16), dpi = 300) +
    scale_colour_manual(values = color_key) +
    theme_cowplot() +
    theme(legend.position = "none") +
    ggtitle("PLS")

p11 <- ggplot(plotdata %>% filter(ccre == "dELS")) +
    ggrastr::rasterise(geom_point(aes(x=mH2A2, y=H3K27ac, color=newclass), size = 0.2, stroke = 0, shape = 16), dpi = 300) +
    scale_colour_manual(values = color_key) +
    theme_cowplot() +
    theme(legend.position = "none") +
    ggtitle("dELS")

p12 <- ggplot(plotdata %>% filter(ccre == "pELS")) +
    ggrastr::rasterise(geom_point(aes(x=mH2A2, y=H3K27ac, color=newclass), size = 0.4, stroke = 0, shape = 16), dpi = 300) +
    scale_colour_manual(values = color_key) +
    theme_cowplot() +
    theme(legend.position = "none") +
    ggtitle("pELS")

pdf(file='HMEC_PLS.pdf', width=3, height=3)
plot(p10)
dev.off()
pdf(file='HMEC_dELS.pdf', width=3, height=3)
plot(p11)
dev.off()
pdf(file='HMEC_pELS.pdf', width=3, height=3)
plot(p12)
dev.off()

################### Supp Fig 2b ####################

# bedtools intersect -loj -a HMEC-ATAC-K4m1-outBL.bed -b HMEC.lily.se.bed > peaks.SE.bed
ses <- read.table("peaks.SE.bed")
colnames(ses) <- c("chr", "start", "end", "peakname", "schr", "sstart", "send", "tag", "SEscore", "strand")
plotdata <- plotdata %>% left_join(ses %>% select(peakname, tag, SEscore))

seplot <- ggplot(plotdata) +
    ggrastr::rasterise(geom_point(data = subset(plotdata, tag != "SE"), aes(x=mH2A2, y=H3K27ac), size = 0.5, stroke = 0, shape = 16, color = "#cccccc"), dpi = 300) +
    new_scale_color() +
    ggrastr::rasterise(geom_point(data = subset(plotdata %>% arrange(SEscore), tag == "SE"), aes(x=mH2A2, y=H3K27ac, color=SEscore), size = 0.8, stroke = 0, shape = 16), dpi = 300) +
    scale_colour_gradient(low = "#b3e2cd", high = "black") + 
    # facet_wrap(vars(newclass)) +
    theme_cowplot()

pdf(file='HMEC_SE.pdf', width=4.5, height=3)
plot(seplot)
dev.off()

################## Supp Fig 2c and 2d (Build ChromHMM model and enrichment analysis) ###

# Download BAM files from Encode
# Make cellmarkfiletable
# Put BED files of each class from classification algorithm run on HMEC into /HMEC/classified/
# java -mx4000M -jar ${ChromHMM_INSTALLDIR}/ChromHMM.jar BinarizeBam ${ChromHMM_INSTALLDIR}/CHROMSIZES/hg19.txt ${BAM_DIR} HMEC.cellmarkfiletable /HMEC/BinarizedBam/
# java -mx4000M -jar ${ChromHMM_INSTALLDIR}/ChromHMM.jar LearnModel /HMEC/BinarizedBam/ /HMEC/out 13 hg19
# tail -n +2 /HMEC/out/HMEC_13_dense.bed | sort -k 1,1 -k2,2n > HMEC_13.bed
# java -mx4000M -jar ${ChromHMM_INSTALLDIR}/ChromHMM.jar OverlapEnrichment /HMEC/out/HMEC_13.bed /HMEC/classified/ /HMEC/out/classenrich

# tail -n +2 /HMEC/out/HMEC_13_dense.bed | sort -k 1,1 -k2,2n > HMEC_13.bed
# bedtools intersect -wao -a HMEC_classes.bed -b HMEC_13.bed > class_wao_CH13.bed

temp <- read.table("class_wao_CH13.bed") %>% select(V4, V7, V11, V2, V3, V9, V10, V17) %>% set_colnames(c("peak", "newclass", "chstate", "S1", "E1", "S2", "E2", "overlap"))
temp <- temp %>% group_by(peak) %>% arrange(S2)
temp <- temp %>% mutate(ovp = round(overlap * 100 / (E1 - S1)))

newclasses <- read.table("/HMEC/classified.bed") %>% select(V7) %>% set_colnames(c("newclass"))
newclass_counts <- data.frame(table(newclasses$newclass)) %>% set_colnames(c("newclass", "total"))

classes <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")
counts <- lapply(classes, function(class){
    return(
        temp %>% 
        filter(newclass == class) %>% 
        select(peak, chstate, ovp) %>% 
        uncount(ovp) %>% 
        group_by(peak) %>% 
        mutate(pos = row_number()) %>%
        mutate(chstate = factor(chstate, levels = 1:13), pos = factor(pos, levels = 1:105)) %>%
        group_by(pos, chstate) %>%
        tally())})
names(counts) <- classes
counts <- bind_rows(counts, .id = "newclass")
counts <- counts %>% left_join(newclass_counts) %>% mutate(prop = n/total, proptotal = n/sum(newclass_counts$total))
counts <- counts %>% mutate(newclass = factor(newclass, levels = c("ATAConly", "H3K4me3", "Active", "mBE", "Inactive")))

colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

pdf(file = "ChromHMM_overlap.pdf", width = 10, height = 5)
ggplot(counts %>% filter(pos == 50)) + 
    geom_bar(aes(fill=chstate, y=prop, x=newclass), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    scale_fill_manual(values = colors) +
    coord_flip() + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_cowplot() + 
    theme(axis.title.y = element_blank(), axis.title.x = element_blank())
dev.off()



















