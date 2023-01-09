# Run classification algorithm separately (Run mBE_pipeline.sh as per documentation) 
# *_classes.bed obtained from classification algorithm 

###########################################################################
# Collect all signals used in classification into one single table

setwd("/mforge/research/labs/experpath/maia/m237371/mBE/figures")
# ENCODE_cCREs_mm10.tab downloaded from UCSC Genome Browser's Table Browser https://genome.ucsc.edu/cgi-bin/hgTables
encode <- read.table("/mforge/research/labs/experpath/maia/m237371/public_data/ENCODE_cCREs_mm10.tab", skip=1, sep="\t") %>%
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

basedir <- "/mforge/research/labs/experpath/maia/m237371/mBE/DF/cCREs"
signal_files <- c("H3K27ac", "H2Az", "H3K27me3", "mH2A1", "mH2A2", "ATAC", "H3K4me1", "H3K4me3_signal", "CTCF")
signals <- c("H3K27ac", "H2Az", "H3K27me3", "mH2A1", "mH2A2", "ATAC", "H3K4me1", "H3K4me3", "CTCF")
loj_file <- "/mforge/research/labs/experpath/maia/m237371/mBE/DF/bedprep/DF_loj_cCREs.bed"
currpeaks_file <- "/mforge/research/labs/experpath/maia/m237371/mBE/DF/bedprep/DF-ATAC-outBL-resized.bed"
origpeaks_file <- "/mforge/research/labs/experpath/maia/m237371/mBE/DF/bedprep/DF-ATAC-outBL.bed"
DF <- get_signals()

setwd("/mforge/research/labs/experpath/maia/m237371/mBE/figures")
plotdata <- DF %>% filter(mark %in% c("H3K27ac", "H3K4me1", "H2Az", "H3K4me3", "H3K27me3", "mH2A1", "mH2A2", "CTCF")) %>%
    mutate(mark = factor(mark, levels = c("H3K27ac", "H3K4me1", "H2Az", "H3K4me3", "H3K27me3", "mH2A1", "mH2A2", "CTCF")))

plotdata1 <- plotdata %>% group_by(newclass, mark) %>% summarize(median_z = median(z)) %>%
    mutate(newclass=as.character(newclass)) %>%
    mutate(newclass=replace(newclass, newclass=="H3K4me3", "APL")) %>%
    mutate(newclass=replace(newclass, newclass=="ATAConly", "ATAC-only")) %>%
    mutate(newclass = factor(newclass, levels = rev(c("Active", "APL", "ATAC-only", "Inactive", "mBE")))) %>%
    mutate(mark = factor(mark, levels = c("H3K4me1", "H3K27ac", "H2Az", "H3K4me3", "H3K27me3", "CTCF", "mH2A1", "mH2A2")))

pdf(file='fig3b_fixed.pdf', width=4, height=3)
p1 <- ggplot(plotdata1 %>% filter(!(mark %in% c("H2Az", "CTCF"))), aes(x = mark, y = newclass)) +
    geom_tile(aes(fill = median_z)) +
    # geom_point(data = plotdata1 %>% filter(notsig == "*"), shape=1, size=1) +
    coord_equal() +
    scale_fill_gradientn(colors = bluered(256), limits = c(-3.2, 3.2)) +
    # scale_fill_gradient2(mid = "white", low = "blue", high = "red", limits = c(-2.05, 2.05)) +
    guides(fill=guide_colorbar(ticks.colour = NA)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), strip.background = element_blank())
p1
dev.off()

setwd("/data/")
write.table(plotdata1 %>% filter(!(mark %in% c("H2Az", "CTCF"))), file = "fig3a_data2.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

plotdata1 <- plotdata %>% filter(mark == "H3K4me1") %>% select(newclass) %>% table() %>% data.frame() %>%
    set_colnames(c("newclass", "Freq")) %>%
    mutate(newclass=as.character(newclass)) %>%
    mutate(newclass=replace(newclass, newclass=="H3K4me3", "APL")) %>%
    mutate(newclass=replace(newclass, newclass=="ATAConly", "ATAC-only")) %>%
    mutate(newclass = factor(newclass, levels = rev(c("Active", "APL", "ATAC-only", "Inactive", "mBE"))))

color_key <- c(Active = "#0f9448", APL = "#e78ac3", `ATAC-only` = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")
pdf(file='fig3b_bar.pdf', width=4, height=1)
p1 <- ggplot(plotdata1, aes(x = Freq, y = "y", fill = newclass, label = Freq)) +
    geom_bar(position="fill", stat="identity") +
    geom_text(size = 4, position = position_fill(vjust = 0.5), color = "white") +
    scale_x_continuous(expand = c(0,0), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), position = "top") + 
    scale_fill_manual(values = color_key) +
    theme_cowplot() + theme(legend.position = "None", axis.title.y = element_blank())
p1
dev.off()

setwd("/data/")
write.table(plotdata1, file = "fig3a_data1.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###################### Fig 3c Reprogramming factor (OSKM) binding sites enrichment #################

# Get TSS annotation for mm9 genome and define 1000 bp around TSS as a region to define "active TSS"
# git clone https://github.com/buenrostrolab/tss-annotation.git
mm9_tss <- read.table("tss-annotation/TSS/mm9.refGene.TSS.bed", sep="\t") %>%
    set_colnames(c("chr", "start", "end", "strand")) %>%
    mutate(start = start - 500, end = end + 500) %>%
    select(chr, start, end)
write.table(mm9_tss, file = "mm9_TSS.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Use RNAseq data generated in lab to determine active TSS
# STAR --runMode alignReads --runThreadN 16 --genomeDir /reference/STAR/mm9 --readFilesIn 1_dKO2_d3_S1_Lall_R1_001.fastq.bz2 --outFileNamePrefix 1_dKO2_ --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand bzcat
# samtools index 1_dKO2_Aligned.sortedByCoord.out.bam
# STAR --runMode alignReads --runThreadN 16 --genomeDir /reference/STAR/mm9 --readFilesIn 2_dKO2_d3_S1_Lall_R1_001.fastq.bz2 --outFileNamePrefix 2_dKO2_ --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand bzcat
# samtools index 2_dKO2_Aligned.sortedByCoord.out.bam
# bedtools multicov -bams 1_dKO2_Aligned.sortedByCoord.out.bam 2_dKO2_Aligned.sortedByCoord.out.bam -bed mm9_TSS.bed -S > DF.mm9_TSS.bed

# samtools view -F 0x904 -c 1_dKO2_Aligned.sortedByCoord.out.bam
# samtools view -F 0x904 -c 2_dKO2_Aligned.sortedByCoord.out.bam

DF_expr <- read.table("DF.mm9_TSS.bed") %>% 
    set_colnames(c("chr", "start", "end", "Rep1", "Rep2")) %>%
    mutate(length = end-start+1)

d <- DGEList(DF_expr %>% select(Rep1, Rep2))
d <- calcNormFactors(d)
RPKM <- rpkm(d, DF_expr$length, log=FALSE)
DF_expr <- DF_expr %>% select(-c(Rep1, Rep2)) %>% bind_cols(RPKM)
DF_expr <- DF_expr %>% mutate(meanExpr = rowMeans(across(c(Rep1, Rep2))))

pdf("DF_meanExpr.pdf", width=6, height=6)
ggplot(DF_expr, aes(x = "", y = meanExpr)) +
    geom_violin() +
    geom_hline(yintercept = 1) +
    scale_y_continuous(trans=pseudolog10_trans) +
    theme_cowplot()
dev.off()
# Define mean expression > 1 to be "active TSS"
write.table(DF_expr %>% filter(meanExpr > 1) %>% select(chr, start, end), file = "DF_activeTSS.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Get OSKM binding sites (ChIPseq peaks) from Chronis et al. Cell 2017. 
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2417nnn/GSM2417130/suppl/GSM2417130_48h_Oct4_OSKM.bed.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2417nnn/GSM2417131/suppl/GSM2417131_48h_Sox2_OSKM.bed.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2417nnn/GSM2417132/suppl/GSM2417132_48h_Klf4_OSKM.bed.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2417nnn/GSM2417133/suppl/GSM2417133_48h_cMyc_OSKM.bed.gz
# gunzip *OSKM.bed.gz
# sort -k 1,1 -k2,2n GSM2417133_48h_cMyc_OSKM.bed | awk '{print $1, $2, $3, "Peak", "0", "+"}' | awk '{$1=$1}1' OFS="\t" > Myc_6col.bed
# sort -k 1,1 -k2,2n GSM2417132_48h_Klf4_OSKM.bed | awk '{print $1, $2, $3, "Peak", "0", "+"}' | awk '{$1=$1}1' OFS="\t" > Klf4_6col.bed
# sort -k 1,1 -k2,2n GSM2417130_48h_Oct4_OSKM.bed | awk '{print $1, $2, $3, "Peak", "0", "+"}' | awk '{$1=$1}1' OFS="\t" > Oct4_6col.bed
# sort -k 1,1 -k2,2n GSM2417131_48h_Sox2_OSKM.bed | awk '{print $1, $2, $3, "Peak", "0", "+"}' | awk '{$1=$1}1' OFS="\t" > Sox2_6col.bed

# Get ChIPseq signals for histone mark H3K27me3 and histone variants mH2A1 and mH2A2 from GSE40813
# and get the signal intensity at the OSKM binding sites
# computeMatrix scale-regions -S DF_H3K27me3_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R DF_activeTSS.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix actTSS_H3K27me3 \
                            # --outFileName actTSS_H3K27me3.tab.gz
# computeMatrix scale-regions -S DF_H3K27me3_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Oct4_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Oct4_H3K27me3 \
                            # --outFileName Oct4_H3K27me3.tab.gz
# computeMatrix scale-regions -S DF_H3K27me3_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Sox2_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Sox2_H3K27me3 \
                            # --outFileName Sox2_H3K27me3.tab.gz
# computeMatrix scale-regions -S DF_H3K27me3_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Klf4_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Klf4_H3K27me3 \
                            # --outFileName Klf4_H3K27me3.tab.gz
# computeMatrix scale-regions -S DF_H3K27me3_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Myc_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Myc_H3K27me3 \
                            # --outFileName Myc_H3K27me3.tab.gz

# computeMatrix scale-regions -S DF_mH2A1_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R DF_activeTSS.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix actTSS_mH2A1 \
                            # --outFileName actTSS_mH2A1.tab.gz
# computeMatrix scale-regions -S DF_mH2A1_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Oct4_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Oct4_mH2A1 \
                            # --outFileName Oct4_mH2A1.tab.gz
# computeMatrix scale-regions -S DF_mH2A1_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Sox2_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Sox2_mH2A1 \
                            # --outFileName Sox2_mH2A1.tab.gz
# computeMatrix scale-regions -S DF_mH2A1_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Klf4_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Klf4_mH2A1 \
                            # --outFileName Klf4_mH2A1.tab.gz
# computeMatrix scale-regions -S DF_mH2A1_wt_m.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Myc_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Myc_mH2A1 \
                            # --outFileName Myc_mH2A1.tab.gz

# computeMatrix scale-regions -S DF_mH2A2_wt_f.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R DF_activeTSS.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix actTSS_mH2A2 \
                            # --outFileName actTSS_mH2A2.tab.gz
# computeMatrix scale-regions -S DF_mH2A2_wt_f.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Oct4_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Oct4_mH2A2 \
                            # --outFileName Oct4_mH2A2.tab.gz
# computeMatrix scale-regions -S DF_mH2A2_wt_f.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Sox2_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Sox2_mH2A2 \
                            # --outFileName Sox2_mH2A2.tab.gz
# computeMatrix scale-regions -S DF_mH2A2_wt_f.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Klf4_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Klf4_mH2A2 \
                            # --outFileName Klf4_mH2A2.tab.gz
# computeMatrix scale-regions -S DF_mH2A2_wt_f.ChIP_MACS2.mm9.q1e-2_FE.bw \
                            # -R Myc_6col.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix Myc_mH2A2 \
                            # --outFileName Myc_mH2A2.tab.gz

# Test enrichment and plot the data
TFs <- c("actTSS", "Myc", "Klf4", "Oct4", "Sox2")
marks <- c("H3K27me3", "mH2A1", "mH2A2")

read_mark1 <- function(mark){
    tab <- read.table(gzfile(mark), skip=1, sep="\t")
    return(tab %>% select(V7))
}

actTSS <- bind_cols(lapply(paste0("actTSS_", marks, ".tab.gz"), read_mark1)) %>% set_colnames(marks) %>% mutate(TF = "actTSS")
Myc <- bind_cols(lapply(paste0("Myc_", marks, ".tab.gz"), read_mark1)) %>% set_colnames(marks) %>% mutate(TF = "Myc")
Klf4 <- bind_cols(lapply(paste0("Klf4_", marks, ".tab.gz"), read_mark1)) %>% set_colnames(marks) %>% mutate(TF = "Klf4")
Oct4 <- bind_cols(lapply(paste0("Oct4_", marks, ".tab.gz"), read_mark1)) %>% set_colnames(marks) %>% mutate(TF = "Oct4")
Sox2 <- bind_cols(lapply(paste0("Sox2_", marks, ".tab.gz"), read_mark1)) %>% set_colnames(marks) %>% mutate(TF = "Sox2")
dat <- bind_rows(actTSS, Myc, Klf4, Oct4, Sox2) %>%
    pivot_longer(H3K27me3:mH2A2, names_to = "mark", values_to = "avg_signal") %>%
    mutate(TF = factor(TF, levels = TFs), mark = factor(mark, levels = marks))

dim(actTSS)
dim(Oct4)
dim(Sox2)
dim(Klf4)
dim(Myc)

nvalues <- data.frame(peakset = c("actTSS", "Oct4", "Sox2", "Klf4", "Myc"), n = c(nrow(actTSS), nrow(Oct4), nrow(Sox2), nrow(Klf4), nrow(Myc)), mark = "H3K27me3")
my_comparisons <- list( c("actTSS", "Oct4"), c("actTSS", "Sox2"), c("actTSS", "Klf4"), c("actTSS", "Myc") )
p1 <- ggplot(dat %>% mutate(TF = factor(TF, levels = rev(c("actTSS", "Oct4", "Sox2", "Klf4", "Myc")))), aes(y = TF, x = avg_signal)) +
    geom_boxplot(aes(fill = TF), outlier.shape = NA) +
    geom_text(data = nvalues, mapping = aes(y = peakset, label = paste0("n = ",n)), x = 1.8, color = "black", size = 3) +
    coord_cartesian(xlim=c(0,2.1)) +
    # stat_compare_means(comparisons = my_comparisons, tip.length = 0) +
    # scale_x_continuous(trans = pseudolog10_trans) + 
    scale_fill_manual(values = rev(c("white", brewer.pal(4, "Set2")))) +
    facet_grid(rows=vars(mark)) +
    theme_cowplot() + theme(legend.position = "none")

setwd("/data/")
write.table(dat, file = "suppfig3a_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

TFs <- c("Oct4", "Sox2", "Klf4", "Myc")
marks <- c("H3K27me3", "mH2A1", "mH2A2")
ins <- expand.grid(TFs, marks)
ins <- split(ins, seq(nrow(ins)))
get_tests <- function(temp){
    x <- dat %>% filter(TF == "actTSS" & mark == as.character(temp$Var2)) %>% select(avg_signal) %>% unlist() %>% unname()
    y <- dat %>% filter(TF == as.character(temp$Var1) & mark == as.character(temp$Var2)) %>% select(avg_signal) %>% unlist() %>% unname()
    p <- wilcox.test(x = x, y = y)$p.value
    return(data.frame(p = p, FC = median(y) / median(x)))
}
tests <- lapply(ins, get_tests)
names(tests) <- lapply(ins, function(temp){return(paste0(as.character(temp$Var1), "_", as.character(temp$Var2)))})
tests <- bind_rows(tests, .id = "TF_mark")
tests <- tests %>% separate(TF_mark, sep = "_", into = c("TF", "mark"))
tests <- tests %>% mutate(p = case_when(p == 0 ~ 1e-300, TRUE ~ p)) %>%
    mutate(neg_log10_pvalue = -log10(p), log2_FC = log2(FC))
tests <- tests %>% mutate(TF = factor(TF, levels = TFs))
# jittered_p <- tests %>% mutate(p = p + (10^(-1 * runif(12, min=235, max=245)))) %>%
    # mutate(neg_log10_pvalue = -log10(p))

p2 <- ggplot(tests, aes(x = log2_FC, y = neg_log10_pvalue)) +
    geom_point(aes(color = TF, shape = mark), size = 4) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    scale_color_brewer(palette = "Set2") +
    # scale_shape_manual(values=c(0, 8, 2)) +
    # coord_cartesian(xlim=c(-0.5, 1.2)) +
    theme_cowplot()
pdf("DF_ChIPSeq_OSKM_volcano.pdf", width=10, height=4)
plot(p1 + p2 + plot_layout(widths = c(3, 2)))
dev.off()

setwd("/data/")
write.table(tests, file = "fig3c_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


########### Fig 3d (GAT enrichment analysis of ChIPseq peaks from Chronis et al. Cell 2017 against classified CRE in DF ###############
## cell_9383_mmc1_48h.bed obtained from Chronis et al. Cell 2017
# grep 48h_OSKM_Brg1 cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_Brg1.bed
# grep 48h_OSKM_Cebpa cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_Cebpa.bed
# grep 48h_OSKM_Cebpb cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_Cebpb.bed
# grep 48h_OSKM_Fra1 cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_Fra1.bed
# grep 48h_OSKM_Hdac1 cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_Hdac1.bed
# grep 48h_OSKM_Klf4 cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_Klf4.bed
# grep 48h_OSKM_Oct4 cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_Oct4.bed
# grep 48h_OSKM_Runx1 cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_Runx1.bed
# grep 48h_OSKM_Sox2 cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_Sox2.bed
# grep 48h_OSKM_cMyc cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_cMyc.bed
# grep 48h_OSKM_p300 cell_9383_mmc1_48h.bed | sort -k1,1 -k2,2n > 48h_OSKM_p300.bed

# sort -k1,1 -k2,2n cell_9383_mmc1_48h.bed > cell_9383_mmc1_48h_sorted.bed
# bedtools merge -i cell_9383_mmc1_48h_sorted.bed > cell_9383_mmc1_48h_merged.bed

# cat <(echo 'track name="Brg1"') 48h_OSKM_Brg1.bed <(echo 'track name="Cebpa"') 48h_OSKM_Cebpa.bed <(echo 'track name="Cebpb"') 48h_OSKM_Cebpb.bed <(echo 'track name="Fra1"') 48h_OSKM_Fra1.bed <(echo 'track name="Hdac1"') 48h_OSKM_Hdac1.bed <(echo 'track name="Klf4"') 48h_OSKM_Klf4.bed <(echo 'track name="Oct4"') 48h_OSKM_Oct4.bed <(echo 'track name="Runx1"') 48h_OSKM_Runx1.bed <(echo 'track name="Sox2"') 48h_OSKM_Sox2.bed <(echo 'track name="cMyc"') 48h_OSKM_cMyc.bed <(echo 'track name="p300"') 48h_OSKM_p300.bed > DF_TFBS_tracks.bed

# gat-run.py --with-segment-tracks --segments=DF_TFBS_tracks.bed --annotations=DF_classified_4col.bed --workspace=cell_9383_mmc1_48h_merged.bed --num-samples=1000 --log=1 1> DF_TFBS_tracks.GAT 2> err &

setwd("/mforge/research/labs/experpath/maia/m237371/mBE/GAT")
color_key <- c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", Inactive = "#f15a2b", mBE = "#2b598b")

out <- read.table("DF_TFBS_tracks.GAT", header = TRUE, na.strings = "na")
out <- out %>% mutate(track = factor(track, levels = out %>% filter(annotation == "mBE") %>% select(track, l2fold) %>% arrange(l2fold) %>% select(track) %>% unlist() %>% unname()), annotation = factor(annotation, levels = c("mBE", "Inactive", "ATAConly", "Active", "H3K4me3")))
pdf(file='gat_DF_TFBS_hm_temp.pdf', width=5, height=2)
p1 <- ggplot(out, aes(x = track, y = annotation)) +
    geom_tile(aes(fill = l2fold)) +
    coord_equal() +
    # scale_fill_gradientn(colors = bluered(256)) +
    scale_fill_gradient2(low = bluered(256)[[1]], mid = "white", high = bluered(256)[[256]], midpoint = 0) +
    scale_x_discrete(position = "top") +
    guides(fill=guide_colorbar(ticks.colour = NA)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
p1
dev.off()

out %<>% mutate(signif = qvalue < 0.05) %>% select(track, annotation, l2fold, signif)
setwd("/data/")
write.table(out, file = "fig3d_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


########### Fig 4c (PCA of H3K37ac signals from breast cancer subtypes) #######

# Download H3K27ac ChIPseq signals (bigwig) from GSE85158 for the 12 breast cancer types

# grep "mBE" HMEC_classes.bed > mBE.bed

# multiBigwigSummary BED-file -b MCF10A.H3K27ac.bw MCF7.H3K27ac.bw ZR751.H3K27ac.bw MB361.H3K27ac.bw UACC812.H3K27ac.bw SKBR3.H3K27ac.bw AU565.H3K27ac.bw HCC1954.H3K27ac.bw MB231.H3K27ac.bw MB436.H3K27ac.bw MB468.H3K27ac.bw HCC1937.H3K27ac.bw -o mBE.npz --BED mBE.bed &
# multiBigwigSummary BED-file -b MCF10A.H3K27ac.bw MCF7.H3K27ac.bw ZR751.H3K27ac.bw MB361.H3K27ac.bw UACC812.H3K27ac.bw SKBR3.H3K27ac.bw AU565.H3K27ac.bw HCC1954.H3K27ac.bw MB231.H3K27ac.bw MB436.H3K27ac.bw MB468.H3K27ac.bw HCC1937.H3K27ac.bw -o ccres.npz --BED HMEC_classes.bed &
# plotPCA -in mBE.npz -o mBE_PCA.pdf -l MCF10A MCF7 ZR751 MB361 UACC812 SKBR3 AU565 HCC1954 MB231 MB436 MB468 HCC1937 --outFileNameData mBE_PCA.tab &
# plotPCA -in ccres.npz -o ccres_PCA.pdf -l MCF10A MCF7 ZR751 MB361 UACC812 SKBR3 AU565 HCC1954 MB231 MB436 MB468 HCC1937 --outFileNameData ccre_PCA.tab &

################ GARFIELD (Fig 4d Enrichment of GWAS variants in classified CRE in HMEC) ###################

# cd /softwares/
# wget https://www.ebi.ac.uk/birney-srv/GARFIELD/package-v2/garfield-v2.tar.gz
# wget -bqc https://www.ebi.ac.uk/birney-srv/GARFIELD/package-v2/garfield-data.tar.gz -o garfield-data.tar.gz &

# ## Create pval data from summary statistics (See documentation)
# wget https://bcac.ccge.medschl.cam.ac.uk/files/oncoarray_bcac_public_release_oct17.txt.gz
# vi /softwares/garfield-v2/garfield-create-input-gwas.sh
# # chrcol=3, poscol=4, pvalcol=10
# # TRAITNAME=BC_GCST004988
# # GWASFILENAME=oncoarray_bcac_public_release_oct17.txt
# /softwares/garfield-v2/garfield-create-input-gwas.sh

# ## Create annotation data encoded in format described in documentation
# mkdir /temp/UK10K_variants 
# for i in {1..22} 'X'
# do
    # awk -v i=$i '{print "chr"i"\t"$1-1"\t"$1}' /softwares/garfield-data/maftssd/chr${i} > /temp/UK10K_variants/chr${i}.bed
# done

# cat /temp/UK10K_variants/chr{1..22}.bed /temp/UK10K_variants/chrX.bed >> /temp/UK10K_variants.bed

# bedtools intersect -loj -a /temp/UK10K_variants.bed -b /temp/HMEC_classified.bed > /temp/HMEC_annotations.bed
# awk '{if ($10 == ".") print $1" "$3" 00000"; else if ($10 == "Active") print $1" "$3" 10000"; else if ($10 == "H3K4me3") print $1" "$3" 01000"; else if ($10 == "ATAConly") print $1" "$3" 00100"; else if ($10 == "Inactive") print $1" "$3" 00010"; else if ($10 == "mBE") print $1" "$3" 00001"; else print "ERROR"} ' /temp/HMEC_annotations.bed > /temp/HMEC_annotations.encoded

# bedtools intersect -loj -a /temp/UK10K_variants.bed -b /temp/MCF7_classified.bed > /temp/MCF7_annotations.bed
# awk '{if ($10 == ".") print $1" "$3" 00000"; else if ($10 == "Active") print $1" "$3" 10000"; else if ($10 == "H3K4me3") print $1" "$3" 01000"; else if ($10 == "ATAConly") print $1" "$3" 00100"; else if ($10 == "Inactive") print $1" "$3" 00010"; else if ($10 == "mBE") print $1" "$3" 00001"; else print "ERROR"} ' /temp/MCF7_annotations.bed > /temp/MCF7_annotations.encoded

# bedtools intersect -loj -a /temp/UK10K_variants.bed -b /temp/231L_classified.bed > /temp/231L_annotations.bed
# awk '{if ($10 == ".") print $1" "$3" 0000"; else if ($10 == "Active") print $1" "$3" 1000"; else if ($10 == "H3K4me3") print $1" "$3" 0100"; else if ($10 == "ATAConly") print $1" "$3" 0010"; else if ($10 == "Inactive") print $1" "$3" 0001"; else print "ERROR"} ' /temp/231L_annotations.bed > /temp/231L_annotations.encoded

# paste -d" " /temp/HMEC_annotations.encoded /temp/MCF7_annotations.encoded /temp/231L_annotations.encoded | awk -F" " '{print $1" "$2" "$3 $6 $9}' > /temp/breast_annotations.encoded

# for i in {1..22} 'X'
# do
    # grep "chr${i} " /temp/breast_annotations.encoded | cut -d" " -f2- > /softwares/garfield-data/annotation-breast/chr${i}
# done
# vi /softwares/garfield-data/annotation-breast/link_file.txt # Create link_file (see documentation)

# vi /softwares/garfield-v2/garfield
# # INPUTNAME=BC_GCST004988
# # DATADIR=/softwares/garfield-data
# /softwares/garfield-v2/garfield
# # Output in /softwares/garfield-data/output/BC_GCST004988/garfield.test.BC_GCST004988.out

setwd("/softwares/garfield-data/output/BC_GCST004988/")
color_key <- c(Active = "#0f9448", APL = "#e78ac3", ATAConly = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")

tab1 <- read.table("garfield.test.BC_GCST004988.out", header = TRUE) %>% 
    separate(col = "Annotation", sep = "_", into = c(NA, "Annotation")) %>%
    mutate(Annotation = factor(Annotation, levels = names(color_key))) %>% 
    mutate(Celltype = factor(Celltype, levels = c("HMEC", "MCF7", "231L"))) %>% 
    mutate(PThresh = factor(PThresh))

pdf(file = "garfield.test.BC_GCST004988.pdf", width = 5, height = 4)
p1 <- ggplot(tab1 %>% filter(PThresh == "5e-08"), aes(x = log2(OR), y = -log10(Pvalue), shape = Celltype, color = Annotation)) +
    geom_point(size = 2) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed") +
    scale_shape_manual(values=c(16, 4, 15)) +
    scale_color_manual(values=color_key) +
    # lims(x = c(-3, 3)) +
    theme_cowplot()
p1
dev.off()

tab1 %<>% filter(PThresh == "5e-08") %>% mutate(log2OR = log2(OR), neg_log10_pval = -log10(Pvalue)) %>% select(c(Celltype, Annotation, log2OR, neg_log10_pval))
setwd("/data/")
write.table(tab1, file = "fig4d_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)






