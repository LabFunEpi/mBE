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

########### Fig 3d (Fisher's exact test - by intersections of ChIPseq peaks from Chronis et al. Cell 2017 against classified CRE in DF ###############
# fisherTestIntersect1c.py --tfBed cell_9383_mmc1_48h.bed --enhancerBed classified_wNoClass.short.bed --sortBy mBE --plotOrder mBE Inactive ATAConly Active H3K4me3 --trimTFNames 48h_OSKM_

########### Fig 4c (PCA of H3K37ac signals from breast cancer subtypes) #######

# Download H3K27ac ChIPseq signals (bigwig) from GSE85158 for the 12 breast cancer types

# grep "mBE" HMEC_classes.bed > mBE.bed

# multiBigwigSummary BED-file -b MCF10A.H3K27ac.bw MCF7.H3K27ac.bw ZR751.H3K27ac.bw MB361.H3K27ac.bw UACC812.H3K27ac.bw SKBR3.H3K27ac.bw AU565.H3K27ac.bw HCC1954.H3K27ac.bw MB231.H3K27ac.bw MB436.H3K27ac.bw MB468.H3K27ac.bw HCC1937.H3K27ac.bw -o mBE.npz --BED mBE.bed &
# multiBigwigSummary BED-file -b MCF10A.H3K27ac.bw MCF7.H3K27ac.bw ZR751.H3K27ac.bw MB361.H3K27ac.bw UACC812.H3K27ac.bw SKBR3.H3K27ac.bw AU565.H3K27ac.bw HCC1954.H3K27ac.bw MB231.H3K27ac.bw MB436.H3K27ac.bw MB468.H3K27ac.bw HCC1937.H3K27ac.bw -o ccres.npz --BED HMEC_classes.bed &
# plotPCA -in mBE.npz -o mBE_PCA.pdf -l MCF10A MCF7 ZR751 MB361 UACC812 SKBR3 AU565 HCC1954 MB231 MB436 MB468 HCC1937 --outFileNameData mBE_PCA.tab &
# plotPCA -in ccres.npz -o ccres_PCA.pdf -l MCF10A MCF7 ZR751 MB361 UACC812 SKBR3 AU565 HCC1954 MB231 MB436 MB468 HCC1937 --outFileNameData ccre_PCA.tab &

############ Fig 4d (Fisher's exact test of presence of breast cancer GWAS variants in classified CRE in HMEC) ##########

# Download GWAS data for breast cancer: NHGRI-EBI GWAS Catalog: https://www.ebi.ac.uk/gwas/home; EFO ID: EFO_0000305
# We use data downloaded on 2021-06-14

library(splitstackshape) 
gwas <- read.table("gwas-association-downloaded_2022-04-11-EFO_0006861-withChildTraits.tsv", skip = 1, sep = "\t", quote = "") %>% 
    select(V22) %>% 
    cSplit('V22', ';') %>%
    unlist() %>% unname()
gwas <- unique(gwas[!is.na(gwas)])
setdiff(gwas, union(str_subset(gwas, "^rs"), str_subset(gwas, "^chr")))
write.table(as.data.frame(str_subset(gwas, "^rs")), file = "gwas_rs.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

gwas_chr <- data.frame(x = str_subset(gwas, "^chr")) %>% 
    separate(x, c("chr", "start")) %>% 
    mutate(start = as.numeric(start)) %>%
    mutate(start = start - 1) %>%
    mutate(end = start + 1) %>% 
    unite("name", c(chr, end), sep = "-", remove = FALSE) %>% 
    relocate(name, .after = end)
write.table(gwas_chr, file = "gwas_chr.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Lift Over to hg19
# https://genome.ucsc.edu/cgi-bin/hgLiftOver
# lifted_over.bed

# Download variant database to get the positions for all variants with variant id
# http://grch37.ensembl.org/biomart/martview/7d7265ec7122686b3585e4b3e33f6bee
# mart_export.txt

gwas_rs <- read.table("mart_export.txt", skip = 1, sep = "\t")
gwas_rs <- gwas_rs %>% filter(!grepl("^H", V3)) %>% 
    select(-V2) %>% set_colnames(c("name", "chr", "start", "end")) %>% 
    mutate(chr = paste0("chr", chr)) %>% relocate(name, .after = end)

gwas_rs_clean <- gwas_rs %>% filter(start > end) %>% relocate(start, .after = "end") %>% set_colnames(c("chr", "start", "end", "name")) %>%
    bind_rows(gwas_rs %>% filter(start == end) %>% mutate(start = start - 1)) %>%
    bind_rows(gwas_rs %>% filter(start < end))

gwas_chr <- read.table("gwas2022_liftedOver.bed", sep = "\t") %>%
    set_colnames(c("chr", "start", "end", "name"))

gwas <- bind_rows(gwas_rs_clean, gwas_chr)
write.table(gwas, file = "gwas1_hg19.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# python fisherTestGWAS.py
# python fisherTestGWASByClass.py --allEnhancers --cellLine "HMEC" --pCrit 0.001
# python fisherTestGWASByClass.py --allEnhancers --cellLine "MCF7" --pCrit 0.001
# python fisherTestGWASByClass.py --allEnhancers --cellLine "231L" --pCrit 0.001
# python plotByGWASVolcano.py --allEnhancers --plotRange 3.5 --pCrit 0.001 --yRange 15

























