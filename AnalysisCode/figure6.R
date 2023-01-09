
wd <- "/mBE/231L/"
setwd(wd)
signal_files <- c("mH2A2", "H3K27ac", "H3K4me1", "H3K27me3", "mH2A1", "ATAC", "m2mH2A2", "m2H3K27ac", "m2H3K4me1", "m2H3K27me3", "m2ATAC", "p300", "Input", "BRD4", "m2p300", "m2Input", "m2BRD4", "BRD4S1", "BRD4S2", "BRD4L")
signals <- c("mH2A2", "H3K27ac", "H3K4me1", "H3K27me3", "mH2A1", "ATAC", "m2mH2A2", "m2H3K27ac", "m2H3K4me1", "m2H3K27me3", "m2ATAC", "p300", "Input", "BRD4", "m2p300", "m2Input", "m2BRD4", "BRD4S1", "BRD4S2", "BRD4L")
rawdatas <- lapply(signal_files, read.table, skip=3)
rawdatas <- lapply(rawdatas, function(x){x %>% mutate(across(everything(), ~replace_na(.x, 0)))})
print(lapply(rawdatas, dim))

# Transformations
datas <- rawdatas
mydata <- bind_cols(lapply(datas, function(x){x %>% rowSums(.)}))
colnames(mydata) <- signals

encode <- read.table("/public_data/ENCODE_cCREs_hg38-1.tab", skip=1, sep="\t") %>%
    select(V4, V11) %>%
    set_colnames(c("name", "ccre"))

hg19list <- read.table("/public_data/ENCODE_cCREs_hg19.bed") %>%
    select(V1:V4) %>%
    set_colnames(c("echr", "estart", "eend", "name")) %>%
    left_join(encode)

lojed <- read.table("/mBE/231L/candidates_wENCODE.bed") %>%
    select(V4,V10) %>%
    set_colnames(c("peakname", "name")) %>%
    left_join(encode) %>%
    group_by(peakname) %>% dplyr::slice(1) %>%
    replace_na(list(ccre = "notccre"))

currpeaks <- read.table("/mBE/231L/candidate_enhancer_centers.bed") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname")) %>% 
    left_join(lojed)
origpeaks <- read.table("/mBE/231L/candidate_enhancers.bed") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname"))

mydata <- mydata %>% bind_cols(currpeaks)

mydata1 <- mydata %>% filter(ccre == "notccre")
mydata2 <- mydata %>% filter(ccre != "notccre")

plotdata <- mydata2
myclasses <- c("Active", "H3K4me3", "ATAConly", "Inactive")
color_key <- c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", Inactive = "#f15a2b")

classified <- read.table("classified.bed") %>% set_colnames(c("chr", "start", "end", "peakname", "name", "ccre", "newclass"))
plotdata <- plotdata %>% bind_cols(classified %>% select(newclass)) %>% mutate(newclass=factor(newclass, levels = myclasses)) %>% select(-c(chr, start, end)) %>% left_join(origpeaks)

###################### Fig 6a ######################

  # Active  H3K4me3 ATAConly Inactive
   # 16341     5194    13921    14875
library(ggpmisc)
theme_set(theme_cowplot())
plotdata1 <- plotdata
p1 <- ggplot(data = plotdata1, mapping = aes(x=ATAC, y=m2ATAC)) +
    ggrastr::rasterise(geom_point(mapping = aes(color=newclass), size = 0.5, stroke = 0, shape = 16), dpi = 400) +
    geom_abline(intercept = 0, slope = 1, linetype=3) +
    # stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + 
    stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), geom = 'text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), label.y = "bottom", label.x = 20000, size = 3) +
    geom_smooth(method='lm', se=FALSE, color="black", fullrange=TRUE, lwd=0.5) +
    scale_colour_manual(values = color_key) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    coord_cartesian(xlim = c(0, 35000), ylim = c(0, 35000)) + 
    theme(legend.position = "none")
p2 <- ggplot(data = plotdata1, mapping = aes(x=mH2A2, y=m2mH2A2)) +
    ggrastr::rasterise(geom_point(mapping = aes(color=newclass), size = 0.5, stroke = 0, shape = 16), dpi = 400) +
    geom_abline(intercept = 0, slope = 1, linetype=3) +
    # stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + 
    stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), geom = 'text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), label.y = "bottom", label.x = 2000, size = 3) +
    geom_smooth(method='lm', se=FALSE, color="black", fullrange=TRUE, lwd=0.5) +
    scale_colour_manual(values = color_key) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    coord_cartesian(xlim = c(0, 5000), ylim = c(0, 5000)) + 
    theme(legend.position = "none")
p3 <- ggplot(data = plotdata1, mapping = aes(x=H3K4me1, y=m2H3K4me1)) +
    ggrastr::rasterise(geom_point(mapping = aes(color=newclass), size = 0.5, stroke = 0, shape = 16), dpi = 400) +
    geom_abline(intercept = 0, slope = 1, linetype=3) +
    # stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + 
    stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), geom = 'text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), label.y = "bottom", label.x = 'right', size = 3) +
    geom_smooth(method='lm', se=FALSE, color="black", fullrange=TRUE, lwd=0.5) +
    scale_colour_manual(values = color_key) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    coord_cartesian(xlim = c(0, 700), ylim = c(0, 700)) + 
    theme(legend.position = "none")
p4 <- ggplot(data = plotdata1, mapping = aes(x=H3K27ac, y=m2H3K27ac)) +
    ggrastr::rasterise(geom_point(mapping = aes(color=newclass), size = 0.5, stroke = 0, shape = 16), dpi = 400) +
    geom_abline(intercept = 0, slope = 1, linetype=3) +
    # stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + 
    stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), geom = 'text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), label.y = "bottom", label.x = 'right', size = 3) +
    geom_smooth(method='lm', se=FALSE, color="black", fullrange=TRUE, lwd=0.5) +
    scale_colour_manual(values = color_key) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    coord_cartesian(xlim = c(0, 820), ylim = c(0, 820)) + 
    theme(legend.position = "none")
p5 <- ggplot(data = plotdata1, mapping = aes(x=p300, y=m2p300)) +
    ggrastr::rasterise(geom_point(mapping = aes(color=newclass), size = 0.5, stroke = 0, shape = 16), dpi = 400) +
    geom_abline(intercept = 0, slope = 1, linetype=3) +
    # stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + 
    stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), geom = 'text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), label.y = "bottom", label.x = 'right', size = 3) +
    geom_smooth(method='lm', se=FALSE, color="black", fullrange=TRUE, lwd=0.5) +
    scale_colour_manual(values = color_key) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    coord_cartesian(xlim = c(0, 400), ylim = c(0, 400)) +
    theme(legend.position = "none")
p6 <- ggplot(data = plotdata1, mapping = aes(x=BRD4, y=m2BRD4)) +
    ggrastr::rasterise(geom_point(mapping = aes(color=newclass), size = 0.5, stroke = 0, shape = 16), dpi = 400) +
    geom_abline(intercept = 0, slope = 1, linetype=3) +
    # stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + 
    stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), geom = 'text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), label.y = "bottom", label.x = 'right', size = 3) +
    geom_smooth(method='lm', se=FALSE, color="black", fullrange=TRUE, lwd=0.5) +
    scale_colour_manual(values = color_key) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    coord_cartesian(xlim = c(0, 900), ylim = c(0, 900)) +
    theme(legend.position = "none")

pdf(file='diff-dotplot.pdf', width=12, height=8)
plot((p1 | p2 | p3) / (p4 | p5 | p6)) 
dev.off()

setwd("/data/")
write.table(plotdata, file = "fig6a_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###################### Fig 6b ######################

plotdata1 <- plotdata %>% 
    mutate(ATAC = m2ATAC/ATAC, mH2A2 = m2mH2A2/mH2A2, H3K4me1 = m2H3K4me1/H3K4me1, H3K27ac = m2H3K27ac/H3K27ac, p300 = m2p300/p300, BRD4 = m2BRD4/BRD4) %>%
    select(ATAC, mH2A2, H3K4me1, H3K27ac, p300, BRD4) %>%
    pivot_longer(everything(), names_to = "mark", values_to = "signal") %>%
    mutate(mark = factor(mark, levels = c("ATAC", "mH2A2", "H3K4me1", "H3K27ac", "p300", "BRD4")))

p1 <- ggplot(plotdata1, mapping = aes(x=mark, y=signal)) +
    coord_cartesian(ylim = c(2^-1.5, 2^2.5)) +
    scale_y_continuous(trans="log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    ylab("log2 (FC)") +
    theme(axis.title.y = element_text(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(file='231L_ChIPSeq_m2OEvsGFP_boxplots.pdf', width=2.5, height=4.5)
p1
dev.off()

plotdata1 <- plotdata %>% 
    mutate(ATAC = m2ATAC/ATAC, mH2A2 = m2mH2A2/mH2A2, H3K4me1 = m2H3K4me1/H3K4me1, H3K27ac = m2H3K27ac/H3K27ac, p300 = m2p300/p300, BRD4 = m2BRD4/BRD4) %>%
    select(ATAC, mH2A2, H3K4me1, H3K27ac, p300, BRD4)
setwd("/data/")
write.table(plotdata1, file = "fig6b_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########################## Fig 6c ##################################

wd <- "/mforge/research/labs/experpath/maia/m237371/mBE/231L/cCREs"

setwd(wd)
signal_files <- c("BRD4_GFP.BRD4CR", "BRD4_m2OE.BRD4CR")
signals <- c("BRD4_GFP", "BRD4_m2OE")
rawdatas <- lapply(signal_files, read.table, skip=3)
rawdatas <- lapply(rawdatas, function(x){x %>% mutate(across(everything(), ~replace_na(.x, 0)))})
print(lapply(rawdatas, dim))

# Transformations
datas <- rawdatas
mydata <- bind_cols(lapply(datas, function(x){x %>% rowSums(.)}))
colnames(mydata) <- signals

peaks <- read.table("BRD4CR.sorted.bed") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname"))
peakstags <- bind_rows(read.table("BRD4SE_CR") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname")) %>%
    mutate(tag = "SE"),
    read.table("BRD4SC_CR") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname")) %>%
    mutate(tag = "Co"),
    read.table("BRD4LC_CR") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname")) %>%
    mutate(tag = "Co"),
    read.table("BRD4LE_CR") %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname")) %>%
    mutate(tag = "LE"))
peaks <- peaks %>% left_join(peakstags %>% select(peakname, tag))

mydata <- mydata %>% bind_cols(peaks)

   # Co    LE    SE
# 21836 27684  2151
library(ggpmisc)
plotdata1 <- mydata
p1 <- ggplot(data = plotdata1, mapping = aes(x=BRD4_GFP, y=BRD4_m2OE, color=tag)) +
    ggrastr::rasterise(geom_point(size = 0.5, stroke = 0, shape = 16), dpi = 400) +
    geom_abline(intercept = 0, slope = 1, linetype=3) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + 
    stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), geom = 'text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), label.y = "bottom", label.x = 'right', size = 3) +
    geom_smooth(method='lm', se=FALSE, fullrange=TRUE, lwd=0.5) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    coord_cartesian(xlim = c(0, 1500), ylim = c(0, 1500)) +
    theme_cowplot() + 
    theme(legend.position = "right")

pdf(file='231L_ChIPseq_BRD4_short_long_common.pdf', width=4.5, height=4)
p1
dev.off()

plotdata1 <- mydata %>% select(tag, BRD4_GFP, BRD4_m2OE)
setwd("/data/")
write.table(plotdata1, file = "fig6c_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########################## Fig 6d ##############################

plotdata1 <- plotdata %>% 
    select(ATAC, m2ATAC, mH2A2, m2mH2A2, H3K4me1, m2H3K4me1, H3K27ac, m2H3K27ac, p300, m2p300, BRD4, m2BRD4, newclass, peakname, chr, start, end, ccre) %>%
    mutate(mH2A2 = case_when(plotdata$mH2A2 == 0 ~ 1, TRUE ~ plotdata$mH2A2)) %>% 
    mutate(BRD4_FC = log2(m2BRD4/BRD4)) %>%
    mutate(BRD4_cl = case_when(BRD4_FC > log2(1.5) ~ "G", BRD4_FC < -log2(1.5) ~ "L", TRUE ~ "N"))

write.table(plotdata1 %>% filter(BRD4_cl == "G") %>% select(chr, start, end, peakname) %>% mutate(group = "Gain"), file = "BRD4_Gain", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(BRD4_cl == "N") %>% select(chr, start, end, peakname) %>% mutate(group = "Neut"), file = "BRD4_Neut", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(BRD4_cl == "L") %>% select(chr, start, end, peakname) %>% mutate(group = "Loss"), file = "BRD4_Loss", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

library(ReMapEnrich) 
BRD4_Loss <- bedToGranges("BRD4_Loss")
remapCatalog <- bedToGranges("/public_data/remap2022_nr_macs2_hg19_v1_0.bed")

BRD4_Loss_en <- enrichment(BRD4_Loss, remapCatalog)
BRD4_Loss_en <- BRD4_Loss_en %>% mutate(category = factor(category, levels = rev(BRD4_Loss_en$category)))

pdf(file='231L_ChIPSeq_BRD4LostPeaks_ReMapEnrich.pdf', width=6, height=4)
p2 <- ggplot(BRD4_Loss_en %>% slice_head(n = 10)) +
    geom_point(aes(x = mapped.peaks.ratio, y = category, size = nb.overlaps, color = q.significance)) +
    scale_color_gradient(low = "#6699ff", high = "#ff5050") +
    theme_cowplot() + theme(axis.title.y = element_blank()) +
    ggtitle("BRD4 Lost Peaks - ReMapEnrich")
p2
dev.off()

setwd("/data/")
write.table(BRD4_Loss_en %>% slice_head(n = 10) %>% select(mapped.peaks.ratio, category, nb.overlaps, q.significance), file = "fig6d_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########################## Supp Fig 6b ##############################
# BRD4=
# m2BRD4=
# ab=1000
# computeMatrix reference-point -S ${BRD4} ${m2BRD4} \
                              # -R BRD4_Loss BRD4_Neut BRD4_Gain \
                              # --referencePoint center \
                              # -b ${ab} -a ${ab} -p 16 \
                              # --missingDataAsZero \
                              # --outFileNameMatrix BRD4gainloss \
                              # --outFileName BRD4gainloss.tab.gz \
                              # --sortRegions keep
# plotHeatmap -m BRD4gainloss.tab.gz --colorList 'white,blue' --sortRegions descend --sortUsing mean --heatmapHeight 14 -out 231L_ChIPseq_BRD4gainloss_heatmap.pdf

########################## Supp Fig 6c ##############################

# ZMYND8=${sharedmBE}/ZMYND8/GSM2913937_X_1Z81.normalized.bigWig
# peaks=/mBE/231L/231L-ATAC-K4m1-outBL.bed
# ab=1000
# computeMatrix reference-point -S ${ZMYND8} \
                              # -R ${peaks} \
                              # --referencePoint center \
                              # -b ${ab} -a ${ab} -p 16 \
                              # --outFileNameMatrix ZMYND8 \
                              # --outFileName ZMYND8.tab.gz \
                              # --outFileSortedRegions ZMYND8_avail.bed

signals <- c("ZMYND8")
rawdatas <- lapply(signals, read.table, skip=3)
rawdatas <- lapply(rawdatas, function(x){x %>% mutate(across(everything(), ~replace_na(.x, 0)))})
print(lapply(rawdatas, dim))

# Transformations
datas <- rawdatas
mydata <- bind_cols(lapply(datas, function(x){x %>% rowSums(.)}))
colnames(mydata) <- signals

ZMYND8peaks <- read.table("/mBE/231L/ZMYND8_avail.bed") %>%
    select(V4) %>%
    set_colnames(c("peakname"))

mydata <- mydata %>% bind_cols(ZMYND8peaks)

plotdata <- plotdata %>% inner_join(mydata)
BRD4_Loss <- read.table(file = "BRD4_Loss", sep = "\t") %>% set_colnames(c("chr", "start", "end", "peakname", "BRD4status"))
BRD4_Gain <- read.table(file = "BRD4_Gain", sep = "\t") %>% set_colnames(c("chr", "start", "end", "peakname", "BRD4status"))
BRD4_Neut <- read.table(file = "BRD4_Neut", sep = "\t") %>% set_colnames(c("chr", "start", "end", "peakname", "BRD4status"))
BRD4_peak_classes <- bind_rows(BRD4_Loss, BRD4_Neut, BRD4_Gain)
plotdata <- plotdata %>% left_join(BRD4_peak_classes %>% select(peakname, BRD4status))

plotdata <- plotdata %>% mutate(across(c("ZMYND8"), function(x) {log(((x*10000)/sum(x))+0.01)}))
plotdata <- plotdata %>% mutate(across(c("ZMYND8"), function(x) {(x-mean(x))/sd(x)}))
pdf(file='231L_ChIPSeq_ZMYND8_boxplot.pdf', width=3, height=6)
p1 <- ggplot(plotdata, aes(x = BRD4status, y = ZMYND8)) +
    geom_boxplot(outlier.alpha = 0.1, outlier.shape = 1, outlier.size = 0.3) +
    # scale_y_continuous(trans=pseudolog10_trans) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1
dev.off()

setwd("/data/")
write.table(plotdata %>% select(BRD4status, ZMYND8), file = "suppfig6c_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########################## Fig 6g ##############################

all_signals <- read.table(file = "/mBE/all_signals.tab")
colnames(all_signals) <- c("CellLine", "newclass", "mark", "z")
plotdata <- all_signals %>% filter(!(mark %in% c("ATAC", "H4K12acVeh", "H4K12acE2", "H4K12ac"))) %>% 
    mutate(mark = factor(mark, levels = c("H3K27ac", "H3K4me1", "H2Az", "H3K4me3", "H3K27me3", "mH2A1", "mH2A2", "CTCF"))) %>%
    mutate(CellLine = factor(CellLine, levels = c("HMEC", "NHM", "MCF7", "231L", "HepG2")))

color_key <- c(Active = "#0f9448", APL = "#e78ac3", `ATAC-only` = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")
plotdata <- all_signals %>% filter(mark %in% c("H4K12ac", "H4K12acVeh")) %>% 
    mutate(mark = replace(mark, mark=="H4K12acVeh", "H4K12ac")) %>%
    mutate(newclass=replace(newclass, newclass=="H3K4me3", "APL")) %>%
    mutate(newclass=replace(newclass, newclass=="ATAConly", "ATAC-only")) %>%
    mutate(CellLine = factor(case_when(CellLine == "HMEC" ~ "HME1", CellLine == "NHM" ~ "Hmel", CellLine == "MCF7" ~ "MCF7", TRUE ~ CellLine), levels = c("HME1", "Hmel", "MCF7"))) %>%
    mutate(newclass = factor(newclass, levels = c("Active", "APL", "ATAC-only", "Inactive", "mBE")))
pdf(file='H4K12ac_boxplots.pdf', width=4, height=4)
ggplot(plotdata, aes(x = newclass, y = z, fill = newclass)) +
    facet_grid(cols = vars(CellLine)) +
    geom_boxplot(outlier.alpha = 0.1, outlier.shape = 1, outlier.size = 0.3) +
    coord_cartesian(ylim = c(-4, 4)) +
    scale_fill_manual(values = color_key) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
dev.off()

setwd("/data/")
write.table(plotdata %>% select(-mark), file = "fig6g_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


