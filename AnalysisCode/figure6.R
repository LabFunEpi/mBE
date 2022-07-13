
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

#######################################

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

pdf(file='diff-dotplot-temp.pdf', width=12, height=8)
plot((p1 | p2 | p3) / (p4 | p5 | p6)) 
dev.off()

theme_set(theme_cowplot() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
plotdata1 <- plotdata %>% 
    select(ATAC, m2ATAC) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal")
p1 <- ggplot(plotdata1, mapping = aes(x=condition, y=signal)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(trans="log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    ggtitle("ATAC")
plotdata1 <- plotdata %>% 
    select(mH2A2, m2mH2A2) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal")
p2 <- ggplot(plotdata1, mapping = aes(x=condition, y=signal)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(trans="log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    ggtitle("mH2A2")
plotdata1 <- plotdata %>% 
    select(H3K4me1, m2H3K4me1) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal")
p3 <- ggplot(plotdata1, mapping = aes(x=condition, y=signal)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(trans="log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    ggtitle("H3K4me1")
plotdata1 <- plotdata %>% 
    select(H3K27ac, m2H3K27ac) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal")
p4 <- ggplot(plotdata1, mapping = aes(x=condition, y=signal)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(trans="log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    ggtitle("H3K27ac")
plotdata1 <- plotdata %>% 
    select(p300, m2p300) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal")
p5 <- ggplot(plotdata1, mapping = aes(x=condition, y=signal)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(trans="log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    ggtitle("p300")
plotdata1 <- plotdata %>% 
    select(BRD4, m2BRD4) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal")
p6 <- ggplot(plotdata1, mapping = aes(x=condition, y=signal)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(trans="log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    ggtitle("BRD4")
pdf(file='diff-boxplot-temp.pdf', width=8, height=4)
(p1 | p2 | p3 | p4 | p5 | p6)
dev.off()

theme_set(theme_cowplot() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
statdata <- plotdata %>% 
    select(ATAC, m2ATAC) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal") %>%
    mutate(mark = "ATAC") %>%
    bind_rows(plotdata %>% 
    select(mH2A2, m2mH2A2) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal") %>%
    mutate(mark = "mH2A2")) %>%
    bind_rows(plotdata %>% 
    select(H3K4me1, m2H3K4me1) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal") %>%
    mutate(mark = "H3K4me1")) %>%
    bind_rows(plotdata %>% 
    select(H3K27ac, m2H3K27ac) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal") %>%
    mutate(mark = "H3K27ac")) %>%
    bind_rows(plotdata %>% 
    select(p300, m2p300) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal") %>%
    mutate(mark = "p300")) %>%
    bind_rows(plotdata %>% 
    select(BRD4, m2BRD4) %>% 
    set_colnames(c("CTRL", "m2OE")) %>%
    pivot_longer(everything(), names_to = "condition", values_to = "signal") %>%
    mutate(mark = "BRD4"))

stats <- compare_means(signal ~ condition, data = statdata, group.by = "mark", paired = TRUE)

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


#######################################

plotdata1 <- plotdata %>% 
    select(ATAC, m2ATAC, mH2A2, m2mH2A2, H3K4me1, m2H3K4me1, H3K27ac, m2H3K27ac, p300, m2p300, BRD4, m2BRD4, newclass, peakname, chr, start, end, ccre) %>%
    mutate(mH2A2 = case_when(plotdata$mH2A2 == 0 ~ 1, TRUE ~ plotdata$mH2A2)) %>% 
    mutate(ATAC_FC = log2(m2ATAC/ATAC), mH2A2_FC = log2(m2mH2A2/mH2A2), H3K4me1_FC = log2(m2H3K4me1/H3K4me1), H3K27ac_FC = log2(m2H3K27ac/H3K27ac), p300_FC = log2(m2p300/p300), BRD4_FC = log2(m2BRD4/BRD4)) %>%
    mutate(ATAC_cl = case_when(ATAC_FC > log2(1.5) ~ "G", ATAC_FC < -log2(1.5) ~ "L", TRUE ~ "N")) %>%
    mutate(mH2A2_cl = case_when(mH2A2_FC > log2(3) ~ "G", mH2A2_FC < -log2(1.5) ~ "L", TRUE ~ "N")) %>%
    mutate(H3K4me1_cl = case_when(H3K4me1_FC > log2(1.5) ~ "G", H3K4me1_FC < -log2(1.5) ~ "L", TRUE ~ "N")) %>%
    mutate(H3K27ac_cl = case_when(H3K27ac_FC > log2(1.5) ~ "G", H3K27ac_FC < -log2(1.5) ~ "L", TRUE ~ "N")) %>%
    mutate(p300_cl = case_when(p300_FC > log2(1.5) ~ "G", p300_FC < -log2(1.5) ~ "L", TRUE ~ "N")) %>%
    mutate(BRD4_cl = case_when(BRD4_FC > log2(1.5) ~ "G", BRD4_FC < -log2(1.5) ~ "L", TRUE ~ "N"))

write.table(plotdata1 %>% filter(ATAC_cl == "G") %>% select(chr, start, end, peakname) %>% mutate(group = "Gain"), file = "ATAC_Gain", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(ATAC_cl == "N") %>% select(chr, start, end, peakname) %>% mutate(group = "Neut"), file = "ATAC_Neut", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(ATAC_cl == "L") %>% select(chr, start, end, peakname) %>% mutate(group = "Loss"), file = "ATAC_Loss", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(mH2A2_cl == "G") %>% select(chr, start, end, peakname) %>% mutate(group = "Gain"), file = "mH2A2_Gain", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(mH2A2_cl == "N") %>% select(chr, start, end, peakname) %>% mutate(group = "Neut"), file = "mH2A2_Neut", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(mH2A2_cl == "L") %>% select(chr, start, end, peakname) %>% mutate(group = "Loss"), file = "mH2A2_Loss", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(H3K4me1_cl == "G") %>% select(chr, start, end, peakname) %>% mutate(group = "Gain"), file = "H3K4me1_Gain", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(H3K4me1_cl == "N") %>% select(chr, start, end, peakname) %>% mutate(group = "Neut"), file = "H3K4me1_Neut", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(H3K4me1_cl == "L") %>% select(chr, start, end, peakname) %>% mutate(group = "Loss"), file = "H3K4me1_Loss", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(H3K27ac_cl == "G") %>% select(chr, start, end, peakname) %>% mutate(group = "Gain"), file = "H3K27ac_Gain", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(H3K27ac_cl == "N") %>% select(chr, start, end, peakname) %>% mutate(group = "Neut"), file = "H3K27ac_Neut", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(H3K27ac_cl == "L") %>% select(chr, start, end, peakname) %>% mutate(group = "Loss"), file = "H3K27ac_Loss", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(p300_cl == "G") %>% select(chr, start, end, peakname) %>% mutate(group = "Gain"), file = "p300_Gain", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(p300_cl == "N") %>% select(chr, start, end, peakname) %>% mutate(group = "Neut"), file = "p300_Neut", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(p300_cl == "L") %>% select(chr, start, end, peakname) %>% mutate(group = "Loss"), file = "p300_Loss", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(BRD4_cl == "G") %>% select(chr, start, end, peakname) %>% mutate(group = "Gain"), file = "BRD4_Gain", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(BRD4_cl == "N") %>% select(chr, start, end, peakname) %>% mutate(group = "Neut"), file = "BRD4_Neut", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(plotdata1 %>% filter(BRD4_cl == "L") %>% select(chr, start, end, peakname) %>% mutate(group = "Loss"), file = "BRD4_Loss", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

library(ReMapEnrich) 
BRD4_Gain <- bedToGranges("BRD4_Gain")
BRD4_Loss <- bedToGranges("BRD4_Loss")
remapCatalog <- bedToGranges("/public_data/remap2022_nr_macs2_hg19_v1_0.bed")

BRD4_Gain_en <- enrichment(BRD4_Gain, remapCatalog)
BRD4_Gain_en <- BRD4_Gain_en %>% mutate(category = factor(category, levels = rev(BRD4_Gain_en$category)))
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
pdf(file='231L_ChIPSeq_ZMYND8_boxplot.pdf', width=6, height=6)
p1 <- ggplot(plotdata, aes(x = newclass, y = ZMYND8)) +
    geom_boxplot(outlier.alpha = 0.1, outlier.shape = 1, outlier.size = 0.3) +
    # scale_y_continuous(trans=pseudolog10_trans) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- ggplot(plotdata, aes(x = BRD4status, y = ZMYND8)) +
    geom_boxplot(outlier.alpha = 0.1, outlier.shape = 1, outlier.size = 0.3) +
    # scale_y_continuous(trans=pseudolog10_trans) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1 | p2
dev.off()






