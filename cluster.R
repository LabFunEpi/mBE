suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(ggpubr))

args <- commandArgs(trailingOnly = TRUE)
candidate_enhancers.bed <- args[1]
candidates_wENCODE.bed <- args[2]
outdir <- args[3]

setwd(outdir)

signals <- c("H3K27ac", "H3K4me1", "H3K4me3", "H2Az", "H3K27me3", "mH2A1", "mH2A2", "CTCF")
rawdatas <- lapply(signals, read.table, skip=3)
rawdatas <- lapply(rawdatas, function(x){x %>% mutate(across(everything(), ~replace_na(.x, 0)))})

# Transformations
datas <- rawdatas
mydata <- bind_cols(lapply(datas, function(x){x %>% rowSums(.)}))
colnames(mydata) <- signals

# Normalization
mydata <- mydata %>% mutate(across(everything(), function(x) {log(((x*10000)/sum(x))+0.01)}))
mydata <- mydata %>% mutate(across(everything(), function(x) {(x-mean(x))/sd(x)}))

candidates_wENCODE <- read.table(candidates_wENCODE.bed) %>%
    select(V4,V8,V14) %>%
    set_colnames(c("peakname", "encode_id", "encode_class")) %>%
    group_by(peakname) %>% dplyr::slice(1)
candidate_enhancers <- read.table(candidate_enhancers.bed) %>%
    select(V1:V4) %>%
    set_colnames(c("chr", "start", "end", "peakname")) %>%
    left_join(candidates_wENCODE)

mydata <- mydata %>% bind_cols(candidate_enhancers) %>% filter(encode_class != ".")

km5 <- mydata %>% select(all_of(signals)) %>% kmeans(centers = 5, nstart = 20)
km5$cluster <- factor(km5$cluster)
mydata <- mydata %>% mutate(cluster=km5$cluster)

write.table(mydata %>% select(chr, start, end, peakname, encode_id, encode_class, cluster), file = "cluster_results.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

p1 <- ggplot(mydata) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=cluster), size = 0.5, stroke = 0, shape = 16) +
    annotate("text", x = km5$centers[, 7], y = km5$centers[, 1], fontface=2, label = table(km5$cluster)) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme_cowplot()

p2 <- ggplot(mydata) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=H3K4me1), size = 0.5, stroke = 0, shape = 16) +
    facet_wrap(vars(cluster)) +
    scale_colour_gradient2(midpoint = 0, low = "red", mid = "grey", high = "blue", limits = NULL) +
    theme_cowplot()
    
p3 <- ggplot(mydata) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=H3K4me3), size = 0.5, stroke = 0, shape = 16) +
    facet_wrap(vars(cluster)) +
    scale_colour_gradient2(midpoint = 0, low = "red", mid = "grey", high = "blue", limits = NULL) +
    theme_cowplot()

p4 <- ggplot(mydata) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=H2Az), size = 0.5, stroke = 0, shape = 16) +
    facet_wrap(vars(cluster)) +
    scale_colour_gradient2(midpoint = 0, low = "red", mid = "grey", high = "blue", limits = NULL) +
    theme_cowplot()
    
p5 <- ggplot(mydata) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=H3K27me3), size = 0.5, stroke = 0, shape = 16) +
    facet_wrap(vars(cluster)) +
    scale_colour_gradient2(midpoint = 0, low = "red", mid = "grey", high = "blue", limits = NULL) +
    theme_cowplot()
    
p6 <- ggplot(mydata) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=mH2A1), size = 0.5, stroke = 0, shape = 16) +
    facet_wrap(vars(cluster)) +
    scale_colour_gradient2(midpoint = 0, low = "red", mid = "grey", high = "blue", limits = NULL) +
    theme_cowplot()
    
p7 <- ggplot(mydata) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=CTCF), size = 0.5, stroke = 0, shape = 16) +
    facet_wrap(vars(cluster)) +
    scale_colour_gradient2(low = "red", mid = "grey", high = "blue", limits = NULL) +
    theme_cowplot()

forlegend <- ggplot(mydata) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=cluster), size = 0.8, stroke = 0, shape = 16) +
    theme_cowplot() +
    guides(colour = guide_legend(override.aes = list(size=5)))
mylegend <- get_legend(forlegend)
p8 <- ggplot(mydata %>% filter(encode_class == "CTCF-only")) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=cluster), size = 0.8, stroke = 0, shape = 16) +
    theme_cowplot() +
    theme(legend.position = "none") +
    ggtitle("CTCF-only")

p9 <- ggplot(mydata %>% filter(encode_class == "DNase-H3K4me3")) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=cluster), size = 0.8, stroke = 0, shape = 16) +
    theme_cowplot() +
    theme(legend.position = "none") +
    ggtitle("DNase-H3K4me3")

p10 <- ggplot(mydata %>% filter(encode_class == "PLS")) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=cluster), size = 0.5, stroke = 0, shape = 16) +
    theme_cowplot() +
    theme(legend.position = "none") +
    ggtitle("PLS")

p11 <- ggplot(mydata %>% filter(encode_class == "dELS")) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=cluster), size = 0.2, stroke = 0, shape = 16) +
    theme_cowplot() +
    theme(legend.position = "none") +
    ggtitle("dELS")

p12 <- ggplot(mydata %>% filter(encode_class == "pELS")) +
    geom_point(aes(x=mH2A2, y=H3K27ac, color=cluster), size = 0.4, stroke = 0, shape = 16) +
    theme_cowplot() +
    theme(legend.position = "none") +
    ggtitle("pELS")

pdf(file='plots.pdf', width=12, height=10)
plot(p1)
plot(p2)
plot(p3)
plot(p4)
plot(p5)
plot(p6)
plot(p7)
plot((p8 | p9 | p10) / (p11 | p12 | as_ggplot(mylegend)))
dev.off()
