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

#########################################  Fig 1c and Supp Fig 1b ###########################################################

setwd("/mBE/")
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

pdf(file='fig1c.pdf', width=8, height=3)
p1 <- ggplot(plotdata1 %>% filter(!(mark %in% c("H2Az", "CTCF"))), aes(x = mark, y = newclass)) +
    facet_grid(cols = vars(CellLine)) +
    geom_tile(aes(fill = median_z)) +
    coord_equal() +
    scale_fill_gradientn(colors = bluered(256), limits = c(-2.05, 2.05)) +
    guides(fill=guide_colorbar(ticks.colour = NA)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), strip.background = element_blank())
p1
dev.off()

write.table(plotdata1 %>% select(-notsig) %>% filter(!(mark %in% c("H2Az", "CTCF"))), file = "fig1c_data1.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

plotdata1 <- plotdata %>% filter(mark == "H3K4me1") %>% select(CellLine, newclass) %>% table() %>% data.frame() %>% filter(CellLine != "231L") %>%
    mutate(newclass=as.character(newclass)) %>%
    mutate(newclass=replace(newclass, newclass=="H3K4me3", "APL")) %>%
    mutate(newclass=replace(newclass, newclass=="ATAConly", "ATAC-only")) %>%
    mutate(newclass = factor(newclass, levels = rev(c("Active", "APL", "ATAC-only", "Inactive", "mBE")))) %>%
    mutate(CellLine = factor(CellLine, levels = rev(c("HMEC", "NHM", "MCF7", "HepG2"))))

color_key <- c(Active = "#0f9448", APL = "#e78ac3", `ATAC-only` = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")
pdf(file='fig1b_bar.pdf', width=9, height=2)
p1 <- ggplot(plotdata1, aes(x = Freq, y = CellLine, fill = newclass, label = Freq)) +
    geom_bar(position="fill", stat="identity") +
    geom_text(size = 4, position = position_fill(vjust = 0.5), color = "white") +
    scale_x_continuous(expand = c(0,0), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), position = "top") + 
    scale_fill_manual(values = color_key) +
    theme_cowplot() + theme(legend.position = "None", axis.title.y = element_blank())
p1
dev.off()

write.table(plotdata1, file = "fig1c_data2.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

plotdata1 <- plotdata %>% filter(CellLine != "231L") %>%
    mutate(newclass=replace(newclass, newclass=="H3K4me3", "APL")) %>%
    mutate(newclass=replace(newclass, newclass=="ATAConly", "ATAC-only")) %>%
    mutate(newclass = factor(newclass, levels = c("Active", "APL", "ATAC-only", "Inactive", "mBE"))) %>%
    mutate(CellLine = factor(CellLine, levels = c("HMEC", "NHM", "MCF7", "HepG2"))) %>%
    mutate(mark = factor(mark, levels = c("H3K4me1", "H3K27ac", "H2Az", "H3K4me3", "H3K27me3", "CTCF", "mH2A1", "mH2A2")))
pdf(file='fig1b-supp.pdf', width=12, height=10)
ggplot(plotdata1, aes(x = mark, y = z)) +
    facet_grid(cols = vars(newclass), rows = vars(CellLine)) +
    geom_boxplot(outlier.alpha = 0.1, outlier.shape = NA) +
    coord_cartesian(ylim = c(-4, 4)) +
    # theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

write.table(plotdata1, file = "suppfig1b_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

############ Fig 1e and Supp Fig 1c ###########################

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
# plotHeatmap -m H3K36me3.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H3K36me3.pdf --samplesLabel H3K36me3
# plotHeatmap -m H2BK12ac.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H2BK12ac.pdf --samplesLabel H2BK12ac
# plotHeatmap -m H2BK120ac.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out H2BK120ac.pdf --samplesLabel H2BK120ac

############ Fig 1d and Fig 3b (GAT) #################

# cd /softwares/homer/4.11/data/genomes/hg19/annotations/basic
# cat protein-coding.ann.txt promoters.ann.txt introns.ann.txt tts.ann.txt intergenic.ann.txt > /mBE/GAT/hg19_annotations_basic.txt
# cd /softwares/homer/4.11/data/genomes/mm9/annotations/basic
# cat protein-coding.ann.txt promoters.ann.txt introns.ann.txt tts.ann.txt intergenic.ann.txt > /mBE/GAT/mm9_annotations_basic.txt

# cd /mBE/GAT/
# cut -f2,3,4,6 hg19_annotations_basic.txt | awk 'BEGIN{OFS="\t";}$2>$3{tmp=$2;$2=$3;$3=tmp} 1' | awk '{ if (($2 >= 0) && ($3 >= 0)) { print } }' | sort -k1,1 -k2,2n | uniq > hg19_annotations_basic.bed
# cut -f2,3,4,6 mm9_annotations_basic.txt | awk 'BEGIN{OFS="\t";}$2>$3{tmp=$2;$2=$3;$3=tmp} 1' | awk '{ if (($2 >= 0) && ($3 >= 0)) { print } }' | sort -k1,1 -k2,2n | uniq > mm9_annotations_basic.bed

# cat hg19_annotations_basic.bed | awk '{count[$4]++} END {for (word in count) print word, count[word]}'
# cat hg19_annotations_basic.bed | awk '{count[$4]+=$3-$2} END {for (word in count) print word, count[word]}'

# cut -f2,3,4,6 /softwares/homer/4.11/data/genomes/hg19/annotations/basic/cpgIsland.ann.txt | awk 'BEGIN{OFS="\t";}$2>$3{tmp=$2;$2=$3;$3=tmp} 1' | awk '{ if (($2 >= 0) && ($3 >= 0)) { print } }' | sort -k1,1 -k2,2n | uniq > cpgIsland.ann.bed
# cut -f2,3,4,6 /softwares/homer/4.11/data/genomes/mm9/annotations/basic/cpgIsland.ann.txt | awk 'BEGIN{OFS="\t";}$2>$3{tmp=$2;$2=$3;$3=tmp} 1' | awk '{ if (($2 >= 0) && ($3 >= 0)) { print } }' | sort -k1,1 -k2,2n | uniq > mm9_cpgIsland.ann.bed

# grep "Active" HMEC_classes.bed > HMEC_Active.bed
# grep "H3K4me3" HMEC_classes.bed > HMEC_APL.bed
# grep "ATAConly" HMEC_classes.bed > HMEC_ATAConly.bed
# grep "Inactive" HMEC_classes.bed > HMEC_Inactive.bed
# grep "mBE" HMEC_classes.bed > HMEC_mBE.bed
# grep "Active" NHM_classes.bed > NHM_Active.bed
# grep "H3K4me3" NHM_classes.bed > NHM_APL.bed
# grep "ATAConly" NHM_classes.bed > NHM_ATAConly.bed
# grep "Inactive" NHM_classes.bed > NHM_Inactive.bed
# grep "mBE" NHM_classes.bed > NHM_mBE.bed
# grep "Active" MCF7_classes.bed > MCF7_Active.bed
# grep "H3K4me3" MCF7_classes.bed > MCF7_APL.bed
# grep "ATAConly" MCF7_classes.bed > MCF7_ATAConly.bed
# grep "Inactive" MCF7_classes.bed > MCF7_Inactive.bed
# grep "mBE" MCF7_classes.bed > MCF7_mBE.bed
# grep "Active" HepG2_classes.bed > HepG2_Active.bed
# grep "H3K4me3" HepG2_classes.bed > HepG2_APL.bed
# grep "ATAConly" HepG2_classes.bed > HepG2_ATAConly.bed
# grep "Inactive" HepG2_classes.bed > HepG2_Inactive.bed
# grep "mBE" HepG2_classes.bed > HepG2_mBE.bed
# grep "Active" DF_classes.bed > DF_Active.bed
# grep "H3K4me3" DF_classes.bed > DF_APL.bed
# grep "ATAConly" DF_classes.bed > DF_ATAConly.bed
# grep "Inactive" DF_classes.bed > DF_Inactive.bed
# grep "mBE" DF_classes.bed > DF_mBE.bed

# # hg19_map_contigs.bed and mm9_map_contigs.bed downloaded from UCSC Genome Browser
# cut -f1,2,3 hg19_map_contigs.bed | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"ws"}' > hg19_ws.bed
# cut -f1,2,3 mm9_map_contigs.bed | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,"ws"}' > mm9_ws.bed

# cat <(echo 'track name="Active"') HMEC_Active.bed <(echo 'track name="APL"') HMEC_APL.bed <(echo 'track name="ATAConly"') HMEC_ATAConly.bed <(echo 'track name="Inactive"') HMEC_Inactive.bed <(echo 'track name="mBE"') HMEC_mBE.bed > HMEC_tracks.bed
# cat <(echo 'track name="Active"') NHM_Active.bed <(echo 'track name="APL"') NHM_APL.bed <(echo 'track name="ATAConly"') NHM_ATAConly.bed <(echo 'track name="Inactive"') NHM_Inactive.bed <(echo 'track name="mBE"') NHM_mBE.bed > NHM_tracks.bed
# cat <(echo 'track name="Active"') MCF7_Active.bed <(echo 'track name="APL"') MCF7_APL.bed <(echo 'track name="ATAConly"') MCF7_ATAConly.bed <(echo 'track name="Inactive"') MCF7_Inactive.bed <(echo 'track name="mBE"') MCF7_mBE.bed > MCF7_tracks.bed
# cat <(echo 'track name="Active"') HepG2_Active.bed <(echo 'track name="APL"') HepG2_APL.bed <(echo 'track name="ATAConly"') HepG2_ATAConly.bed <(echo 'track name="Inactive"') HepG2_Inactive.bed <(echo 'track name="mBE"') HepG2_mBE.bed > HepG2_tracks.bed
# cat <(echo 'track name="Active"') DF_Active.bed <(echo 'track name="APL"') DF_APL.bed <(echo 'track name="ATAConly"') DF_ATAConly.bed <(echo 'track name="Inactive"') DF_Inactive.bed <(echo 'track name="mBE"') DF_mBE.bed > DF_tracks.bed

# gat-run.py --with-segment-tracks --segments=HMEC_tracks.bed --annotations=hg19_annotations_basic.bed --workspace=hg19_ws.bed --num-samples=1000 --isochore-file=cpgIsland.ann.bed --log=1 1> HMEC_tracks.GAT 2> err &
# gat-run.py --with-segment-tracks --segments=NHM_tracks.bed --annotations=hg19_annotations_basic.bed --workspace=hg19_ws.bed --num-samples=1000 --isochore-file=cpgIsland.ann.bed --log=2 1> NHM_tracks.GAT 2> err &
# gat-run.py --with-segment-tracks --segments=MCF7_tracks.bed --annotations=hg19_annotations_basic.bed --workspace=hg19_ws.bed --num-samples=1000 --isochore-file=cpgIsland.ann.bed --log=3 1> MCF7_tracks.GAT 2> err &
# gat-run.py --with-segment-tracks --segments=HepG2_tracks.bed --annotations=hg19_annotations_basic.bed --workspace=hg19_ws.bed --num-samples=1000 --isochore-file=cpgIsland.ann.bed --log=4 1> HepG2_tracks.GAT 2> err &

# gat-run.py --with-segment-tracks --segments=DF_tracks.bed --annotations=mm9_annotations_basic.bed --workspace=mm9_ws.bed --num-samples=1000 --isochore-file=mm9_cpgIsland.ann.bed --log=4 1> DF_tracks.GAT 2> err &

setwd("/mBE/GAT")
color_key <- c(Active = "#0f9448", APL = "#e78ac3", ATAConly = "#E0AC69", Inactive = "#f15a2b", mBE = "#2b598b")

files <- paste0(c("HMEC", "NHM", "MCF7", "HepG2"), "_tracks.GAT")
temp <- lapply(files, function(x){read.table(x, header = TRUE, na.strings = "na")})
names(temp) <- c("HMEC", "NHM", "MCF7", "HepG2")
out <- bind_rows(temp, .id = "cellline")
out <- out %>% mutate(class = factor(track, levels = names(color_key))) %>%
    mutate(annotation = factor(recode_factor(annotation, !!!c(`N` = "Intergenic", `E` = "Exon", `3UTR` = "3UTR", `I` = "Intron", `TTS` = "TTS", `5UTR` = "5UTR", `P` = "Promoter")), levels = c("Promoter", "5UTR", "Exon", "Intron", "3UTR", "TTS", "Intergenic"))) %>%
    mutate(negLogQval = -log10(qvalue), signif = negLogQval > -log10(0.05), cellline = factor(cellline, c("HMEC", "NHM", "MCF7", "HepG2")))

pdf(file = "gat.pdf", width = 5, height = 4)
ggplot(out, aes(x = annotation, y = l2fold, color = class, alpha = signif)) +
    geom_point(size = 2) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_manual(values=color_key) +
    facet_grid(cols = vars(cellline)) +
    coord_cartesian(ylim = c(-2.4, 1.8)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), legend.position = "None")
dev.off()

write.table(out, file = "fig1d_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

### DF
setwd("/mBE/GAT")
color_key <- c(Active = "#0f9448", APL = "#e78ac3", ATAConly = "#E0AC69", Inactive = "#f15a2b", mBE = "#2b598b")

out <- read.table("DF_tracks.GAT", header = TRUE, na.strings = "na")
out <- out %>% mutate(class = factor(track, levels = names(color_key))) %>%
    mutate(annotation = factor(recode_factor(annotation, !!!c(`N` = "Intergenic", `E` = "Exon", `3UTR` = "3UTR", `I` = "Intron", `TTS` = "TTS", `5UTR` = "5UTR", `P` = "Promoter")), levels = c("Promoter", "5UTR", "Exon", "Intron", "3UTR", "TTS", "Intergenic"))) %>%
    mutate(negLogQval = -log10(qvalue), signif = negLogQval > -log10(0.05))
pdf(file = "gat_DF.pdf", width = 2, height = 2.5)
ggplot(out, aes(x = annotation, y = l2fold, color = class, alpha = signif)) +
    geom_point(size = 2) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_manual(values=color_key) +
    # coord_cartesian(ylim = c(-2.4, 1.8)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), legend.position = "None")
dev.off()

write.table(out, file = "fig3b_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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

setwd("/mBE/HMEC/cCREs")
logRPKM <- readRDS("logRPKM.rds")
logRPKM <- logRPKM %>% filter(dataset %in% c("enc.polyA", "enc.total"))
level_key <- c(enc.polyA = "polyA", enc.total = "total")
enhExprData <- logRPKM %>% mutate(RPKM = 2^(logRPKM), dataset = recode(dataset, !!!level_key), cellline = "HMEC") %>% select(-logRPKM)

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

setwd("/mBE/NHM/cCREs")
logRPKM <- readRDS("logRPKM.rds")
logRPKM <- logRPKM %>% filter(dataset %in% c("S1"))
level_key <- c(S1 = "total")
enhExprData %<>% bind_rows(logRPKM %>% mutate(RPKM = 2^(logRPKM), dataset = recode(dataset, !!!level_key), cellline = "NHM") %>% select(-logRPKM))


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

setwd("/mBE/MCF7/cCREs")
logRPKM <- readRDS("logRPKM.rds")
logRPKM <- logRPKM %>% filter(dataset %in% c("enc.polyA"))
level_key <- c(enc.polyA = "total")
enhExprData %<>% bind_rows(logRPKM %>% mutate(RPKM = 2^(logRPKM), dataset = recode(dataset, !!!level_key), cellline = "MCF7") %>% select(-logRPKM))


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

setwd("/mBE/HepG2/cCREs")
logRPKM <- readRDS("logRPKM.rds")
logRPKM <- logRPKM %>% filter(dataset %in% c("enc.polyA", "enc.total"))
level_key <- c(enc.polyA = "polyA", enc.total = "total")
enhExprData %<>% bind_rows(logRPKM %>% mutate(RPKM = 2^(logRPKM), dataset = recode(dataset, !!!level_key), cellline = "HepG2") %>% select(-logRPKM))

level_key <- c(H3K4me3 = "APL")
enhExprData <- enhExprData %>% mutate(newclass = recode(newclass, !!!level_key))

setwd("/data/")
write.table(enhExprData %>% filter(dataset == "total") %>% select(-dataset), file = "fig2a_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(enhExprData %>% filter(dataset == "polyA") %>% select(-dataset), file = "suppfig2c_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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

setwd("/mBE/HMEC1/cCREs")
TCGAlogRPKM <- readRDS("TCGAlogRPKM.rds")
level_key <- c(H3K4me3 = "APL")
TCGAlogRPKM %<>% mutate(RPKM = 2^(logRPKM), newclass = recode(newclass, !!!level_key)) %>% select(-c(logRPKM, meanRPKM, peakname))

setwd("/data/")
write.table(TCGAlogRPKM, file = "fig2b_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


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

###################### Supp Fig 1d and 2d (Pairwise overlap plots with Intervene) #################

## Supp Fig 1d
# mkdir HMEC
# grep "Active" HMEC_classes.bed > HMEC/A_Active
# grep "H3K4me3" HMEC_classes.bed > HMEC/B_APL.bed
# grep "ATAConly" HMEC_classes.bed > HMEC/C_ATAConly.bed
# grep "Inactive" HMEC_classes.bed > HMEC/D_Inactive.bed
# grep "mBE" HMEC_classes.bed > HMEC/E_mBE.bed
# cp HMEC.lily.se.bed HMEC/F_SE.bed
# intervene pairwise -i HMEC/* --filenames --compute frac --htype color --output HMEC/ --figsize 10 10

# mkdir NHM
# grep "Active" NHM_classes.bed > NHM/A_Active.bed
# grep "H3K4me3" NHM_classes.bed > NHM/B_APL.bed
# grep "ATAConly" NHM_classes.bed > NHM/C_ATAConly.bed
# grep "Inactive" NHM_classes.bed > NHM/D_Inactive.bed
# grep "mBE" NHM_classes.bed > NHM/E_mBE.bed
# cp NHM.lily.se.bed NHM/F_SE.bed
# intervene pairwise -i NHM/* --filenames --compute frac --htype color --output NHM/ --figsize 10 10

# mkdir MCF7
# grep "Active" MCF7_classes.bed > MCF7/A_Active.bed
# grep "H3K4me3" MCF7_classes.bed > MCF7/B_APL.bed
# grep "ATAConly" MCF7_classes.bed > MCF7/C_ATAConly.bed
# grep "Inactive" MCF7_classes.bed > MCF7/D_Inactive.bed
# grep "mBE" MCF7_classes.bed > MCF7/E_mBE.bed
# cp MCF7.lily.se.bed MCF7/F_SE.bed
# intervene pairwise -i MCF7/* --filenames --compute frac --htype color --output MCF7/ --figsize 10 10

# mkdir HepG2
# grep "Active" HepG2_classes.bed > HepG2/A_Active.bed
# grep "H3K4me3" HepG2_classes.bed > HepG2/B_APL.bed
# grep "ATAConly" HepG2_classes.bed > HepG2/C_ATAConly.bed
# grep "Inactive" HepG2_classes.bed > HepG2/D_Inactive.bed
# grep "mBE" HepG2_classes.bed > HepG2/E_mBE.bed
# cp HepG2.lily.se.bed HepG2/F_SE.bed
# intervene pairwise -i HepG2/* --filenames --compute frac --htype color --output HepG2/ --figsize 10 10

## Supp Fig 2d
# intervene pairwise -i pairwise/* --filenames --compute frac --htype color --output pairwise/ --figsize 10 10

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

setwd("/mBE/HMEC1/cCREs")
proc_tbl <- readRDS("proc_tbl.rds")

setwd("/data/")
write.table(proc_tbl, file = "fig2c_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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

# intervene upset -i FINAL_CLASSES/Active/*.bed --names=NHM,HMEC,MCF7,HepG2 --output FINAL_CLASSES/Active/ --figsize 5 4
# intervene upset -i FINAL_CLASSES/Inactive/*.bed --names=NHM,HMEC,MCF7,HepG2 --output FINAL_CLASSES/Inactive/ --figsize 5 4
# intervene upset -i FINAL_CLASSES/ATAConly/*.bed --names=NHM,HMEC,MCF7,HepG2 --output FINAL_CLASSES/ATAConly/ --figsize 5 4
# intervene upset -i FINAL_CLASSES/CRE/*.bed --names=NHM,HMEC,MCF7,HepG2 --output FINAL_CLASSES/CRE/ --figsize 5 4

# intervene upset -i FINAL_CLASSES/mBE/*.bed --names=NHM,HMEC,MCF7,HepG2 --output FINAL_CLASSES/mBE/ --figsize 5 3 --order degree
# intervene upset -i FINAL_CLASSES/APL/*.bed --names=NHM,HMEC,MCF7,HepG2 --output FINAL_CLASSES/APL/ --figsize 5 3

setwd("/mBE/")
library("UpSetR")

pdf("upset_mBE.pdf", width=5, height=3, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('HepG2'=11439,'MCF7'=7172,'MCF7&HepG2'=605,'HMEC'=15721,'HMEC&HepG2'=569,'HMEC&MCF7'=1542,'HMEC&MCF7&HepG2'=210,'NHM'=5243,'NHM&HepG2'=154,'NHM&MCF7'=196,'NHM&MCF7&HepG2'=20,'NHM&HMEC'=648,'NHM&HMEC&HepG2'=33,'NHM&HMEC&MCF7'=96,'NHM&HMEC&MCF7&HepG2'=15)
upset(fromExpression(expressionInput), nsets=4, sets = rev(c("HMEC", "MCF7", "HepG2", "NHM")), nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size", keep.order = TRUE, order.by = "freq")
invisible(dev.off())

pdf("upset_APL.pdf", width=5, height=3, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('HepG2'=6185,'MCF7'=4937,'MCF7&HepG2'=1803,'HMEC'=4328,'HMEC&HepG2'=1335,'HMEC&MCF7'=2070,'HMEC&MCF7&HepG2'=1676,'NHM'=5679,'NHM&HepG2'=123,'NHM&MCF7'=38,'NHM&MCF7&HepG2'=27,'NHM&HMEC'=117,'NHM&HMEC&HepG2'=33,'NHM&HMEC&MCF7'=28,'NHM&HMEC&MCF7&HepG2'=18)
upset(fromExpression(expressionInput), nsets=4, sets = rev(c("HMEC", "MCF7", "HepG2", "NHM")), nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size", keep.order = TRUE, order.by = "freq")
invisible(dev.off())

pdf("upset_Active.pdf", width=4, height=3, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('HepG2'=17760,'MCF7'=6013,'MCF7&HepG2'=877,'HMEC'=14078,'HMEC&HepG2'=684,'HMEC&MCF7'=1667,'HMEC&MCF7&HepG2'=376,'NHM'=6598,'NHM&HepG2'=236,'NHM&MCF7'=133,'NHM&MCF7&HepG2'=24,'NHM&HMEC'=816,'NHM&HMEC&HepG2'=70,'NHM&HMEC&MCF7'=155,'NHM&HMEC&MCF7&HepG2'=31)
upset(fromExpression(expressionInput), nsets=4, sets = rev(c("HMEC", "MCF7", "HepG2", "NHM")), nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size", keep.order = TRUE, order.by = "freq")
invisible(dev.off())

pdf("upset_ATAConly.pdf", width=4, height=3, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('HepG2'=12232,'MCF7'=5634,'MCF7&HepG2'=233,'HMEC'=11668,'HMEC&HepG2'=126,'HMEC&MCF7'=236,'HMEC&MCF7&HepG2'=9,'NHM'=5734,'NHM&HepG2'=44,'NHM&MCF7'=9,'NHM&MCF7&HepG2'=1,'NHM&HMEC'=497,'NHM&HMEC&HepG2'=6,'NHM&HMEC&MCF7'=1,'NHM&HMEC&MCF7&HepG2'=0)
upset(fromExpression(expressionInput), nsets=4, sets = rev(c("HMEC", "MCF7", "HepG2", "NHM")), nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size", keep.order = TRUE, order.by = "freq")
invisible(dev.off())

pdf("upset_Inactive.pdf", width=4, height=3, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('HepG2'=9509,'MCF7'=5314,'MCF7&HepG2'=304,'HMEC'=4785,'HMEC&HepG2'=564,'HMEC&MCF7'=1741,'HMEC&MCF7&HepG2'=137,'NHM'=8960,'NHM&HepG2'=22,'NHM&MCF7'=61,'NHM&MCF7&HepG2'=1,'NHM&HMEC'=43,'NHM&HMEC&HepG2'=0,'NHM&HMEC&MCF7'=9,'NHM&HMEC&MCF7&HepG2'=0)
upset(fromExpression(expressionInput), nsets=4, sets = rev(c("HMEC", "MCF7", "HepG2", "NHM")), nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size", keep.order = TRUE, order.by = "freq")
invisible(dev.off())

pdf("upset_CRE.pdf", width=4, height=3, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('HepG2'=33500,'MCF7'=15807,'MCF7&HepG2'=6602,'HMEC'=36780,'HMEC&HepG2'=4235,'HMEC&MCF7'=9955,'HMEC&MCF7&HepG2'=6729,'NHM'=24164,'NHM&HepG2'=1732,'NHM&MCF7'=986,'NHM&MCF7&HepG2'=489,'NHM&HMEC'=5390,'NHM&HMEC&HepG2'=622,'NHM&HMEC&MCF7'=1422,'NHM&HMEC&MCF7&HepG2'=712)
upset(fromExpression(expressionInput), nsets=4, sets = rev(c("HMEC", "MCF7", "HepG2", "NHM")), nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size", keep.order = TRUE, order.by = "freq")
invisible(dev.off())

venndata1 <- data.frame(
mBE = c('HepG2'=11439,'MCF7'=7172,'MCF7&HepG2'=605,'HMEC'=15721,'HMEC&HepG2'=569,'HMEC&MCF7'=1542,'HMEC&MCF7&HepG2'=210,'NHM'=5243,'NHM&HepG2'=154,'NHM&MCF7'=196,'NHM&MCF7&HepG2'=20,'NHM&HMEC'=648,'NHM&HMEC&HepG2'=33,'NHM&HMEC&MCF7'=96,'NHM&HMEC&MCF7&HepG2'=15), 
APL = c('HepG2'=6185,'MCF7'=4937,'MCF7&HepG2'=1803,'HMEC'=4328,'HMEC&HepG2'=1335,'HMEC&MCF7'=2070,'HMEC&MCF7&HepG2'=1676,'NHM'=5679,'NHM&HepG2'=123,'NHM&MCF7'=38,'NHM&MCF7&HepG2'=27,'NHM&HMEC'=117,'NHM&HMEC&HepG2'=33,'NHM&HMEC&MCF7'=28,'NHM&HMEC&MCF7&HepG2'=18)
)
venndata2 <- data.frame(
Active = c('HepG2'=17760,'MCF7'=6013,'MCF7&HepG2'=877,'HMEC'=14078,'HMEC&HepG2'=684,'HMEC&MCF7'=1667,'HMEC&MCF7&HepG2'=376,'NHM'=6598,'NHM&HepG2'=236,'NHM&MCF7'=133,'NHM&MCF7&HepG2'=24,'NHM&HMEC'=816,'NHM&HMEC&HepG2'=70,'NHM&HMEC&MCF7'=155,'NHM&HMEC&MCF7&HepG2'=31),
ATAConly = c('HepG2'=12232,'MCF7'=5634,'MCF7&HepG2'=233,'HMEC'=11668,'HMEC&HepG2'=126,'HMEC&MCF7'=236,'HMEC&MCF7&HepG2'=9,'NHM'=5734,'NHM&HepG2'=44,'NHM&MCF7'=9,'NHM&MCF7&HepG2'=1,'NHM&HMEC'=497,'NHM&HMEC&HepG2'=6,'NHM&HMEC&MCF7'=1,'NHM&HMEC&MCF7&HepG2'=0),
Inactive = c('HepG2'=9509,'MCF7'=5314,'MCF7&HepG2'=304,'HMEC'=4785,'HMEC&HepG2'=564,'HMEC&MCF7'=1741,'HMEC&MCF7&HepG2'=137,'NHM'=8960,'NHM&HepG2'=22,'NHM&MCF7'=61,'NHM&MCF7&HepG2'=1,'NHM&HMEC'=43,'NHM&HMEC&HepG2'=0,'NHM&HMEC&MCF7'=9,'NHM&HMEC&MCF7&HepG2'=0),
CRE = c('HepG2'=33500,'MCF7'=15807,'MCF7&HepG2'=6602,'HMEC'=36780,'HMEC&HepG2'=4235,'HMEC&MCF7'=9955,'HMEC&MCF7&HepG2'=6729,'NHM'=24164,'NHM&HepG2'=1732,'NHM&MCF7'=986,'NHM&MCF7&HepG2'=489,'NHM&HMEC'=5390,'NHM&HMEC&HepG2'=622,'NHM&HMEC&MCF7'=1422,'NHM&HMEC&MCF7&HepG2'=712)
)

setwd("/data/")
write.table(venndata1, file = "fig2d_data.tab", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(venndata2, file = "suppfig2e_data.tab", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

##################### Fig 2e ######################

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

cistromedata <- bind_rows(list(mBE = mBE, APL = APL), .id = "newclass") %>% select(-rowname)
setwd("/data/")
write.table(cistromedata, file = "fig2e_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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

setwd("/data/")
cstat_files <- c("/mBE/HMEC/cCREs/cstats.txt", 
    "/mBE/NHM/cCREs/cstats.txt", 
    "/mBE/MCF7/cCREs/cstats.txt", 
    "/mBE/HepG2/cCREs/cstats.txt")
cell_lines <- c("HMEC", "NHM", "MCF7", "HepG2")
cstats <- lapply(cstat_files, read.table)
names(cstats) <- cell_lines
cstats <- bind_rows(cstats, .id = "cell.line") %>% 
    set_colnames(c("cell.line", "k", "avg_silhouette_score", "CH_score")) %>%
    mutate(k = factor(k, levels = as.character(2:11))) %>%
    mutate(cell.line = factor(cell.line, levels = cell_lines)) %>%
    select(-CH_score)
write.table(cstats, file = "suppfig1a_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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

counts %<>% filter(pos == 50) %>% ungroup() %>% select(c(chstate, newclass, prop))
setwd("/data/")
write.table(counts, file = "suppfig2b_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

















