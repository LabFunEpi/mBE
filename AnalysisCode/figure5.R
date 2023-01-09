########### Fig 5a ##############

setwd("/mBE/MCF7/cCREs")
library(ReMapEnrich) 
Active <- bedToGranges("Active_3col.bed")
Inactive <- bedToGranges("Inactive_3col.bed")
H3K4me3 <- bedToGranges("H3K4me3_3col.bed")
ATAConly <- bedToGranges("ATAConly_3col.bed")
mBE <- bedToGranges("mBE_3col.bed")
remapCatalog <- bedToGranges("/public_data/remap2022/remap2022_MCF-7_all_macs2_hg19_v1_0.bed")

en_Active <- enrichment(Active, remapCatalog, chromSizes = loadChromSizes("hg19"), byChrom = TRUE)
en_Inactive <- enrichment(Inactive, remapCatalog, chromSizes = loadChromSizes("hg19"), byChrom = TRUE)
en_H3K4me3 <- enrichment(H3K4me3, remapCatalog, chromSizes = loadChromSizes("hg19"), byChrom = TRUE)
en_ATAConly <- enrichment(ATAConly, remapCatalog, chromSizes = loadChromSizes("hg19"), byChrom = TRUE)
en_mBE <- enrichment(mBE, remapCatalog, chromSizes = loadChromSizes("hg19"), byChrom = TRUE)

en_list <- list(en_Active, en_H3K4me3, en_ATAConly, en_Inactive, en_mBE)
names(en_list) <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")
en <- bind_rows(en_list, .id = "class")
en %<>% filter(q.value < 0.05) %>% 
    mutate(l2eff = log2(effect.size)) %>%
    select(class, category, effect.size) %>% pivot_wider(names_from = "class", values_from = "effect.size") %>% 
    drop_na() %>%
    pivot_longer(!category, names_to = "class", values_to = "effect.size")
plotorder <- en %>% select(class, category, effect.size) %>% pivot_wider(names_from = "class", values_from = "effect.size") %>%
    arrange(-mBE) %>% select(category) %>% unlist() %>% unname()
en %<>% mutate(category = factor(category, levels = plotorder), class = factor(class, levels = c("H3K4me3", "Active", "ATAConly", "Inactive", "mBE")))

write.table(en, file = "ReMapEnrich_results.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
en <- read.table("ReMapEnrich_results.tsv", sep = "\t", header = TRUE)
plotorder <- en %>% select(class, category, effect.size) %>% pivot_wider(names_from = "class", values_from = "effect.size") %>%
    arrange(-mBE) %>% select(category) %>% unlist() %>% unname()
en %<>% mutate(category = factor(category, levels = plotorder), class = factor(class, levels = c("H3K4me3", "Active", "ATAConly", "Inactive", "mBE")))

pdf(file='MCF7_scATAC_ReMapEnrich.pdf', width=28, height=6)
p1 <- ggplot(en, aes(x = category, y = class)) +
    geom_tile(aes(fill = effect.size)) +
    scale_fill_viridis_c(option = "inferno") +
    scale_x_discrete(position = "top") +
    guides(fill=guide_colorbar(ticks.colour = NA)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
p1
dev.off()

######## Cell Ranger Runs #########

# Running Cell Ranger with MCF7 specific CRE as the input peak set ...
# grep -v -E 'chr7\_KI270803v1\_alt' classified_liftedto_hg38.bed > classified_liftedto_hg38_clean.bed
# sort -k1,1 -k2,2n classified_liftedto_hg38_clean.bed > classified_liftedto_hg38_clean_sorted.bed
# cut -f1,2,3 classified_liftedto_hg38_clean_sorted.bed > classified_hg38_3col_clean_sorted.bed
# readsdir=Fastq
# cratac=scATAC_count_stencil.sh
# cratacaggr=scATAC_aggr_stencil.sh
# crout=scATAC_hg38_stencil
# peaks=classified_hg38_3col_clean_sorted.bed

# qsub -N Rep1_WT -v id=Rep1_WT,sample=MCF7_0_5_E2_scATAC,readsdir=${readsdir},outdir=${crout},peaks=${peaks} ${cratac}
# qsub -N Rep1_KO -v id=Rep1_KO,sample=MCF7_2_1_E2_scATAC,readsdir=${readsdir},outdir=${crout},peaks=${peaks} ${cratac}
# qsub -N Rep2_WT -v id=Rep2_WT,sample=0_5_E2,readsdir=${readsdir},outdir=${crout},peaks=${peaks} ${cratac}
# qsub -N Rep2_KO -v id=Rep2_KO,sample=2_1_E2,readsdir=${readsdir},outdir=${crout},peaks=${peaks} ${cratac}
# qsub -N aggr -hold_jid Rep1_WT,Rep1_KO,Rep2_WT,Rep2_KO -v id=aggr,libraries=${crout}/aggr.csv,outdir=${crout},peaks=${peaks} ${cratacaggr}

# Running Cell Ranger with Cell Ranger peak calling ... 
# readsdir=Fastq
# cratac=scATAC_count.sh
# cratacaggr=scATAC_aggr.sh
# crout=scATAC_hg38

# qsub -N Rep1_WT -v id=Rep1_WT,sample=MCF7_0_5_E2_scATAC,readsdir=${readsdir},outdir=${crout} ${cratac}
# qsub -N Rep1_KO -v id=Rep1_KO,sample=MCF7_2_1_E2_scATAC,readsdir=${readsdir},outdir=${crout} ${cratac}
# qsub -N Rep2_WT -v id=Rep2_WT,sample=0_5_E2,readsdir=${readsdir},outdir=${crout} ${cratac}
# qsub -N Rep2_KO -v id=Rep2_KO,sample=2_1_E2,readsdir=${readsdir},outdir=${crout} ${cratac}
# qsub -N aggr -hold_jid Rep1_WT,Rep1_KO,Rep2_WT,Rep2_KO -v id=aggr,libraries=${crout}/aggr.csv,outdir=${crout} ${cratacaggr}
# qsub -N aggrWT -v id=aggrWT,libraries=${crout}/aggrWT.csv,outdir=${crout} ${cratacaggr}
# qsub -N aggrKO -v id=aggrKO,libraries=${crout}/aggrKO.csv,outdir=${crout} ${cratacaggr}

######## Get Cell Ranger Summary for QC supplementary table ########
setwd("/mBE/MCF7/scATAC_hg38_stencil/")
samplenames <- c("Rep1_WT", "Rep1_KO", "Rep2_WT", "Rep2_KO")

cr_summary <- lapply(paste0(samplenames, "/outs/summary.csv"), function(x){read.csv(x) %>% t() %>% t() %>% data.frame()})
names(cr_summary) <- samplenames
cr_summary <- bind_rows(cr_summary) %>% t()

write.table(cr_summary, file = "summary.table", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(read.csv("aggr/outs/summary.csv") %>% t() %>% data.frame(), file = "aggr.summary.table", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

###################### Prepare Seurat Object - followed by QC filters and sub-sampling of cells ##################
h5 <- "/mBE/MCF7/scATAC_hg38_stencil/aggr/outs/filtered_peak_bc_matrix.h5"
frag.file <- "/mBE/MCF7/scATAC_hg38_stencil/aggr/outs/fragments.tsv.gz"
metadata_csv_1 <- "/mBE/MCF7/scATAC_hg38_stencil/Rep1_WT/outs/singlecell.csv"
metadata_csv_2 <- "/mBE/MCF7/scATAC_hg38_stencil/Rep1_KO/outs/singlecell.csv"
metadata_csv_3 <- "/mBE/MCF7/scATAC_hg38_stencil/Rep2_WT/outs/singlecell.csv"
metadata_csv_4 <- "/mBE/MCF7/scATAC_hg38_stencil/Rep2_KO/outs/singlecell.csv"

atac_counts <- Read10X_h5(h5)
metadata_1 <- read.csv(file = metadata_csv_1, header = TRUE, row.names = 1) %>% filter(is__cell_barcode == 1) %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", NA), sep = "-") %>% mutate(barcode = paste0(barcode, "-1"))
metadata_2 <- read.csv(file = metadata_csv_2, header = TRUE, row.names = 1) %>% filter(is__cell_barcode == 1) %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", NA), sep = "-") %>% mutate(barcode = paste0(barcode, "-2"))
metadata_3 <- read.csv(file = metadata_csv_3, header = TRUE, row.names = 1) %>% filter(is__cell_barcode == 1) %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", NA), sep = "-") %>% mutate(barcode = paste0(barcode, "-3"))
metadata_4 <- read.csv(file = metadata_csv_4, header = TRUE, row.names = 1) %>% filter(is__cell_barcode == 1) %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", NA), sep = "-") %>% mutate(barcode = paste0(barcode, "-4"))
metadata <- data.frame(barcode = colnames(atac_counts)) %>% left_join(bind_rows(metadata_1, metadata_2, metadata_3, metadata_4)) %>% column_to_rownames("barcode")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(":", "-"), genome = 'hg38', fragments = frag.file, min.cells = 10, annotation = annotations)
sobj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", meta.data = metadata)

sobj <- NucleosomeSignal(object = sobj)
sobj <- TSSEnrichment(object = sobj, fast = FALSE)

DefaultAssay(sobj) <- "ATAC"
sobj <- sobj %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:50)

mapping <- data.frame(sampleno = c("1", "2", "3", "4"), sample = c("Rep1_WT", "Rep1_KO", "Rep2_WT", "Rep2_KO"))
all.metadata <- sobj@meta.data %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", "sampleno"), sep = "-") %>% left_join(mapping)
sobj$sample <- all.metadata$sample
saveRDS(sobj, "aggr_seuratobj.rds")

sobj_filtered <- subset(sobj, subset = nFeature_ATAC > 200 & TSS.enrichment > 1)
ncells <- min(table(sobj_filtered$sample))
sub_cells <- c(
    sample(sobj_filtered@meta.data %>% filter(sample == "Rep1_WT") %>% rownames(), size = ncells, replace=F),
    sample(sobj_filtered@meta.data %>% filter(sample == "Rep1_KO") %>% rownames(), size = ncells, replace=F),
    sample(sobj_filtered@meta.data %>% filter(sample == "Rep2_WT") %>% rownames(), size = ncells, replace=F),
    sample(sobj_filtered@meta.data %>% filter(sample == "Rep2_KO") %>% rownames(), size = ncells, replace=F)
)
sobj_sub <- sobj_filtered[,sub_cells]

DefaultAssay(sobj_sub) <- "ATAC"
sobj_sub <- sobj_sub %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:20) %>% FindNeighbors(reduction = 'lsi', dims = 2:20)
sobj_sub <- sobj_sub %>% FindClusters(algorithm = 3, resolution = 0.15)

saveRDS(sobj_sub, "aggr_filtered_sub.rds")

##################################

setwd("/mBE/MCF7/scATAC_hg38_stencil/")
datadir <- "/mBE/MCF7/scATAC_hg38_stencil/"
samplenames <- c("Rep1_WT", "Rep1_KO", "Rep2_WT", "Rep2_KO")
library(ggallin)
library(scico)
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg38"
# saveRDS(annotations, "/references/annotations.rds")
annotations <- readRDS("/references/annotations.rds")
theme_set(theme_cowplot())
genome <- BSgenome.Hsapiens.UCSC.hg38
keepBSgenomeSequences <- function(genome, seqnames)
{
    stopifnot(all(seqnames %in% seqnames(genome)))
    genome@user_seqnames <- setNames(seqnames, seqnames)
    genome@seqinfo <- genome@seqinfo[seqnames]
    genome
}
sequences_to_keep <- paste0("chr", c(1:22, "X", "Y"))
genome <- keepBSgenomeSequences(genome, sequences_to_keep)
library(BiocParallel)
register(SerialParam())
theme_set(theme_cowplot())

sobj_sub <- readRDS("aggr_filtered_sub.rds")

############ Supp Fig 5d (UMAPS) ###########
DefaultAssay(sobj_sub) <- "ATAC"
Idents(sobj_sub) <- "seurat_clusters"
sobj_sub@meta.data <- sobj_sub@meta.data %>% separate(col = "sample", into = c("replicate", "genotype"), sep = "_", remove = FALSE)
pdf(file='MCF7_scATAC_umap.pdf', width=5, height=2.5)
ggplot(data.frame(sobj_sub[["umap"]][[]], genotype = factor(sobj_sub$genotype, levels = c("WT", "KO")), clust = Idents(sobj_sub)), aes(x=UMAP_1, y=UMAP_2, color = clust)) + 
    ggrastr::rasterise(geom_point(size=0.3, stroke=0.1, shape=16, show.legend = TRUE), dpi = 400) +
    facet_grid(cols = vars(genotype)) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.title = element_blank(), strip.background = element_rect(fill = "white"))
ggplot(data.frame(sobj_sub[["umap"]][[]], genotype = factor(sobj_sub$genotype, levels = c("WT", "KO")), sample = factor(sobj_sub$sample, levels = c("Rep1_WT", "Rep2_WT", "Rep1_KO", "Rep2_KO"))), aes(x=UMAP_1, y=UMAP_2, color = sample)) + 
    ggrastr::rasterise(geom_point(size=0.3, stroke=0.1, shape=16, show.legend = TRUE), dpi = 400) +
    facet_grid(cols = vars(genotype)) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.title = element_blank(), strip.background = element_rect(fill = "white"))
dev.off()

setwd("/data/")
tab <- data.frame(sobj_sub[["umap"]][[]], genotype = factor(sobj_sub$genotype, levels = c("WT", "KO")), clust = Idents(sobj_sub))
write.table(tab, file = "suppfig5d_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

############ Fig 5e (Bubble Plot) #######

# Get cell barcodes that remain after filtering (for use in following Python code)
WT <- subset(sobj_sub, subset = genotype == "WT")
KO <- subset(sobj_sub, subset = genotype == "KO")

WT_cells <- as.data.frame(Cells(WT)) %>% set_colnames(c("Barcode")) %>% 
    separate(col = "Barcode", into = c("Barcode", "index"), sep = "-") %>%
    mutate(Barcode = case_when(index == 1 ~ paste0(Barcode, "-1"), index == 3 ~ paste0(Barcode, "-2"))) %>%
    select(Barcode)
KO_cells <- as.data.frame(Cells(KO)) %>% set_colnames(c("Barcode")) %>% 
    separate(col = "Barcode", into = c("Barcode", "index"), sep = "-") %>%
    mutate(Barcode = case_when(index == 2 ~ paste0(Barcode, "-1"), index == 4 ~ paste0(Barcode, "-2"))) %>%
    select(Barcode)
write.table(WT_cells, file = "WT_Cells.csv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(KO_cells, file = "KO_Cells.csv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Python code to intersect the WT and KO specific peaks (from called peaks, not stencil peaks) with MCF7 (Estrogen & Control) ReMap peaks, and generate bubble plot
# cd /mBE/
# rm -rf BubblePlotBase
# mkdir BubblePlotBase
# mkdir BubblePlotBase/WT
# mkdir BubblePlotBase/KO
# basedir=/mBE/MCF7/scATAC_hg38
# cp -r ${basedir}/aggrWT/outs/filtered_peak_bc_matrix/ /mBE/BubblePlotBase/WT/
# cp -r ${basedir}/aggrWT/outs/filtered_peak_bc_matrix.h5 /mBE/BubblePlotBase/WT/
# cp -r ${basedir}/aggrWT/outs/peaks.bed /mBE/BubblePlotBase/WT/
# cp -r ${basedir}/aggrKO/outs/filtered_peak_bc_matrix/ /mBE/BubblePlotBase/KO/
# cp -r ${basedir}/aggrKO/outs/filtered_peak_bc_matrix.h5 /mBE/BubblePlotBase/KO/
# cp -r ${basedir}/aggrKO/outs/peaks.bed /mBE/BubblePlotBase/KO/
# cp -r /mBE/MCF7/scATAC_hg38_stencil/WT_Cells.csv /mBE/BubblePlotBase/
# cp -r /mBE/MCF7/scATAC_hg38_stencil/KO_Cells.csv /mBE/BubblePlotBase/

# cd /mBE/
# python runBEDMergePlot2.py --buildStats --natCommCells --jointTFs --intersectE2Control --intersectPeaks --outfile bubblePlotStats.pkl --csvOut bubblePlotStats.csv
# # python runBEDMergePlot2.py --buildStats --natCommCells --intersectPeaks --outfile bubblePlotStats.pkl --csvOut bubblePlotStats.csv
# python runBEDMergePlot2.py --infile bubblePlotStats.pkl --leaveNegatives

############ Supp Fig 5a,b,c (QC plots) #############
myElbowPlot <- function(object, i, j, reduction = 'pca') {
  data.use <- Stdev(object = object, reduction = reduction)
  if (length(x = data.use) == 0) {
    stop(paste("No standard deviation info stored for", reduction))
  }
  if (j > length(x = data.use)) {
    warning("The object only has information for ", length(x = data.use), " reductions")
    j <- length(x = data.use)
  }
  stdev <- 'Standard Deviation'
  p <- ggplot(data = data.frame(dims = i:j, stdev = data.use[i:j])) +
    geom_point(mapping = aes_string(x = 'dims', y = 'stdev')) +
    labs(
      x = gsub(
        pattern = '_$',
        replacement = '',
        x = Key(object = object[[reduction]])
      ),
      y = stdev
    ) +
    theme_cowplot()
  # write.table(data.frame(dims = i:j, stdev = data.use[i:j]), file = "suppfig5c_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  return(list(p))
}

Idents(sobj_sub) <- ""
pdf(file='MCF7_scATAC_qc1.pdf', width=6, height=5)
p1 <- wrap_plots(lapply(c("Rep1_WT", "Rep2_WT", "Rep1_KO", "Rep2_KO"),  function(x){TSSPlot(subset(sobj_sub, sample == x)) + theme(axis.title = element_blank()) + NoLegend() + ggtitle(x)}), nrow = 1)
p2 <- wrap_plots(lapply(c("Rep1_WT", "Rep2_WT", "Rep1_KO", "Rep2_KO"),  function(x){FragmentHistogram(subset(sobj_sub, sample == x)) + theme(axis.title = element_blank()) + NoLegend() + ggtitle(x)}), nrow = 1) 
p1 / p2
dev.off()


md <- sobj_sub@meta.data %>% mutate(genotype = factor(genotype, levels = c("WT", "KO")))
pdf("MCF7_scATAC_qc_by_cluster.pdf", width = 7, height = 6)
p1 <- ggplot(md, aes(x = seurat_clusters, y = total)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))
p2 <- ggplot(md, aes(x = seurat_clusters, y = pct_reads_in_peaks)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))
p3 <- ggplot(md, aes(x = seurat_clusters, y = peak_region_fragments)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))
p4 <- ggplot(md, aes(x = seurat_clusters, y = TSS.enrichment)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))
p5 <- ggplot(md, aes(x = seurat_clusters, y = nucleosome_signal)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))

p6 <- ggplot(md, aes(x = genotype, y = total)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))
p7 <- ggplot(md, aes(x = genotype, y = pct_reads_in_peaks)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))
p8 <- ggplot(md, aes(x = genotype, y = peak_region_fragments)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))
p9 <- ggplot(md, aes(x = genotype, y = TSS.enrichment)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))
p10 <- ggplot(md, aes(x = genotype, y = nucleosome_signal)) + geom_boxplot(outlier.size = 0.1) + stat_compare_means(aes(label = ..p.signif..))

plot(((p1 | p2 | p3 | p4 | p5) / (p6 | p7 | p8 | p9 | p10)) & theme(legend.position = "None", axis.title.x = element_blank()))
dev.off()

md %<>% select(seurat_clusters, genotype, total, pct_reads_in_peaks, peak_region_fragments, TSS.enrichment, nucleosome_signal)
setwd("/data/")
write.table(md, file = "suppfig5b_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

pdf(file='MCF7_scATAC_qc2.pdf', width=4.5, height=5)
myElbowPlot(sobj_sub, i = 2, j = 30, reduction = "lsi")
dev.off()

##### Cutsite distributions (Fig 5d), Cicero (Supp Fig 5e) #############

WT <- subset(sobj_sub, subset = genotype == "WT")
KO <- subset(sobj_sub, subset = genotype == "KO")

WT.counts <- GetAssayData(WT[["ATAC"]], slot = "counts")
KO.counts <- GetAssayData(KO[["ATAC"]], slot = "counts")

peaks <- read.table("/mBE/MCF7/cCREs/classified_hg38_4col.bed") %>% set_colnames(c("chr", "start", "end", "newclass"))

filtered_peaks_WT <- data.frame(rname=rownames(WT.counts)) %>% separate(rname, c("chr", "start", "end"), convert = TRUE) %>% left_join(peaks)
Active_len <- filtered_peaks_WT %>% filter(newclass == "Active") %>% mutate(length = end-start) %>% select(length) %>% sum()
H3K4me3_len <- filtered_peaks_WT %>% filter(newclass == "H3K4me3") %>% mutate(length = end-start) %>% select(length) %>% sum()
ATAConly_len <- filtered_peaks_WT %>% filter(newclass == "ATAConly") %>% mutate(length = end-start) %>% select(length) %>% sum()
Inactive_len <- filtered_peaks_WT %>% filter(newclass == "Inactive") %>% mutate(length = end-start) %>% select(length) %>% sum()
mBE_len <- filtered_peaks_WT %>% filter(newclass == "mBE") %>% mutate(length = end-start) %>% select(length) %>% sum()
WT.sums <- t(aggregate.Matrix(WT.counts, filtered_peaks_WT$newclass)) %>%
    data.frame() %>%
    mutate(Active = (Active * 1000000) / Active_len) %>%
    mutate(H3K4me3 = (H3K4me3 * 1000000) / H3K4me3_len) %>%
    mutate(ATAConly = (ATAConly * 1000000) / ATAConly_len) %>%
    mutate(Inactive = (Inactive * 1000000) / Inactive_len) %>%
    mutate(mBE = (mBE * 1000000) / mBE_len) %>%
    pivot_longer(everything(), names_to = "newclass", values_to = "cutsites") %>%
    mutate(newclass = factor(newclass, levels = c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE"))) %>%
    mutate(bio = "WT")

filtered_peaks_KO <- data.frame(rname=rownames(KO.counts)) %>% separate(rname, c("chr", "start", "end"), convert = TRUE) %>% left_join(peaks)
Active_len <- filtered_peaks_KO %>% filter(newclass == "Active") %>% mutate(length = end-start) %>% select(length) %>% sum()
H3K4me3_len <- filtered_peaks_KO %>% filter(newclass == "H3K4me3") %>% mutate(length = end-start) %>% select(length) %>% sum()
ATAConly_len <- filtered_peaks_KO %>% filter(newclass == "ATAConly") %>% mutate(length = end-start) %>% select(length) %>% sum()
Inactive_len <- filtered_peaks_KO %>% filter(newclass == "Inactive") %>% mutate(length = end-start) %>% select(length) %>% sum()
mBE_len <- filtered_peaks_KO %>% filter(newclass == "mBE") %>% mutate(length = end-start) %>% select(length) %>% sum()
KO.sums <- t(aggregate.Matrix(KO.counts, filtered_peaks_KO$newclass)) %>%
    data.frame() %>%
    mutate(Active = (Active * 1000000) / Active_len) %>%
    mutate(H3K4me3 = (H3K4me3 * 1000000) / H3K4me3_len) %>%
    mutate(ATAConly = (ATAConly * 1000000) / ATAConly_len) %>%
    mutate(Inactive = (Inactive * 1000000) / Inactive_len) %>%
    mutate(mBE = (mBE * 1000000) / mBE_len) %>%
    pivot_longer(everything(), names_to = "newclass", values_to = "cutsites") %>%
    mutate(newclass = factor(newclass, levels = c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE"))) %>%
    mutate(bio = "KO")

plotdata <- bind_rows(WT.sums, KO.sums) %>% mutate(bio = factor(bio, levels = c("WT", "KO")))
plotdata <- plotdata %>% bind_cols(plotdata %>% group_by(bio, newclass) %>% select(group_cols()) %>% unite("group", everything()))
plotdata <- plotdata %>% mutate(group = factor(group, levels = c("WT_Active", "KO_Active", "WT_H3K4me3", "KO_H3K4me3", "WT_ATAConly", "KO_ATAConly", "WT_Inactive", "KO_Inactive", "WT_mBE", "KO_mBE")))

my_comparisons <- list( c("WT_Active", "KO_Active"), c("WT_H3K4me3", "KO_H3K4me3"), c("WT_ATAConly", "KO_ATAConly"), c("WT_Inactive", "KO_Inactive"), c("WT_mBE", "KO_mBE") )

comp1 <- list( c("WT_Active", "KO_Active"),  c("WT_mBE", "KO_mBE") )
comp2 <- list( c("WT_H3K4me3", "KO_H3K4me3"), c("WT_ATAConly", "KO_ATAConly"), c("WT_Inactive", "KO_Inactive") )

means <- aggregate(cutsites ~  group, plotdata, mean) %>%
    rename(means = cutsites) %>%
    left_join(aggregate(cutsites ~  group, plotdata, max)) %>%
    rename(maxs = cutsites) %>%
    mutate(bio = rep(c("WT", "KO"), 5))
means1 <- means %>% filter(group %in% c("WT_Active", "KO_Active", "WT_mBE", "KO_mBE"))
means2 <- means %>% filter(!(group %in% c("WT_Active", "KO_Active", "WT_mBE", "KO_mBE")))

# auto_cicero_MCF7.R (Running this separately because it's very time consuming; Run this before continuing)
conns.list <- readRDS("conns.list.rds")
WT_conns <- conns.list[["WT"]] %>% drop_na(coaccess) %>% filter(coaccess > 0)
KO_conns <- conns.list[["KO"]] %>% drop_na(coaccess) %>% filter(coaccess > 0)
conns_per_peak <- WT_conns %>% count(Peak1) %>% mutate(dataset = "WT") %>%
    bind_rows(KO_conns %>% count(Peak1) %>% mutate(dataset = "KO"))
mu <- ddply(conns_per_peak, "dataset", summarise, grp.mean=mean(n))

# Plot all violins
p1 <- ggplot(plotdata %>% filter(group %in% c("WT_Active", "KO_Active", "WT_mBE", "KO_mBE")), aes(group, cutsites, fill = bio)) + 
    geom_violin(draw_quantiles = c(0.5)) + 
    stat_compare_means(comparisons = comp1) +
    # geom_jitter(size = 0.3, stroke = 0, shape = 16) +
	scale_fill_manual(values=c("#F68824", "#8F4B9D")) + 
    theme_cowplot() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none", axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    geom_text(data = means1, aes(label = sprintf("%0.1f", round(means, digits = 1)), y = means + 30))
p2 <- ggplot(plotdata %>% filter(!(group %in% c("WT_Active", "KO_Active", "WT_mBE", "KO_mBE"))), aes(group, cutsites, fill = bio)) + 
    geom_violin(draw_quantiles = c(0.5)) + 
    stat_compare_means(comparisons = comp2) +
    # geom_jitter(size = 0.3, stroke = 0, shape = 16) +
	scale_fill_manual(values=c("#F68824", "#8F4B9D")) + 
    theme_cowplot() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none", axis.ticks.x = element_blank(), axis.title = element_blank()) +
    geom_text(data = means2, aes(label = sprintf("%0.1f", round(means, digits = 1)), y = means + 30))
p3 <- ggplot(conns_per_peak %>% mutate(dataset = factor(dataset, levels = c("WT", "KO"))), aes(x=dataset, y=n, fill=dataset)) + 
    geom_violin(draw_quantiles = c(0.5)) + 
    stat_compare_means() +
    # scale_fill_manual(values=c("#F5C168", "#7D6CB9")) + 
    scale_fill_manual(values=c("#F68824", "#8F4B9D")) + 
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5), axis.title.x = element_blank()) +
    geom_text(data = mu, aes(label = sprintf("%0.1f", round(grp.mean, digits = 1)), y = grp.mean + 3))
pdf("MCF7_scATAC_violins.pdf", width=8, height=5)
plot((p1 | p2 | p3) + plot_layout(widths = c(3, 3, 2)))
dev.off()

setwd("/data/")
write.table(plotdata, file = "fig5d_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(conns_per_peak, file = "suppfig5e_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########## Interaction analysis ##############

setwd("/mBE/MCF7/scATAC_hg38_stencil/")
conns.list <- readRDS("conns.list.rds")
WT_conns <- conns.list[["WT"]] %>% drop_na(coaccess) %>% filter(coaccess > 0.1) %>%
    separate(col = Peak1, sep = "-", into = c(NA, "start1", NA), remove = FALSE, convert = TRUE) %>%
    separate(col = Peak2, sep = "-", into = c(NA, "start2", NA), remove = FALSE, convert = TRUE)
WT_conns_1 <- WT_conns %>%
    filter(start2 < start1) %>%
    mutate(Peak1t = Peak2, Peak2 = Peak1) %>%
    select(Peak1 = Peak1t, Peak2, coaccess) %>%
    bind_rows(WT_conns %>%
    filter(start2 > start1) %>%
    select(Peak1, Peak2, coaccess)) %>%
    distinct()

KO_conns <- conns.list[["KO"]] %>% drop_na(coaccess) %>% filter(coaccess > 0.1) %>%
    separate(col = Peak1, sep = "-", into = c(NA, "start1", NA), remove = FALSE, convert = TRUE) %>%
    separate(col = Peak2, sep = "-", into = c(NA, "start2", NA), remove = FALSE, convert = TRUE)
KO_conns_1 <- KO_conns %>%
    filter(start2 < start1) %>%
    mutate(Peak1t = Peak2, Peak2 = Peak1) %>%
    select(Peak1 = Peak1t, Peak2, coaccess) %>%
    bind_rows(KO_conns %>%
    filter(start2 > start1) %>%
    select(Peak1, Peak2, coaccess)) %>%
    distinct()

peaks <- read.table("/mBE/MCF7/cCREs/classified_liftedto_hg38_clean_sorted.bed") %>%
    select(V1, V2, V3, V7) %>%
    set_colnames(c("chr", "start", "end", "class")) %>%
    mutate(start = as.character(start), end = as.character(end))

WT_conns_2 <- WT_conns_1 %>% 
    separate(col = Peak1, sep = "-", into = c("chr", "start", "end"), remove = TRUE) %>%
    left_join(peaks) %>% rename(c1 = class) %>% rename(chr1 = chr, start1 = start, end1 = end) %>%
    separate(col = Peak2, sep = "-", into = c("chr", "start", "end"), remove = TRUE) %>%
    left_join(peaks) %>% rename(c2 = class) %>% rename(chr2 = chr, start2 = start, end2 = end)
WT_conns_3 <- WT_conns_2 %>% select(-c(c1, c2)) %>% bind_cols(as.data.frame(t(apply(WT_conns_2 %>% select(c(c1, c2)), 1, sort)))
)
KO_conns_2 <- KO_conns_1 %>% 
    separate(col = Peak1, sep = "-", into = c("chr", "start", "end"), remove = TRUE) %>%
    left_join(peaks) %>% rename(c1 = class) %>% rename(chr1 = chr, start1 = start, end1 = end) %>%
    separate(col = Peak2, sep = "-", into = c("chr", "start", "end"), remove = TRUE) %>%
    left_join(peaks) %>% rename(c2 = class) %>% rename(chr2 = chr, start2 = start, end2 = end)
KO_conns_3 <- KO_conns_2 %>% select(-c(c1, c2)) %>% bind_cols(as.data.frame(t(apply(KO_conns_2 %>% select(c(c1, c2)), 1, sort)))
)

############# Fig 5f (Circos - Chord Diagram) ################

summ <- plyr::count(WT_conns_3 %>% select(V1, V2)) %>% rename(WT = freq) %>%
    left_join(plyr::count(KO_conns_3 %>% select(V1, V2)) %>% rename(KO = freq)) %>%
    pivot_longer(c(WT, KO), names_to = "genotype", values_to = "count") %>%
    rename(from = V2, to = V1) %>% mutate(from = factor(from, levels = c("mBE", "H3K4me3", "Active", "Inactive", "ATAConly"))) %>%
    mutate(to = factor(to, levels = c("mBE", "H3K4me3", "Active", "Inactive", "ATAConly")))

library(circlize)
pdf("MCF7_scATAC_cicero_circos.pdf", width=4, height=4)
grid.col = c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")
temp <- chordDiagram(summ %>% select(from, to, count) %>% set_colnames(c("from", "to", "value")), grid.col = grid.col, transparency = ifelse(summ$genotype == "KO", 0.2, 0.7))
dev.off()

setwd("/data/")
write.table(summ, file = "fig5f_data1.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

myclasses <- c("Active", "H3K4me3", "ATAConly", "Inactive", "mBE")
summ <- plyr::count(WT_conns_3 %>% select(V1, V2)) %>% rename(WT = freq) %>%
    left_join(plyr::count(KO_conns_3 %>% select(V1, V2)) %>% rename(KO = freq)) %>%
    unite("type", V1:V2, remove = TRUE, sep = "-") %>%
    mutate(type = case_when(type == "ATAConly-Active" ~ "Active-ATAConly", type == "ATAConly-H3K4me3" ~ "H3K4me3-ATAConly", type == "ATAConly-Inactive" ~ "Inactive-ATAConly", type == "ATAConly-mBE" ~ "mBE-ATAConly", TRUE ~ type)) %>%
    mutate(type = factor(type, levels = c("Active-Active", "Active-H3K4me3", "Active-Inactive", "Active-mBE", "Active-ATAConly", "H3K4me3-H3K4me3", "H3K4me3-Inactive", "H3K4me3-mBE", "H3K4me3-ATAConly", "Inactive-Inactive", "Inactive-mBE", "Inactive-ATAConly", "mBE-mBE", "mBE-ATAConly", "ATAConly-ATAConly"))) %>%
    mutate(FC = KO/WT) %>% mutate(type = fct_reorder(type, FC))

pdf("MCF7_scATAC_cicero_FC.pdf", width=4, height=5)
ggplot(summ, aes(y = type, x = FC)) +
    geom_col(position = position_dodge()) +
    scale_x_continuous(expand = expansion(mult = c(0, .1))) +
    theme_cowplot() + theme(axis.title.y = element_blank())
dev.off()

setwd("/data/")
write.table(summ, file = "fig5f_data2.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

############# Supp Fig 6f (Comparison of interactions; WT vs KO) ################

WT_conns_4 <- WT_conns_2 %>% 
    filter((c1 == "H3K4me3" & c2 != "H3K4me3") | (c1 != "H3K4me3" & c2 == "H3K4me3")) %>%
    unite("bing", chr1:end1, sep = "-") %>% unite("bong", chr2:end2, sep = "-") %>% 
    unite("intertype", c1:c2, sep = "-") %>% select(-coaccess)
WT_conns_5 <- WT_conns_4 %>% filter(intertype %in% c("H3K4me3-Inactive", "H3K4me3-mBE", "H3K4me3-Active", "H3K4me3-ATAConly")) %>%
    bind_rows(WT_conns_4 %>% filter(intertype == "ATAConly-H3K4me3") %>%
    mutate(temp = bing, bing = bong) %>% mutate(bong = temp) %>% select(-temp) %>% mutate(intertype = "H3K4me3-ATAConly")) %>%
    bind_rows(WT_conns_4 %>% filter(intertype == "Active-H3K4me3") %>%
    mutate(temp = bing, bing = bong) %>% mutate(bong = temp) %>% select(-temp) %>% mutate(intertype = "H3K4me3-Active")) %>%
    bind_rows(WT_conns_4 %>% filter(intertype == "Inactive-H3K4me3") %>%
    mutate(temp = bing, bing = bong) %>% mutate(bong = temp) %>% select(-temp) %>% mutate(intertype = "H3K4me3-Inactive")) %>%
    bind_rows(WT_conns_4 %>% filter(intertype == "mBE-H3K4me3") %>%
    mutate(temp = bing, bing = bong) %>% mutate(bong = temp) %>% select(-temp) %>% mutate(intertype = "H3K4me3-mBE"))

KO_conns_4 <- KO_conns_2 %>% 
    filter((c1 == "H3K4me3" & c2 != "H3K4me3") | (c1 != "H3K4me3" & c2 == "H3K4me3")) %>%
    unite("bing", chr1:end1, sep = "-") %>% unite("bong", chr2:end2, sep = "-") %>% 
    unite("intertype", c1:c2, sep = "-") %>% select(-coaccess)
KO_conns_5 <- KO_conns_4 %>% filter(intertype %in% c("H3K4me3-Inactive", "H3K4me3-mBE", "H3K4me3-Active", "H3K4me3-ATAConly")) %>%
    bind_rows(KO_conns_4 %>% filter(intertype == "ATAConly-H3K4me3") %>%
    mutate(temp = bing, bing = bong) %>% mutate(bong = temp) %>% select(-temp) %>% mutate(intertype = "H3K4me3-ATAConly")) %>%
    bind_rows(KO_conns_4 %>% filter(intertype == "Active-H3K4me3") %>%
    mutate(temp = bing, bing = bong) %>% mutate(bong = temp) %>% select(-temp) %>% mutate(intertype = "H3K4me3-Active")) %>%
    bind_rows(KO_conns_4 %>% filter(intertype == "Inactive-H3K4me3") %>%
    mutate(temp = bing, bing = bong) %>% mutate(bong = temp) %>% select(-temp) %>% mutate(intertype = "H3K4me3-Inactive")) %>%
    bind_rows(KO_conns_4 %>% filter(intertype == "mBE-H3K4me3") %>%
    mutate(temp = bing, bing = bong) %>% mutate(bong = temp) %>% select(-temp) %>% mutate(intertype = "H3K4me3-mBE"))

conns_new <- setdiff(WT_conns_5, KO_conns_5) %>% mutate(genotype = "WTonly") %>%
    bind_rows(intersect(WT_conns_5, KO_conns_5) %>% mutate(genotype = "both")) %>%
    bind_rows(setdiff(KO_conns_5, WT_conns_5) %>% mutate(genotype = "KOonly"))
conns_final <- conns_new %>% bind_cols(ClosestFeature(sobj_sub, conns_new$bing) %>% select(gene_name, distance))

summ <- conns_final %>% filter(distance == 0) %>% select(intertype, genotype) %>% table() %>% data.frame() %>% mutate(genotype = factor(genotype, levels = rev(c("WTonly", "both", "KOonly")))) %>% separate("intertype", into = c(NA, "newclass"), sep = "-") %>% mutate(newclass = factor(newclass, levels = c("Active", "ATAConly", "Inactive", "mBE")))

pdf("MCF7_scATAC_interactions_summary.pdf", width=3, height=4)
ggplot(summ, aes(y = genotype, x = Freq, fill = genotype)) +
    geom_col(color = "black") +
    facet_grid(rows = vars(newclass), switch="both") +
    scale_fill_grey() +
    scale_x_continuous(expand = expansion(mult = c(0, .1))) +
    theme_cowplot() + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

setwd("/data/")
write.table(summ, file = "suppfig5f_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

######### Supp Fig 6g (Gene Ontology analysis of promoter-enhancer interactions) ###########

summ2 <- conns_final %>% filter(genotype == "KOonly", intertype == "H3K4me3-mBE", distance == 0) %>% select(gene_name) %>% 
    unlist() %>% unname() %>% table() %>% data.frame() %>% arrange(Freq) %>% set_colnames(c("gene", "mBE")) %>%
    full_join(conns_final %>% filter(genotype == "KOonly", intertype == "H3K4me3-Active", distance == 0) %>% select(gene_name) %>% 
    unlist() %>% unname() %>% table() %>% data.frame() %>% arrange(Freq) %>% set_colnames(c("gene", "Active"))) %>%
    full_join(conns_final %>% filter(genotype == "KOonly", intertype == "H3K4me3-Inactive", distance == 0) %>% select(gene_name) %>% 
    unlist() %>% unname() %>% table() %>% data.frame() %>% arrange(Freq) %>% set_colnames(c("gene", "Inactive"))) %>%
    full_join(conns_final %>% filter(genotype == "KOonly", intertype == "H3K4me3-ATAConly", distance == 0) %>% select(gene_name) %>% 
    unlist() %>% unname() %>% table() %>% data.frame() %>% arrange(Freq) %>% set_colnames(c("gene", "ATAConly"))) %>%
    replace_na(list(mBE = 0, Active = 0, Inactive = 0, ATAConly = 0)) %>% mutate(total=Reduce("+",.[2:5])) %>% arrange(-total) %>%
    head(n = 25) %>% select(-total) %>% pivot_longer(!gene, names_to = "newclass", values_to = "count")

library(gprofiler2)
gost_results <- gost(conns_final %>% filter(genotype == "KOonly", intertype == "H3K4me3-mBE", distance == 0) %>% select(gene_name) %>% unlist() %>% unname() %>% unique())
gost_mBE <- gost_results$result %>% 
    filter(significant == TRUE, source == "GO:BP") %>%
    select(term_name, p_value, intersection_size, term_size, query_size) %>%
    slice_min(p_value, n = 6)
gost_results <- gost(conns_final %>% filter(genotype == "KOonly", intertype == "H3K4me3-Active", distance == 0) %>% select(gene_name) %>% unlist() %>% unname() %>% unique())
gost_Active <- gost_results$result %>% 
    filter(significant == TRUE, source == "GO:BP") %>%
    select(term_name, p_value, intersection_size, term_size, query_size) %>%
    slice_min(p_value, n = 6)
gost_results <- gost(conns_final %>% filter(genotype == "KOonly", intertype == "H3K4me3-Inactive", distance == 0) %>% select(gene_name) %>% unlist() %>% unname() %>% unique())
gost_Inactive <- gost_results$result %>% 
    filter(significant == TRUE, source == "GO:BP") %>%
    select(term_name, p_value, intersection_size, term_size, query_size) %>%
    slice_min(p_value, n = 6)
gost_results <- gost(conns_final %>% filter(genotype == "KOonly", intertype == "H3K4me3-ATAConly", distance == 0) %>% select(gene_name) %>% unlist() %>% unname() %>% unique())
gost_ATAConly <- gost_results$result %>% 
    filter(significant == TRUE, source == "GO:BP") %>%
    select(term_name, p_value, intersection_size, term_size, query_size) %>%
    slice_min(p_value, n = 6)
gost_all <- gost_Active %>% mutate(newclass = "Active") %>%
    bind_rows(gost_ATAConly %>% mutate(newclass = "ATAConly")) %>%
    bind_rows(gost_Inactive %>% mutate(newclass = "Inactive")) %>%
    bind_rows(gost_mBE %>% mutate(newclass = "mBE")) %>%
    mutate(term_name = factor(term_name, levels = c("cellular response to nitrogen compound", "positive regulation of protein localization", "response to endogenous stimulus", "response to nitrogen compound", "response to organonitrogen compound", "response to oxygen-containing compound", "osteoblast development", "regulation of RNA metabolic process", "regulation of nucleobase-containing compound metabolic process", "regulation of primary metabolic process", "regulation of cellular metabolic process", "regulation of nitrogen compound metabolic process", "cellular developmental process", "positive regulation of cellular metabolic process", "positive regulation of metabolic process", "regulation of metabolic process", "cellular protein modification process", "macromolecule modification", "protein modification process", "mitochondrial fusion", "maintenance of cell number", "stem cell population maintenance")), newclass = factor(newclass, levels = c("Active", "ATAConly", "Inactive", "mBE")))

pdf("MCF7_scATAC_interactions_GO.pdf", width=7, height=6)
ggplot(gost_all, aes(x = newclass, y = term_name, color = -log10(p_value))) +
    geom_point(aes(size = intersection_size)) +
    scale_y_discrete(position = "right") +
    scale_color_viridis() +
    theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank())
dev.off()

gost_all <- gost_Active %>% mutate(newclass = "Active") %>%
    bind_rows(gost_ATAConly %>% mutate(newclass = "ATAConly")) %>%
    bind_rows(gost_Inactive %>% mutate(newclass = "Inactive")) %>%
    bind_rows(gost_mBE %>% mutate(newclass = "mBE")) %>%
    relocate(newclass, .before = term_name)
write.table(gost_all, file = "gost_results.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

gost_all <- read.table("gost_results.tsv", sep = "\t", header = TRUE)
setwd("/data/")
write.table(gost_all %>% mutate(neg_log10_pval = -log10(p_value)) %>% select(newclass, term_name, neg_log10_pval, intersection_size), file = "suppfig5g_data.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

######### Supp Fig 6h - interactions (UCSC Browser Tracks) ###########

WT_conns_4 <- WT_conns_2 %>% filter((c1 == "H3K4me3" & c2 != "H3K4me3") | (c1 != "H3K4me3" & c2 == "H3K4me3")) %>%
    mutate(start1 = as.numeric(start1), end1 = as.numeric(end1), start2 = as.numeric(start2), end2 = as.numeric(end2)) %>%
    select(-c(coaccess, c1, c2))
WT_conns_5 <- WT_conns_4 %>% filter(start2 > start1) %>%
    bind_rows(WT_conns_4 %>% filter(start2 < start1) %>%
    mutate(chrtemp = chr1, starttemp = start1, endtemp = end1, chr1 = chr2, start1 = start2, end1 = end2) %>%
    mutate(chr2 = chrtemp, start2 = starttemp, end2 = endtemp) %>% select(-c(chrtemp, starttemp, endtemp)))
KO_conns_4 <- KO_conns_2 %>% filter((c1 == "H3K4me3" & c2 != "H3K4me3") | (c1 != "H3K4me3" & c2 == "H3K4me3")) %>%
    mutate(start1 = as.numeric(start1), end1 = as.numeric(end1), start2 = as.numeric(start2), end2 = as.numeric(end2)) %>%
    select(-c(coaccess, c1, c2))
KO_conns_5 <- KO_conns_4 %>% filter(start2 > start1) %>%
    bind_rows(KO_conns_4 %>% filter(start2 < start1) %>%
    mutate(chrtemp = chr1, starttemp = start1, endtemp = end1, chr1 = chr2, start1 = start2, end1 = end2) %>%
    mutate(chr2 = chrtemp, start2 = starttemp, end2 = endtemp) %>% select(-c(chrtemp, starttemp, endtemp)))
conns_merged <- setdiff(WT_conns_5, KO_conns_5) %>% mutate(genotype = "WTonly") %>%
    bind_rows(intersect(WT_conns_5, KO_conns_5) %>% mutate(genotype = "both")) %>%
    bind_rows(setdiff(KO_conns_5, WT_conns_5) %>% mutate(genotype = "KOonly")) %>%
    mutate(chrom = chr1, chromStart = start1, chromEnd = end2, name = ".", score = 0, value = 1, exp = ".", color = case_when(genotype == "WTonly" ~ "#008000", genotype == "both" ~ "#989798", genotype == "KOonly" ~ "#FF0000"), sourceChrom = chr1, sourceStart = start1, sourceEnd = end1, sourceName = ".", sourceStrand = ".", targetChrom = chr2, targetStart = start2, targetEnd = end2, targetName = ".", targetStrand = ".") %>% select(-c(chr1:genotype))
write.table(conns_merged, file = "merged_interactions.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Add the following header lines to the tsv file, and add Custom Track to UCSC Genome Browser for visualization: 
# track type=interact name="MCF7 WT and KO" description="MCF7 WT and KO" maxHeightPixels=200:100:50 visibility=full
# browser position chr20:50,631,851-50,853,470

color_key <- c(Active = "#0f9448", H3K4me3 = "#e78ac3", ATAConly = "#E0AC69", mBE = "#2b598b", Inactive = "#f15a2b")
peaks_track <- peaks %>% mutate(chrom = chr, chromStart = start, chromEnd = end, name = ".", score = 0, strand = ".", thickStart = start, thickEnd = end, itemRgb = case_when(class == "Active" ~ "15,148,72", class == "H3K4me3" ~ "231,138,195", class == "ATAConly" ~ "224,172,105", class == "mBE" ~ "43,89,139", class == "Inactive" ~ "241,90,43")) %>% select(-c(chr:class))
write.table(peaks_track, file = "peaks_track.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Add the following header lines to the tsv file, and add Custom Track to UCSC Genome Browser for visualization: 
# track name="MCF7 CRE" description="MCF7 CRE" visibility=2 itemRgb="On" visibility=dense
# browser position chr20:50,631,851-50,853,470

