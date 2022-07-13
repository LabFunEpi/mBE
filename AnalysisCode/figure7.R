######## Cell Ranger Runs #########

# crmult=scMULT_count_mouse.sh
# crmultaggr=scMULT_aggr_mouse.sh
# crout=/mBE/MULT

# qsub -N D92 -v id=D92,libraries=${crout}/D92.csv,outdir=${crout} ${crmult}
# qsub -N D93 -v id=D93,libraries=${crout}/D93.csv,outdir=${crout} ${crmult}
# qsub -N D92-2 -v id=D92-2,libraries=${crout}/D92-2.csv,outdir=${crout} ${crmult}
# qsub -N D93-2 -v id=D93-2,libraries=${crout}/D93-2.csv,outdir=${crout} ${crmult}
# qsub -N aggr -v id=aggr,libraries=${crout}/aggr.csv,outdir=${crout} ${crmultaggr}

######## Get Cell Ranger Summary ########
setwd("/mBE/MULT/")
samplenames <- c("D92", "D92-2", "D93", "D93-2")

cr_summary <- lapply(paste0(samplenames, "/outs/summary.csv"), function(x){read.csv(x) %>% t() %>% t() %>% data.frame()})
names(cr_summary) <- samplenames
cr_summary <- bind_rows(cr_summary) %>% t()

write.table(cr_summary, file = "summary.table", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(read.csv("aggr/outs/summary.csv") %>% t() %>% data.frame(), file = "aggr.summary.table", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)


#################################################################

setwd("/mBE/MULT/")
datadir <- "/mBE/MULT/"
samplenames <- c("D92", "D92-2", "D93", "D93-2")
library(ggallin)
library(scico)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "mm10"
# saveRDS(annotations, "/references/annotations_mm10.rds")
annotations <- readRDS("/references/annotations_mm10.rds")
theme_set(theme_cowplot())
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
keepBSgenomeSequences <- function(genome, seqnames)
{
    stopifnot(all(seqnames %in% seqnames(genome)))
    genome@user_seqnames <- setNames(seqnames, seqnames)
    genome@seqinfo <- genome@seqinfo[seqnames]
    genome
}
sequences_to_keep <- paste0("chr", c(1:19, "X", "Y"))
genome <- keepBSgenomeSequences(genome, sequences_to_keep)
library(BiocParallel)
register(SerialParam())


###################### Cell Ranger - Batch Correction (Aggr) ##################

h5 <- "/mBE/MULT/aggr/outs/filtered_feature_bc_matrix.h5"
frag.file <- "/mBE/MULT/aggr/outs/atac_fragments.tsv.gz"
metadata_csv_1 <- "/mBE/MULT/D92/outs/per_barcode_metrics.csv"
metadata_csv_2 <- "/mBE/MULT/D92-2/outs/per_barcode_metrics.csv"
metadata_csv_3 <- "/mBE/MULT/D93/outs/per_barcode_metrics.csv"
metadata_csv_4 <- "/mBE/MULT/D93-2/outs/per_barcode_metrics.csv"

inputdata.10x <- Read10X_h5(h5)
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
metadata_1 <- read.csv(file = metadata_csv_1, header = TRUE, row.names = 1) %>% filter(is_cell == 1) %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", NA), sep = "-") %>% mutate(barcode = paste0(barcode, "-1"))
metadata_2 <- read.csv(file = metadata_csv_2, header = TRUE, row.names = 1) %>% filter(is_cell == 1) %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", NA), sep = "-") %>% mutate(barcode = paste0(barcode, "-2"))
metadata_3 <- read.csv(file = metadata_csv_3, header = TRUE, row.names = 1) %>% filter(is_cell == 1) %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", NA), sep = "-") %>% mutate(barcode = paste0(barcode, "-3"))
metadata_4 <- read.csv(file = metadata_csv_4, header = TRUE, row.names = 1) %>% filter(is_cell == 1) %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", NA), sep = "-") %>% mutate(barcode = paste0(barcode, "-4"))
metadata <- data.frame(barcode = colnames(rna_counts)) %>% left_join(bind_rows(metadata_1, metadata_2, metadata_3, metadata_4)) %>% column_to_rownames("barcode")

sobj <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
sobj[["percent.mito"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = frag.file,
    min.cells = 10,
    annotation = annotations
)
sobj[["ATAC"]] <- chrom_assay
DefaultAssay(sobj) <- "ATAC"
sobj <- NucleosomeSignal(object = sobj)
sobj <- TSSEnrichment(object = sobj, fast = FALSE)

DefaultAssay(sobj) <- "RNA"
sobj <- sobj %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100) %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DefaultAssay(sobj) <- "ATAC"
sobj <- sobj %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") %>% FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50)) %>% RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>% FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE)

mapping <- data.frame(sampleno = c("1", "2", "3", "4"), sample = c("D92", "D92-2", "D93", "D93-2"))
all.metadata <- sobj@meta.data %>% rownames_to_column("barcode") %>% separate(col = "barcode", into = c("barcode", "sampleno"), sep = "-") %>% left_join(mapping)
sobj$sample <- all.metadata$sample
saveRDS(sobj, "aggr_seuratobj.rds")

sobj_filtered <- subset(sobj, subset = percent.mito < 35 & nFeature_RNA > 200 & nFeature_ATAC > 200 & TSS.enrichment > 1)
DefaultAssay(sobj_filtered) <- "RNA"
sobj_filtered <- sobj_filtered %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100) %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DefaultAssay(sobj_filtered) <- "ATAC"
sobj_filtered <- sobj_filtered %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") %>% FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50)) %>% RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>% FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE)


ncells <- min(table(sobj_filtered$sample))
sub_cells <- c(
    sample(sobj_filtered@meta.data %>% filter(sample == "D92") %>% rownames(), size = ncells, replace=F),
    sample(sobj_filtered@meta.data %>% filter(sample == "D92-2") %>% rownames(), size = ncells, replace=F),
    sample(sobj_filtered@meta.data %>% filter(sample == "D93") %>% rownames(), size = ncells, replace=F),
    sample(sobj_filtered@meta.data %>% filter(sample == "D93-2") %>% rownames(), size = ncells, replace=F)
)
sobj_sub <- sobj_filtered[,sub_cells]
DefaultAssay(sobj_sub) <- "RNA"
sobj_sub <- sobj_sub %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100) %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DefaultAssay(sobj_sub) <- "ATAC"
sobj_sub <- sobj_sub %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") %>% FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50)) %>% RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") 
sobj_sub <- sobj_sub %>% FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.1)
sobj_sub$clust0.1 <- sobj_sub$seurat_clusters
sobj_sub <- sobj_sub %>% FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.2)
sobj_sub$clust0.2 <- sobj_sub$seurat_clusters
sobj_sub <- sobj_sub %>% FindClusters(graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.12)
sobj_sub$clust0.12 <- sobj_sub$seurat_clusters

DefaultAssay(sobj_sub) <- "ATAC"
sobj_sub <- sobj_sub %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") %>% FindNeighbors(reduction = "lsi", dims = 2:50, annoy.metric = "cosine")
sobj_sub <- sobj_sub %>% FindClusters(algorithm = 3, verbose = FALSE, resolution = 0.1)
sobj_sub$clustatac0.1 <- sobj_sub$seurat_clusters

saveRDS(sobj_sub, "aggr_filtered_sub.rds")

##############

DefaultAssay(sobj_sub) <- "ATAC" # Adding assay for Mouse Motifs
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 10090, all_versions = TRUE))
sobj_sub <- AddMotifs(sobj_sub, genome = genome, pfm = pwm_set)
sobj_sub <- RunChromVAR(
  object = sobj_sub,
  genome = genome,
  assay = "ATAC"
)

sobj_sub[["hATAC"]] <- sobj_sub[["ATAC"]] # Adding assay for Human Motifs
DefaultAssay(sobj_sub) <- "hATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE)) # Adding Human Motifs
sobj_sub <- AddMotifs(sobj_sub, genome = genome, pfm = pwm_set)
sobj_sub <- RunChromVAR(
  object = sobj_sub,
  genome = genome,
  assay = "hATAC",
  new.assay.name = "chromvarhuman"
)


##################

#########

sobj_sub <- readRDS("aggr_filtered_sub.rds")

DefaultAssay(sobj_sub) <- "RNA"
sobj_sub$celltype <- factor(recode_factor(sobj_sub$clust0.1, !!!c(`1` = "Basal", `0` = "Luminal progenitors", `2` = "Luminal ER+")), levels = c("Basal", "Luminal progenitors", "Luminal ER+"))
saveRDS(sobj_sub, "aggr_filtered_sub.rds")

theme_set(theme_cowplot())
umaps <- function(obj){
    p1 <- ggplot(data.frame(obj[["umap.rna"]][[]], label = factor(obj$sample), clust = Idents(obj)), aes(x=rnaUMAP_1, y=rnaUMAP_2, color = label)) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16, show.legend = FALSE), dpi = 400)
    p2 <- ggplot(data.frame(obj[["umap.atac"]][[]], label = factor(obj$sample), clust = Idents(obj)), aes(x=atacUMAP_1, y=atacUMAP_2, color = label)) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16, show.legend = FALSE), dpi = 400)
    p3 <- ggplot(data.frame(obj[["wnn.umap"]][[]], label = factor(obj$sample), clust = Idents(obj)), aes(x=wnnUMAP_1, y=wnnUMAP_2, color = label)) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16, show.legend = TRUE), dpi = 400) +
        guides(color = guide_legend(override.aes = list(size=5)))
    return(p1 | p2 | p3)
}
umaps1 <- function(obj){
    p1 <- ggplot(data.frame(obj[["umap.rna"]][[]], label = factor(obj$sample), clust = Idents(obj)), aes(x=rnaUMAP_1, y=rnaUMAP_2, color = clust)) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16, show.legend = FALSE), dpi = 400)
    p2 <- ggplot(data.frame(obj[["umap.atac"]][[]], label = factor(obj$sample), clust = Idents(obj)), aes(x=atacUMAP_1, y=atacUMAP_2, color = clust)) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16, show.legend = FALSE), dpi = 400)
    p3 <- ggplot(data.frame(obj[["wnn.umap"]][[]], label = factor(obj$sample), clust = Idents(obj)), aes(x=wnnUMAP_1, y=wnnUMAP_2, color = clust)) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16, show.legend = TRUE), dpi = 400) +
        guides(color = guide_legend(override.aes = list(size=5)))
    return(p1 | p2 | p3)
}

celltype.table1 <- sobj_sub@meta.data %>% select(celltype, status) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=levels(sobj_sub$celltype)))
p2 <- ggplot(celltype.table1) + 
    geom_bar(aes(fill=celltype, x=count, y=factor(status, levels = c("KO", "WT"))), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_continuous(expand = expansion(mult = c(0, .08)), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    xlab("Cell type proportion") +
    theme_cowplot() + 
    theme(legend.position = "none", legend.title=element_blank(), axis.title.y = element_blank(), axis.line.x=element_blank())
pdf(file='MaSC_scMULT_umaps.pdf', width=12, height=8)
Idents(sobj_sub) <- "celltype"
(umaps(sobj_sub) / umaps1(sobj_sub) / p2) + plot_layout(heights = c(7, 7, 1), nrow = 3)
dev.off()

###########

sobj_sub <- readRDS("aggr_filtered_sub.rds")
DefaultAssay(sobj_sub) <- "RNA"
Idents(sobj_sub) <- "celltype"
types <- levels(Idents(sobj_sub))
sobj_sub$celltype.status <- paste(Idents(sobj_sub), sobj_sub$status, sep = ".")
Idents(sobj_sub) <- "celltype.status"
DefaultAssay(sobj_sub) <- "RNA"
statuses <- c("KO", "WT")
degs <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- levels(Idents(sobj_sub))
        if (!(paste0(celltype, ".", statuses[[1]]) %in% curridents) | !(paste0(celltype, ".", statuses[[2]]) %in% curridents)) return(data.frame())
        if (length(WhichCells(sobj_sub, idents = paste0(celltype, ".", statuses[[1]]))) < 3 | length(WhichCells(sobj_sub, idents = paste0(celltype, ".", statuses[[2]]))) < 3) return(data.frame())
        temp1 <- FindMarkers(sobj_sub, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(degs) <- types
degs <- bind_rows(degs, .id = "celltype")
degs <- degs %>% filter(p_val_adj < 0.05) %>% mutate(neg_log10_adj_pval = -log10(p_val_adj))
write.table(degs, file = "degs.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
degs <- read.table(file = "degs.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "neg_log10_adj_pval")) %>%
    mutate(celltype = factor(celltype, levels = levels(sobj_sub$celltype)))

mygenes <- c("Brd4", "Jarid2")
pdf(file = "MaSC_scMULT_degs.pdf", width = 8, height = 8)
p0 <- ggplot(degs, aes(x = avg_log2FC)) +
    geom_density(aes(color = celltype)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
p1 <- ggplot(degs, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
    geom_point(aes(color = celltype), size = 2, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    # annotate("text", label = "Adj. P value = 0.05", x = 3.5, y = -log10(0.05)+5, size = 3) +
    geom_label_repel(data = degs %>% mutate(gene = case_when(gene %in% mygenes ~ gene, TRUE ~ "")), aes(label=gene, color=celltype), size = 5, max.overlaps = Inf, show.legend = FALSE) +
    geom_text_repel(data = degs %>% mutate(gene = case_when(neg_log10_adj_pval > 75 & !(gene %in% mygenes) ~ gene, TRUE ~ "")), aes(label=gene, color=celltype), size = 4, max.overlaps = Inf, show.legend = FALSE) +
    # coord_cartesian(xlim = c(-4, 4)) +
    # scale_color_manual(values = colors.mn, aesthetics = c("color", "segment.color")) +
    labs(x = "log2 (Fold change of average expression)", y = "-log10 (Adjusted P value)") +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title=element_blank())
p3 <- ggplot(degs, aes(y = neg_log10_adj_pval)) +
    geom_density(aes(color = celltype)) +
    theme_cowplot() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
plot(p0 + guide_area() + p1  + plot_spacer() + plot_layout(nrow = 2, ncol = 2, widths = c(10, 2), heights = c(3, 15), guides = "collect") + plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))))
dev.off()

####### Cicero #################################################
# auto_cicero_MaSC.R  (Running this separately because it's very time consuming; Run this before continuing)
conns.list <- readRDS("conns.list.rds")
WT_conns <- conns.list[["WT"]] %>% drop_na(coaccess) %>% filter(coaccess > 0)
KO_conns <- conns.list[["KO"]] %>% drop_na(coaccess) %>% filter(coaccess > 0)
conns_per_peak <- WT_conns %>% count(Peak1) %>% mutate(dataset = "WT") %>%
    bind_rows(KO_conns %>% count(Peak1) %>% mutate(dataset = "KO"))
mu <- ddply(conns_per_peak, "dataset", summarise, grp.mean=mean(n))
p3 <- ggplot(conns_per_peak %>% mutate(dataset = factor(dataset, levels = c("WT", "KO"))), aes(x=dataset, y=n, fill=dataset)) + 
    geom_violin(draw_quantiles = c(0.5)) + 
    stat_compare_means() +
    # scale_fill_manual(values=c("#F5C168", "#7D6CB9")) + 
    scale_fill_manual(values=c("#F68824", "#8F4B9D")) + 
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5), axis.title.x = element_blank()) +
    geom_text(data = mu, aes(label = sprintf("%0.1f", round(grp.mean, digits = 1)), y = grp.mean + 3))

    # KO     WT
# 217727 199052
pdf("MaSC_scMULT_conns_violin.pdf", width=4, height=5)
p3
dev.off()

####################################

Idents(sobj_sub) <- "clust0.1"
types <- levels(Idents(sobj_sub))
sobj_sub$celltype.status <- paste(Idents(sobj_sub), sobj_sub$status, sep = ".")
Idents(sobj_sub) <- "celltype.status"

lapply(levels(Idents(sobj_sub)), function(id){
    obj_temp <- subset(sobj_sub, subset = celltype.status == id)
    n <- 10
    thresh <- (ncol(obj_temp) / 100) * n
    peaks <- data.frame(peak = rownames(obj_temp[["ATAC"]])[rowSums(obj_temp[["ATAC"]]@counts > 0) > thresh]) %>%
        separate(peak, sep = "-", into = c("chr", "start", "end"))
    write.table(peaks, file = paste0(id, ".bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
})


# ls *.WT.bed | xargs -I{} mv {} /mBE/bedFiles/mus3Clusters.{}
# ls *.KO.bed | xargs -I{} mv {} /mBE/bedFiles/mus3Clusters.{}
# cd /mBE/
# ./fisherTestMusB.py --extreme 32 --tfLimit 65 --pvalueThreshold 0.05 --nClusters 3 --tfBed remap2022_nr_macs2_mm10_v1_0.bed 1> Clust3.log 2> 2 &
# cp /mBE/fishers_3Clusters.csv /mBE/MULT/Fishers/

setwd("/mBE/MULT/Fishers")
generated_file <- "/mBE/MULT/Fishers/fishers_3Clusters.csv"
fishers <- read.csv(file = generated_file, header = TRUE) %>% 
    set_colnames(c("TF", "0_OR", "1_OR", "2_OR", "0_LOR", "1_LOR", "2_LOR", "0_pval", "1_pval", "2_pval")) %>%
    # set_colnames(c("TF", "0_OR", "1_OR", "2_OR", "3_OR", "0_LOR", "1_LOR", "2_LOR", "3_LOR", "0_pval", "1_pval", "2_pval", "3_pval")) %>%
    pivot_longer(!TF, names_sep = "_", names_to = c("celltype", "valuetype"), values_to = "value") %>% 
    pivot_wider(names_from = "valuetype", values_from = "value") %>%
    mutate(neg_log10_pval = -log10(pval), celltype = factor(celltype)) %>%
    filter(pval < 0.05) %>% 
    mutate(celltype = factor(recode_factor(celltype, !!!c(`1` = "Basal", `0` = "Luminal progenitors", `2` = "Luminal ER+")), levels = c("Basal", "Luminal progenitors", "Luminal ER+")))
    # mutate(celltype = factor(recode_factor(celltype, !!!c(`0` = "Luminal progenitor", `1` = "Basal", `2` = "Alveolar precursor", `3` = "Mature luminal")), levels = c("Basal", "Luminal progenitor", "Alveolar precursor", "Mature luminal")))

mygenes <- c("FOXA1", "SOX10", "TRPS1")
highlight <- c("BRD4")
pdf(file = "MaSC_scMULT_ReMap_Fishers.pdf", width = 8, height = 8)
p0 <- ggplot(fishers, aes(x = LOR)) +
    geom_density(aes(color = celltype)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
p1 <- ggplot(fishers, aes(x = LOR, y = neg_log10_pval)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    # annotate("text", label = "P value = 0.05", x = 0.5, y = -log10(0.05)+0.02, size = 3) +
    geom_label_repel(data = fishers %>% mutate(TF = case_when((TF %in% highlight) ~ TF, TRUE ~ "")), aes(label=TF, colour = celltype), size = 5, max.overlaps = Inf, show.legend = FALSE, force = 1) +
    # geom_label_repel(data = fishers %>% mutate(TF = case_when((neg_log10_pval > -log10(0.05)) ~ TF, TRUE ~ "")), aes(label=TF, colour = celltype), size = 5, max.overlaps = Inf, show.legend = FALSE, force = 1) +
    geom_text_repel(data = fishers %>% mutate(TF = case_when((neg_log10_pval > 10 & !(TF %in% highlight) | TF %in% mygenes) ~ TF, TRUE ~ "")), aes(label=TF, colour = celltype), size = 4, max.overlaps = Inf, show.legend = FALSE) +
    # coord_cartesian(xlim = c(-4, 4)) +
    # scale_color_manual(values = colors.mn, aesthetics = c("color", "segment.color")) +
    labs(x = "Log2 Odds Ratio", y = "-log10 (P value)") +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title=element_blank())
p3 <- ggplot(fishers, aes(y = neg_log10_pval)) +
    geom_density(aes(color = celltype)) +
    theme_cowplot() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
plot(p0 + guide_area() + p1 + plot_spacer() + plot_layout(nrow = 2, ncol = 2, widths = c(10, 2), heights = c(3, 15), guides = "collect") + plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))))
dev.off()

./fisherTestMusB.py --extreme 32 --tfLimit 65 --pvalueThreshold 0.05 --nClusters 3 --tfBed remap2022_nr_macs2_mm10_v1_0_minusTSS.bed 1> Clust3_minusTSS.log 2> 2 &
cp /mBE/fishers_3Clusters.csv /mBE/MULT/Fishers/fishers_3Clusters_minusTSS.csv

setwd("/mBE/MULT/Fishers")
generated_file <- "/mBE/MULT/Fishers/fishers_3Clusters_minusTSS.csv"
fishers <- read.csv(file = generated_file, header = TRUE) %>% 
    set_colnames(c("TF", "0_OR", "1_OR", "2_OR", "0_LOR", "1_LOR", "2_LOR", "0_pval", "1_pval", "2_pval")) %>%
    pivot_longer(!TF, names_sep = "_", names_to = c("celltype", "valuetype"), values_to = "value") %>% 
    pivot_wider(names_from = "valuetype", values_from = "value") %>%
    mutate(neg_log10_pval = -log10(pval), celltype = factor(celltype)) %>%
    filter(pval < 0.05) %>% 
    mutate(celltype = factor(recode_factor(celltype, !!!c(`1` = "Basal", `0` = "Luminal progenitors", `2` = "Luminal ER+")), levels = c("Basal", "Luminal progenitors", "Luminal ER+")))

mygenes <- c("FOXA1", "SOX10", "TRPS1")
highlight <- c("BRD4")
pdf(file = "MaSC_scMULT_ReMap_Fishers_minusTSS.pdf", width = 8, height = 8)
p0 <- ggplot(fishers, aes(x = LOR)) +
    geom_density(aes(color = celltype)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
p1 <- ggplot(fishers, aes(x = LOR, y = neg_log10_pval)) +
    geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    geom_label_repel(data = fishers %>% mutate(TF = case_when((TF %in% highlight) ~ TF, TRUE ~ "")), aes(label=TF, colour = celltype), size = 5, max.overlaps = Inf, show.legend = FALSE, force = 1) +
    geom_text_repel(data = fishers %>% mutate(TF = case_when((neg_log10_pval > 10 & !(TF %in% highlight) | TF %in% mygenes) ~ TF, TRUE ~ "")), aes(label=TF, colour = celltype), size = 4, max.overlaps = Inf, show.legend = FALSE) +
    labs(x = "Log2 Odds Ratio", y = "-log10 (P value)") +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title=element_blank())
p3 <- ggplot(fishers, aes(y = neg_log10_pval)) +
    geom_density(aes(color = celltype)) +
    theme_cowplot() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
plot(p0 + guide_area() + p1 + plot_spacer() + plot_layout(nrow = 2, ncol = 2, widths = c(10, 2), heights = c(3, 15), guides = "collect") + plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))))
dev.off()


#########################


setwd("/mBE/MULT/")
sobj_sub <- readRDS("aggr_filtered_sub.rds")

TFs_mouse <- c("Sox10", "Sox9", "Sox4", "Nf1", "Trp63", "Tead4", "Elf5", "Foxa1", "Fosl1", "Atf3", "Fosl2", "Batf", "Jun", "Grhl2", "Zmynd8")
TFs_human <- c("SOX10", "SOX9", "SOX4", "NFIX", "TP63", "TEAD4", "ELF5", "FOXA1", "FOSL1", "ATF3", "FOSL2", "BATF", "JUN", "GRHL2")
TFs_remap <- c("SOX10", "SOX9", "TEAD4", "ELF5", "FOXA1", "FOSL1", "ATF3", "FOSL2", "BATF", "JUN", "GRHL2", "ZMYND8")
DefaultAssay(sobj_sub) <- "hATAC"
motifs <- ConvertMotifID(sobj_sub, name = TFs_human)
DefaultAssay(sobj_sub) <- "ATAC"
motifsm <- ConvertMotifID(sobj_sub, id = colnames(Motifs(sobj_sub)@data))

DefaultAssay(sobj_sub) <- "RNA"
sobj_sub$celltype.status <- paste(sobj_sub$celltype, sobj_sub$status, sep = ".")
Idents(sobj_sub) <- "celltype.status"
types <- levels(sobj_sub$celltype)
statuses <- c("KO", "WT")
degs <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(sobj_sub))
        temp1 <- FindMarkers(sobj_sub, features = TFs_mouse, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(degs) <- types
degs <- bind_rows(degs, .id = "celltype")
degs <- degs %>% mutate(enrichment = case_when(avg_log2FC > 0 ~ -log10(p_val_adj), TRUE ~ log10(p_val_adj)), type_gene = paste0(celltype, "_", gene))
degs <- degs %>% bind_rows(
    expand.grid(types, TFs_mouse) %>% 
    set_colnames(c("celltype", "gene")) %>% 
    mutate(enrichment = 0, type_gene = paste0(celltype, "_", gene)) %>% 
    filter(!(type_gene %in% degs$type_gene)))
degs <- degs %>% mutate(gene = factor(gene, levels = rev(TFs_mouse)), celltype = factor(celltype, levels(sobj_sub$celltype)))

DefaultAssay(sobj_sub) <- 'chromvarhuman'
dams <- lapply(
    types, 
    FUN = function(celltype){
        curridents <- unique(Idents(sobj_sub))
        temp1 <- FindMarkers(sobj_sub, features = motifs, ident.1 = paste0(celltype, ".", statuses[[1]]), ident.2 = paste0(celltype, ".", statuses[[2]]), mean.fxn = rowMeans, fc.name = "avg_diff") %>% rownames_to_column("gene")
        if(is.na(temp1)){ return(data.frame()) } else { return(temp1) }
    }
)
names(dams) <- types
dams <- bind_rows(dams, .id = "celltype")
dams <- dams %>% mutate(enrichment = case_when(avg_diff > 0 ~ -log10(p_val_adj), TRUE ~ log10(p_val_adj)), type_gene = paste0(celltype, "_", gene))
dams <- dams %>% bind_rows(
    expand.grid(types, motifs) %>% 
    set_colnames(c("celltype", "gene")) %>% 
    mutate(enrichment = 0, type_gene = paste0(celltype, "_", gene)) %>% 
    filter(!(type_gene %in% dams$type_gene)))
dams <- dams %>% mutate(gene = factor(gene, levels = rev(motifs)), celltype = factor(celltype, levels(sobj_sub$celltype)))

pdf(file='MaSC_scMULT_selectTFs_WTvsKO.pdf', width=12, height=6)
degs1 <- degs %>% mutate(enrichment = case_when(enrichment > 0 ~ enrichment, TRUE ~ 0))
p1 <- ggplot(degs1, aes(x = celltype, y = gene)) +
    geom_tile(aes(fill = enrichment), color = "grey", size=0.5) +
    scale_fill_gradient(low = "white", high = "red") +
    guides(fill=guide_colorbar(ticks.colour = NA)) +
    theme_cowplot() + ggtitle("Gene Expression") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), plot.title = element_text(size = 12), axis.ticks.length=unit(0, "cm"))
degs1 <- degs %>% mutate(enrichment = case_when(enrichment < 0 ~ -(enrichment), TRUE ~ 0))
p2 <- ggplot(degs1, aes(x = celltype, y = gene)) +
    geom_tile(aes(fill = enrichment), color = "grey", size=0.5) +
    scale_fill_gradient(low = "white", high = "blue") +
    guides(fill=guide_colorbar(ticks.colour = NA)) +
    theme_cowplot() + ggtitle("Gene Expression") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), plot.title = element_text(size = 12), axis.ticks.length=unit(0, "cm"))
dams1 <- dams %>% mutate(enrichment = case_when(enrichment > 0 ~ enrichment, TRUE ~ 0))
p3 <- ggplot(dams1, aes(x = celltype, y = gene)) +
    geom_tile(aes(fill = enrichment), color = "grey", size=0.5) +
    scale_fill_gradient(low = "white", high = "red") +
    guides(fill=guide_colorbar(ticks.colour = NA)) +
    scale_y_discrete(breaks=rev(motifs), labels=rev(TFs_human)) +
    theme_cowplot() + ggtitle("Motif accessibility (Human)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), plot.title = element_text(size = 12), axis.ticks.length=unit(0, "cm"))
dams1 <- dams %>% mutate(enrichment = case_when(enrichment < 0 ~ -(enrichment), TRUE ~ 0))
p4 <- ggplot(dams1, aes(x = celltype, y = gene)) +
    geom_tile(aes(fill = enrichment), color = "grey", size=0.5) +
    scale_fill_gradient(low = "white", high = "blue") +
    guides(fill=guide_colorbar(ticks.colour = NA)) +
    scale_y_discrete(breaks=rev(motifs), labels=rev(TFs_human)) +
    theme_cowplot() + ggtitle("Motif accessibility (Human)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), plot.title = element_text(size = 12), axis.ticks.length=unit(0, "cm"))
(p1 | p2 | p3 | p4) + plot_annotation(caption = 'Red = Enriched in KO; Blue = Enriched in WT')
dev.off()


# pdf(file='MaSC_scMULT_selectTFs_WTvsKO.pdf', width=6, height=6)
# p1 <- ggplot(degs, aes(x = celltype, y = gene)) +
    # geom_tile(aes(fill = enrichment), color = "grey", size=0.5) +
    # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    # theme_cowplot() + ggtitle("Gene Expression") +
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), plot.title = element_text(size = 12), axis.ticks.length=unit(0, "cm"))
# p2 <- ggplot(dams, aes(x = celltype, y = gene)) +
    # geom_tile(aes(fill = enrichment), color = "grey", size=0.5) +
    # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    # scale_y_discrete(breaks=rev(motifs), labels=rev(TFs_human)) +
    # theme_cowplot() + ggtitle("Motif accessibility (Human)") +
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), plot.title = element_text(size = 12), axis.ticks.length=unit(0, "cm"), axis.text.y = element_text(face="italic"))
# (p1 | p2) + plot_annotation(
  # caption = 'Red = Enriched in KO; Blue = Enriched in WT'
# )
# dev.off()


##########################


DefaultAssay(sobj_sub) <- "RNA"
Idents(sobj_sub) <- "celltype"
all.markers.seurat <- FindAllMarkers(object = sobj_sub, only.pos = FALSE)
all.markers.seurat <- all.markers.seurat %>% filter(p_val_adj < 0.05)
top.markers <- all.markers.seurat %>% group_by(cluster) %>% slice_head(n = 5)
pdf(file='MaSC_scMULT_allmarkers.pdf', width=15, height=8)
plotmarkers <- c(unique(rev(top.markers$gene)))
p1 <- DotPlot(sobj_sub, assay = "RNA", features = plotmarkers, cluster.idents = FALSE, dot.scale = 7) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_radius(breaks = c(20, 40, 60, 80)) +
    coord_flip() +
    theme(axis.text.y = element_text(size=20), axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=20))
plotmarkers <- c("Krt14", "Itgb1", "Krt8", "Krt18", "Pgr", "Esr1", "Gata3", "Foxa1", "Ltf", "Csn1s2a", "Csn2", "Csn3", "Cd14", "Cd55", "Aldh1a3", "Axl")
p2 <- DotPlot(sobj_sub, assay = "RNA", features = rev(plotmarkers), cluster.idents = FALSE, dot.scale = 7) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_radius(breaks = c(20, 40, 60, 80)) +
    coord_flip() +
    theme(axis.text.y = element_text(size=20), axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=20))
plotmarkers <- c("Areg", "Elf5", "Krt19", "Csn2", "Acta2", "Krt14", "Cxcl4", "Myh11", "Jund", "Irx5", "Sox4", "Igfbp2", "Cd55")
p3 <- DotPlot(sobj_sub, assay = "RNA", features = rev(plotmarkers), cluster.idents = FALSE, dot.scale = 7) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_radius(breaks = c(20, 40, 60, 80)) +
    coord_flip() +
    theme(axis.text.y = element_text(size=20), axis.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=20))
p1 | p2 | p3
dev.off()

##########

pdf("MaSC_scMULT_markers_featurePlots.pdf", width=16, height=8)
DefaultAssay(sobj_sub) <- "RNA"
genes <- c("Krt14", "Axl", "Csn3", "Cd14", "Krt8", "Elf5", "Pgr", "Esr1")
FeaturePlot(sobj_sub, features = genes, reduction = 'wnn.umap', pt.size = 0.4, raster=FALSE, order = TRUE, ncol = 4) & scale_color_scico(palette = 'lajolla')
dev.off()

##############







