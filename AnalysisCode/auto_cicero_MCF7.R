setwd("/mBE/MCF7/scATAC/")

filtered <- readRDS("aggr_filtered_sub.rds")

library(BSgenome.Hsapiens.UCSC.hg38)
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
genome_lengths <- seqlengths(genome)
genome.df <- data.frame("chr" = names(genome_lengths), "length" = genome_lengths)

statuses <- c("WT", "KO")

conns.list <- lapply(statuses, function(x){
    MULT_sub <- subset(filtered, subset = genotype == x)
    DefaultAssay(MULT_sub) <- "ATAC"
    cds <- as.cell_data_set(MULT_sub)
    cicero.obj <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$UMAP)
    conns <- run_cicero(cicero.obj, genomic_coords = genome.df, sample_num = 100)
    return(conns)
})
names(conns.list) <- statuses

saveRDS(conns.list, "conns.list.rds")

