
args = commandArgs(trailingOnly=TRUE)
wd <- args[1]
filtered_rds_filename <- args[2]
cluster_labels <- args[3]
reduction <- args[4]
prefix <- args[5]
assay <- if (length(args) < 6) "mATAC" else args[6];

setwd(wd)
filtered <- readRDS(filtered_rds_filename)

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

# library(cicero)
# library(monocle)
# detach(package:cicero, unload=TRUE)
# detach(package:monocle, unload=TRUE)
library(cicero)
library(monocle3)

Idents(filtered) <- cluster_labels
DefaultAssay(filtered) <- assay
cds <- as.cell_data_set(filtered)
cicero.obj <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)[[reduction]])
conns <- run_cicero(cicero.obj, genomic_coords = genome.df, sample_num = 100)

saveRDS(conns, paste0(prefix, "_conns.rds"))

