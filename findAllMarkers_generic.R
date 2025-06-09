
args = commandArgs(trailingOnly=TRUE)
wd <- args[1]
filtered_rds_filename <- args[2]
cluster_labels <- args[3]
assay <- args[4]
prefix <- args[5]

setwd(wd)
filtered <- readRDS(filtered_rds_filename)

Idents(filtered) <- cluster_labels
DefaultAssay(filtered) <- assay
markers <- FindAllMarkers(filtered, only.pos = TRUE)

saveRDS(markers, paste0(prefix, "_markers.rds"))
write.table(markers, file = paste0(prefix, "_markers.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

