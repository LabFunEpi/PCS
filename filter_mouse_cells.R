source("init.R")
runid <- "Run10"
setwd("/working")

library(EnsDb.Mmusculus.v79)
annotations_mm10 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotations_mm10) <- paste0('chr', seqlevels(annotations_mm10))
genome(annotations_mm10) <- "mm10"

##########################################################################

make_mult_obj_mm10 <- function(sampledir){
    h5 <- paste0(sampledir, "/outs/filtered_feature_bc_matrix.h5")
    frag.file <- paste0(sampledir, "/outs/atac_fragments.tsv.gz")
    metadata_csv <- paste0(sampledir, "/outs/per_barcode_metrics.csv")
    
    inputdata.10x <- Read10X_h5(h5)
    rna_counts <- inputdata.10x$`Gene Expression`
    atac_counts <- inputdata.10x$Peaks
    metadata <- read.csv(
		file = metadata_csv,
		header = TRUE,
		row.names = 1
	)
	sobj <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
    
    grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]

    chrom_assay <- CreateChromatinAssay(
        counts = atac_counts,
        sep = c(":", "-"),
        genome = 'mm10',
        fragments = frag.file,
        annotation = annotations_mm10
    )
    sobj[["ATAC"]] <- chrom_assay
    saveRDS(sobj, paste0(sampledir, "/outs/seurat_object.rds"))
    return(sobj)
}
samplenames <- c("PH080_mm10", "PH235_C_mm10", "PH235_T_mm10", "PH27_C_mm10", "PH27_T_mm10", "PH626_C_mm10", "PH626_T_mm10")
objs <- lapply(paste0(datadir, samplenames), make_mult_obj_mm10)
names(objs) <- samplenames

# ########################## Merge objs ##############################

# Add sample names
add.samplenames <- function(obj, samplename){
    obj$sample <- samplename
    return(obj)
}
objs <- mapply(add.samplenames, objs, names(objs))
objs <- lapply(objs, function(x){RenameCells(x, add.cell.id = x$sample[[1]])})

print(paste0(Sys.time(), ":: Merging RNA part ... "))

objs <- lapply(objs, function(x){DefaultAssay(x) <- "RNA"; return(DietSeurat(x, assays = "RNA"))})
objs_merged <- merge(objs[[1]], y = objs[2:length(objs)], project = "OvCa")
objs_merged <- JoinLayers(objs_merged)
saveRDS(objs_merged, paste0(runid, "_merged_rna_mm10.rds"))

print(paste0(Sys.time(), ":: Creating ATAC object list ... "))

peakslist <- lapply(paste0(datadir, samplenames), FUN = function(x){
    peaks <- read.table(paste0(x, "/outs/atac_peaks.bed"), col.names = c("chr", "start", "end"));
    return(peaks)
})
allpeaks <- makeGRangesFromDataFrame(bind_rows(peakslist))
combined.peaks <- GenomicRanges::reduce(allpeaks)
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks <- combined.peaks[as.vector(seqnames(combined.peaks) %in% standardChromosomes(combined.peaks))]

plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

atacobjlist <- lapply(paste0(datadir, samplenames), FUN = function(x){
    md <- read.table(
        file = paste0(x, "/outs/per_barcode_metrics.csv"),
        stringsAsFactors = FALSE,
        sep = ",",
        header = TRUE,
        row.names = 1
    )[-1, ]
    md <- md %>% filter(is_cell == 1)
    fragments <- CreateFragmentObject(paste0(x, "/outs/atac_fragments.tsv.gz"), cells = rownames(md));
    counts <- FeatureMatrix(fragments = fragments, features = combined.peaks, cells = rownames(md));
    assay <- CreateChromatinAssay(counts, fragments = fragments)
    obj <- CreateSeuratObject(assay, assay = "mATAC")
    return(obj)
})
names(atacobjlist) <- samplenames
add.samplenames <- function(obj, samplename){
    obj$sample <- samplename
    return(obj)
}
atacobjlist <- mapply(add.samplenames, atacobjlist, names(atacobjlist))
atacobjlist <- lapply(atacobjlist, function(x){
    x <- x %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 20) %>% RunSVD();
    return(x)
})
atacobjlist <- lapply(atacobjlist, function(x){RenameCells(x, add.cell.id = x$sample[[1]])})

print(paste0(Sys.time(), ":: Merging ATAC objects ... "))

merged_atac <- merge(x = atacobjlist[[1]], y = atacobjlist[2:length(atacobjlist)])
saveRDS(merged_atac, paste0(runid, "_merged_atac_mm10.rds"))

mult_integrated <- objs_merged
mult_integrated[["mATAC"]] <- merged_atac[["mATAC"]]

Annotation(mult_integrated[["mATAC"]]) <- annotations_mm10
mult_integrated <- NucleosomeSignal(object = mult_integrated, assay = "mATAC")
mult_integrated <- TSSEnrichment(object = mult_integrated, assay = "mATAC")

mult_integrated[["percent.mt"]] <- PercentageFeatureSet(mult_integrated, pattern = "^mt-")
saveRDS(mult_integrated, paste0(runid, "_merged_mm10.rds"))

print(paste0(Sys.time(), ":: Applying QC filters and running WNN on filtered object ... "))
# mult_integrated <- readRDS(paste0(runid, "_merged_mm10.rds"))

pdf(file=paste0(runid, "_qc_mm10.pdf"), width=10, height=10)
VlnPlot(mult_integrated, features = c("percent.mt", "nFeature_RNA", "nFeature_mATAC", "TSS.enrichment"), ncol = 4, raster = TRUE, pt.size = 0)
dev.off()

filtered <- subset(mult_integrated, subset = percent.mt < 50 & nFeature_RNA > 200 & nFeature_mATAC > 200 & TSS.enrichment > 1)

DefaultAssay(filtered) <- "RNA"
filtered <- filtered %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA(npcs = 50)
filtered <- RunUMAP(filtered, reduction = "pca", dims = 1:50, reduction.name = "rnaumap", reduction.key = "rnaumap_")

DefaultAssay(filtered) <- "mATAC"
filtered <- filtered %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 20) %>% RunSVD()
filtered <- RunUMAP(filtered, reduction = 'lsi', dims = 2:50, reduction.name = "atacumap", reduction.key = "atacumap_")

filtered <- filtered %>% FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), prune.SNN = 0)
filtered <- RunUMAP(filtered, nn.name = "weighted.nn", reduction.name = "wnnumap", reduction.key = "wnnumap_", n.epochs = 500, n.neighbors = 30)

saveRDS(filtered, paste0(runid, "_filtered_mm10.rds"))
######################################################################

filtered_mm10 <- readRDS("Run10_filtered_mm10.rds")
filtered_mm10 <- subset(filtered_mm10, sample %in% setdiff(unique(filtered_mm10$sample), "PH080_mm10"))
filtered_hg38 <- readRDS("RunPDX1_filtered.rds")

DefaultAssay(filtered_hg38) <- "RNA"
filtered_hg38 <- filtered_hg38 %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA(npcs = 50)
filtered_hg38 <- RunUMAP(filtered_hg38, reduction = "pca", dims = 1:50, reduction.name = "rnaumap", reduction.key = "rnaumap_")

DefaultAssay(filtered_hg38) <- "mATAC"
filtered_hg38 <- filtered_hg38 %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 20) %>% RunSVD()
filtered_hg38 <- RunUMAP(filtered_hg38, reduction = 'lsi', dims = 2:50, reduction.name = "atacumap", reduction.key = "atacumap_")

filtered_hg38 <- filtered_hg38 %>% FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), prune.SNN = 0)
filtered_hg38 <- RunUMAP(filtered_hg38, nn.name = "weighted.nn", reduction.name = "wnnumap", reduction.key = "wnnumap_", n.epochs = 500, n.neighbors = 30)

filtered_hg38 %<>% FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:10, 2:10), k.nn = 5)
filtered_hg38 %<>% FindClusters(graph.name = "wsnn", resolution = 0.06, verbose = FALSE)

pdf(file=paste0(runid, "_umaps_PDX_wnnclust.pdf"), width=18, height=6)
umaps_wrapped(filtered_hg38, rightcolors = "seurat_clusters")
dev.off()

# DefaultAssay(filtered_hg38) <- "RNA"
# Idents(filtered_hg38) <- "seurat_clusters"
# markers <- FindMarkers(object = filtered_hg38, ident.1 = 4)

#########

MouseNonOrthologs <- read.table("MouseNonOrthologs.txt", sep="\t", skip=1) %>% pull(V1)
MouseOrthologs <- read.table("MouseOrthologs.txt", sep="\t", skip=1) %>% pull(V1)
MouseNonOrthologs <- intersect(MouseNonOrthologs, rownames(filtered_mm10[["RNA"]]))
MouseOrthologs <- intersect(MouseOrthologs, rownames(filtered_mm10[["RNA"]]))
MouseNonOrthologs <- setdiff(MouseNonOrthologs, MouseOrthologs)
MouseOrthologs <- setdiff(MouseOrthologs, MouseNonOrthologs)

HumanNonOrthologs <- read.table("HumanNonOrthologs.txt", sep="\t", skip=1) %>% pull(V1)
HumanOrthologs <- read.table("HumanOrthologs.txt", sep="\t", skip=1) %>% pull(V1)
HumanNonOrthologs <- intersect(HumanNonOrthologs, rownames(filtered_hg38[["RNA"]]))
HumanOrthologs <- intersect(HumanOrthologs, rownames(filtered_hg38[["RNA"]]))
HumanNonOrthologs <- setdiff(HumanNonOrthologs, HumanOrthologs)
HumanOrthologs <- setdiff(HumanOrthologs, HumanNonOrthologs)

filtered_mm10_cells <- data.frame(bc = colnames(filtered_mm10), sample = filtered_mm10$sample) %>%
    separate_wider_delim(cols="bc", delim="_", names=c(NA, NA, NA, "bc"), too_few = "align_end") %>%
    mutate(sample = str_sub(sample, start=1L, end=-6L)) %>%
    tidyr::unite("bc", c(sample, bc), sep="_") %>% pull(bc)

filtered_hg38_cells <- filtered_hg38@meta.data %>% rownames()

cells <- intersect(filtered_hg38_cells, filtered_mm10_cells)

filtered_mm10_data <- GetAssayData(filtered_mm10, assay = "RNA", layer = "data")
colnames(filtered_mm10_data) <- filtered_mm10_cells
filtered_mm10_data <- filtered_mm10_data[union(MouseOrthologs, MouseNonOrthologs), cells]
filtered_hg38_data <- GetAssayData(filtered_hg38, assay = "RNA", layer = "data")
filtered_hg38_data <- filtered_hg38_data[union(HumanOrthologs, HumanNonOrthologs), cells]

OrthologsData <- data.frame(
    H_genes = filtered_hg38_data %>% colSums(),
    M_genes = filtered_mm10_data %>% colSums(),
    H_Orthologs = filtered_hg38_data[HumanOrthologs,] %>% colSums(),
    M_Orthologs = filtered_mm10_data[MouseOrthologs,] %>% colSums(),
    H_NonOrthologs = filtered_hg38_data[HumanNonOrthologs,] %>% colSums(),
    M_NonOrthologs = filtered_mm10_data[MouseNonOrthologs,] %>% colSums())
OrthologsData %<>% mutate(cluster = factor(OrthologsData %>% dplyr::select(H_genes, M_genes) %>% kmeans(3) %$% cluster, levels = c(1, 2, 3)))
OrthologsData %<>% rownames_to_column("cell")

OrthologsData %<>% left_join(data.frame(cell = colnames(filtered_hg38), seurat_clusters = filtered_hg38$seurat_clusters))
# OrthologsData %<>% mutate(species = case_when(M_genes > H_genes ~ "Mouse", TRUE ~ "Human"))

pdf(file=paste0(runid, "_ortholog_qc_mm10.pdf"), width=12, height=10)
ggplot(OrthologsData, aes(x = H_genes, y = M_genes, color=seurat_clusters)) +
    ggrastr::rasterise(geom_point(size=1.5, stroke=0.1, shape=16), dpi = 400) +
    geom_density2d(color = "black") +
	scale_color_manual(values = getPalette(length(unique(filtered_hg38$seurat_clusters)))) +
	guides(color = guide_legend(override.aes = list(size=5))) +
	coord_cartesian(xlim = c(350, 5600), ylim = c(350, 5600)) +
    theme_cowplot()
p1 <- ggplot(OrthologsData, aes(x = H_Orthologs, y = H_NonOrthologs, color=seurat_clusters)) +
    ggrastr::rasterise(geom_point(size=0.4, stroke=0.1, shape=16), dpi = 400) +
    geom_density2d(color = "black") +
    theme_cowplot() + theme(legend.position = "none")
p2 <- ggplot(OrthologsData, aes(x = M_Orthologs, y = M_NonOrthologs, color=seurat_clusters)) +
    ggrastr::rasterise(geom_point(size=0.4, stroke=0.1, shape=16), dpi = 400) +
    geom_density2d(color = "black") +
    theme_cowplot() + theme(legend.position = "none")
p3 <- ggplot(OrthologsData, aes(x = H_Orthologs, y = M_Orthologs, color=seurat_clusters)) +
    ggrastr::rasterise(geom_point(size=0.4, stroke=0.1, shape=16), dpi = 400) +
    geom_density2d(color = "black") +
    theme_cowplot() + theme(legend.position = "none")
p4 <- ggplot(OrthologsData, aes(x = H_NonOrthologs, y = M_NonOrthologs, color=seurat_clusters)) +
    ggrastr::rasterise(geom_point(size=0.4, stroke=0.1, shape=16), dpi = 400) +
    geom_density2d(color = "black") +
    theme_cowplot() + theme(legend.position = "none")
wrap_plots(p1, p2, p3, p4, ncol=2)
dev.off()

M_cells <- OrthologsData %>% filter(seurat_clusters == 4) %>% pull(cell)
write.table(data.frame(cells=M_cells), "MouseCells.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

################# Filtering out Mouse cells ###################

M_cells <- read.table("MouseCells.tsv") %>% pull(V1)
H_cells <- setdiff(colnames(filtered_hg38), M_cells)

filtered_hg38$CellName <- colnames(filtered_hg38)
filtered <- subset(filtered_hg38, subset = CellName %in% H_cells)
