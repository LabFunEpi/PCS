############ Run Cell Ranger on each sample (example) ###############
# ref=refdata-cellranger-arc-GRCh38-2020-A-2.0.0
# libraries=PH1214.csv
# cellranger-arc count --id=PH1214 --reference=${ref} --libraries=${libraries} --localcores=16 --localmem=115
# ...

############ Process Cell Ranger output with Seurat and Signac ####################

source("init.R")

setwd("/working/")
runid <- "Run10"
datadir <- "/datadir/"
samplenames <- c("FT243", "FT244", "FT245", "FT281", "FT283", "PH1199N", "PH1230", "PH1239", "PH1243", "PH1302", "PH1214", "PH1238", "PH1252", "1303_R")
FT_samplenames <- c("FT243", "FT244", "FT245", "FT281", "FT283")
samplelabels <- c(`FT243` = "FT243", `FT244` = "FT244", `FT245` = "FT245", `FT281` = "FT281", `FT283` = "FT283", `PH1199N` = "PH1199", `PH1230` = "PH1230", `PH1239` = "PH1239", `PH1243` = "PH1243", `PH1302` = "PH1302", `PH1214` = "PH1214", `PH1238` = "PH1238", `PH1252` = "PH1252", `1303_R` = "PH1303")

# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# genome(annotations) <- "hg38"
# saveRDS(annotations, "annotations.rds")

annotations <- readRDS("annotations.rds")

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
pwm_set <- readRDS("/softwares/R_libs/FigR/data/cisBP_human_pfms_2021.rds")

colors.hpca_sub <- c(`T_cells` = "#99CCCC", `NK_cell` = "#E2E8A6", `B_cell` = "#D68645", `Monocyte` = "#B4AEE5", `Macrophage` = "#B81254", `DC` = "#67A5E5", `Neutrophils` = "#3B69B7", `Endothelial_cells` = "#ffd92f", `Fibroblasts` = "#48B758", `Epithelial_cells` = "#e78ac3")
epicolors <- c(`Ciliated` = "#E41A1C", `Secretory` = "#449B75", `Cluster1` = "#C66764", `Cluster2` = "#999999", `Cluster3` = "#AC5782", `Cycling` = "#FFE528")
colors.hpca_sub2 <- c(`T_cell:CD8+` = "#A9E0AB", `T_cell:CD4+` = "#665299", `T_cell:Other` = "#31c4e9", `NK_cell` = "#E2E8A6", `B_cell` = "#D68645", `Monocyte` = "#B4AEE5", `Macrophage` = "#B81254", `DC` = "#67A5E5", `Neutrophils` = "#3B69B7", `Endothelial_cells` = "#ffd92f", `Fibroblasts` = "#48B758", `Epithelial_cells` = "#e78ac3")
colors.customtype <- c(`T_cell:CD8+` = "#A9E0AB", `T_cell:CD4+` = "#665299", `T_cell:Other` = "#31c4e9", `NK_cell` = "#E2E8A6", `B_cell` = "#D68645", `Monocyte` = "#B4AEE5", `Macrophage` = "#B81254", `DC` = "#67A5E5", `Neutrophils` = "#3B69B7", `Endothelial_cells` = "#ffd92f", `Fibroblasts` = "#48B758", `Epithelial_cells` = "#e78ac3", `Tumor` = "#e01c1a")
colors.customtype2 <- c(`T_cell:CD8+` = "#A9E0AB", `T_cell:CD4+` = "#665299", `T_cell:Other` = "#31c4e9", `NK_cell` = "#E2E8A6", `B_cell` = "#D68645", `Monocyte` = "#B4AEE5", `Macrophage` = "#B81254", `DC` = "#67A5E5", `Neutrophils` = "#3B69B7", `Endothelial_cells` = "#ffd92f", `Fibroblasts` = "#48B758", `Epithelial_cells` = "#e78ac3")

colors.subgroup2 <- c(`FT` = "#D37536", `Naive` = "#40B5B1", `NACT` = "#2F5C66")
colors.subgroup3 <- c(`BRCAwt` = "#E7A06E", `BRCAmt` = "#B36935", `Sensitive` = "#50C878", `Resistant` = "#FE6F5E", `NACT` = "#347A7B")

sample_levels <- c("FT243", "FT244", "FT245", "FT281", "FT283", "PH1199", "PH1230", "PH1302", "PH1239", "PH1243", "PH1303", "PH1214", "PH1238", "PH1252")
subgroup2_levels <- c("FT", "Naive", "NACT")
subgroup3_levels <- c("BRCAwt", "BRCAmt", "Sensitive", "Resistant", "NACT")
infotable <- data.frame(
    sample = sample_levels,
	subgroup2 = factor(c("FT", "FT", "FT", "FT", "FT", "Naive", "Naive", "Naive", "Naive", "Naive", "Naive", "NACT", "NACT", "NACT"), levels = subgroup2_levels),
    subgroup3 = factor(c("BRCAwt", "BRCAwt", "BRCAwt", "BRCAmt", "BRCAmt", "Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Resistant", "NACT", "NACT", "NACT"), levels = subgroup3_levels)
)

# ############################ MAKING SEURAT OBJECTS - QC - CLUSTERING - CELLTYPING #######################

make_mult_obj <- function(sampledir){
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
        genome = 'hg38',
        fragments = frag.file,
        annotation = annotations
    )
    sobj[["ATAC"]] <- chrom_assay
    saveRDS(sobj, paste0(sampledir, "/outs/seurat_object.rds"))
    return(sobj)
}

# objs <- lapply(paste0(datadir, samplenames), make_mult_obj)
# names(objs) <- samplenames

objs <- lapply(paste0(datadir, samplenames, "/outs/seurat_object.rds"), readRDS)
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
saveRDS(objs_merged, paste0(runid, "_merged_rna.rds"))

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
saveRDS(merged_atac, paste0(runid, "_merged_atac.rds"))

mult_integrated <- objs_merged
mult_integrated[["mATAC"]] <- merged_atac[["mATAC"]]

Annotation(mult_integrated[["mATAC"]]) <- annotations
mult_integrated <- NucleosomeSignal(object = mult_integrated, assay = "mATAC")
mult_integrated <- TSSEnrichment(object = mult_integrated, assay = "mATAC")

saveRDS(mult_integrated, paste0(runid, "_merged.rds"))
mult_integrated <- readRDS(paste0(runid, "_merged.rds"))

print(paste0(Sys.time(), ":: Applying QC filters and running WNN on filtered object ... "))

################### Plotting QC #######################

bettercolors <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#000000')

snames <- unname(samplelabels)
group <- c(rep("FT", 5), rep("Naive", 5), rep("NACT", 3), "Naive")
grouptbl <- data.frame(sample = snames, group = factor(group, levels = c("FT", "Naive", "NACT")))
grouptbl %<>% mutate(sample = factor(sample, levels = grouptbl %>% arrange(group) %>% pull(sample))) 

plotdata <- data.frame(
    ident = "all", 
    nFeature_RNA = mult_integrated$nFeature_RNA, 
    nCount_RNA = mult_integrated$nCount_RNA,
    percent.mito = mult_integrated$percent.mt,
    nFeature_mATAC = mult_integrated$nFeature_mATAC, 
    nCount_mATAC = mult_integrated$nCount_mATAC,
    TSS.enrichment = mult_integrated$TSS.enrichment,
    nucleosome_signal = mult_integrated$nucleosome_signal,
    sample = factor(case_when(mult_integrated$sample == "1303_R" ~ "PH1303", mult_integrated$sample == "PH1199N" ~ "PH1199", TRUE ~ mult_integrated$sample), levels = levels(grouptbl$sample))
) %>% left_join(grouptbl)
# plotdata1 <- data.frame(ident = "all", cutsites = rowSums(mult_integrated[['mATAC']]@counts), cells = rowSums(mult_integrated[['mATAC']]@counts > 0))
saveRDS(plotdata, 'Run10_qc.rds')

forlegend <- ggplot(plotdata, aes(x=ident, y=nFeature_RNA)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.05)) + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, aes(color = sample)), dpi = 400) +
        scale_y_continuous(trans = pseudolog10_trans) +
        scale_color_manual(values = bettercolors) +
        geom_hline(yintercept = 200) +
        ylab("No. of unique genes detected (RNA)") +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot()
mylegend <- get_legend(forlegend)

p1 <- ggplot(plotdata, aes(x=sample, y=nFeature_RNA)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5)) + 
        # ggrastr::rasterise(geom_jitter(aes(color = sample), size = 0.1, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        geom_hline(yintercept = 200) +
        facet_grid(cols = vars(group), scales = "free", space = "free") +
		ylab("No. of unique genes detected (GEX)") +
        theme_cowplot()
p2 <- ggplot(plotdata, aes(x=sample, y=nCount_RNA)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5)) + 
        # ggrastr::rasterise(geom_jitter(aes(color = sample), size = 0.1, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        facet_grid(cols = vars(group), scales = "free", space = "free") +
		ylab("Total no. of molecules detected (GEX)") +
        theme_cowplot()
p3 <- ggplot(plotdata, aes(x=sample, y=nFeature_mATAC)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5)) + 
        # ggrastr::rasterise(geom_jitter(aes(color = sample), size = 0.1, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        geom_hline(yintercept = 200) +
        facet_grid(cols = vars(group), scales = "free", space = "free") +
		ylab("No. of unique peaks detected (ATAC)") +
        theme_cowplot()
p4 <- ggplot(plotdata, aes(x=sample, y=nCount_mATAC)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5)) + 
        # ggrastr::rasterise(geom_jitter(aes(color = sample), size = 0.1, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(trans = pseudolog10_trans, breaks = log_breaks(n = 10, base = 10), labels = label_number(accuracy = 1, big.mark = "")) +
        facet_grid(cols = vars(group), scales = "free", space = "free") +
		ylab("Total no. of cutsites detected (ATAC)") +
        theme_cowplot()
p5 <- ggplot(plotdata, aes(x=sample, y=percent.mito)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5)) + 
        # ggrastr::rasterise(geom_jitter(aes(color = sample), size = 0.1, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        scale_y_continuous(breaks = seq(0, 100, 10)) +
        geom_hline(yintercept = 20) +
        facet_grid(cols = vars(group), scales = "free", space = "free") +
		ylab("Percent reads mapping to mito (GEX)") +
        theme_cowplot()
p6 <- ggplot(plotdata, aes(x=sample, y=TSS.enrichment)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5)) + 
        # ggrastr::rasterise(geom_jitter(aes(color = sample), size = 0.1, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        coord_cartesian(ylim = c(0, 20)) +
        scale_y_continuous(breaks = seq(0, 20, 5)) +
        geom_hline(yintercept = 1) +
        facet_grid(cols = vars(group), scales = "free", space = "free") +
		ylab("TSS enrichment score (ATAC)") +
        theme_cowplot()
p7 <- ggplot(plotdata, aes(x=sample, y=nucleosome_signal)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5)) + 
        # ggrastr::rasterise(geom_jitter(aes(color = sample), size = 0.1, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        scale_color_manual(values = bettercolors) +
        # scale_y_continuous(breaks = seq(0, 100, 10)) +
        facet_grid(cols = vars(group), scales = "free", space = "free") +
		ylab("Nucleosome Signal (ATAC)") +
        theme_cowplot()

layout <- "
ABCD
EFGH
"
pdf(file='Run10_qc.pdf', width=16, height=8)
plot(p3 + p4 + p1 + p2 + p6 + p7 + p5 + plot_spacer() + plot_layout(design = layout) & theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

####################### Filtering cells using QC cutoffs ####################

filtered <- subset(mult_integrated, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_mATAC > 200 & TSS.enrichment > 1)

DefaultAssay(filtered) <- "RNA"
filtered <- filtered %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA(npcs = 50)
filtered <- RunUMAP(filtered, reduction = "pca", dims = 1:50, reduction.name = "rnaumap", reduction.key = "rnaumap_")

DefaultAssay(filtered) <- "mATAC"
filtered <- filtered %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 20) %>% RunSVD()
filtered <- RunUMAP(filtered, reduction = 'lsi', dims = 2:50, reduction.name = "atacumap", reduction.key = "atacumap_")

filtered <- filtered %>% FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), k.nn = 50)
filtered <- RunUMAP(filtered, nn.name = "weighted.nn", reduction.name = "wnnumap", reduction.key = "wnnumap_")

saveRDS(filtered, paste0(runid, "_filtered.rds"))

#################### Call cell types ###########################

print(paste0(Sys.time(), ":: Calling cell types using SingleR ... "))
# filtered <- readRDS(paste0(runid, "_filtered.rds"))

hpca <- celldex::HumanPrimaryCellAtlasData()
# write.table(hpca %>% colData() %>% data.frame() %>% distinct(), "hpca_celltypes.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

DefaultAssay(filtered) <- "RNA"
diet_GEX <- DietSeurat(filtered, assays = "RNA")
GEX.se <- as.SingleCellExperiment(diet_GEX)

pred.GEX.hpca <- SingleR(test = GEX.se, ref = hpca, labels = hpca$label.main, BPPARAM=MulticoreParam(8))
filtered$hpca <- pred.GEX.hpca$labels

hpca_sub <- hpca[, hpca$label.main %in% c("B_cell", "DC", "Endothelial_cells", "Epithelial_cells", "Fibroblasts", "Macrophage", "Monocyte", "NK_cell", "Neutrophils", "T_cells")]
pred.GEX.hpca_sub <- SingleR(test = GEX.se, ref = hpca_sub, labels = hpca_sub$label.main, BPPARAM=MulticoreParam(8))
pred.GEX.hpca_sub_fine <- SingleR(test = GEX.se, ref = hpca_sub, labels = hpca_sub$label.fine, BPPARAM=MulticoreParam(8))

filtered$hpca_sub <- factor(pred.GEX.hpca_sub$labels, levels = names(colors.hpca_sub))
filtered$hpca_sub2 <- NULL

newmetadata <- filtered@meta.data %>% mutate(hpca_sub_fine = pred.GEX.hpca_sub_fine$labels) %>%
    mutate(hpca_sub_temp = case_when(hpca_sub == "T_cells" ~ hpca_sub_fine, TRUE ~ hpca_sub)) %>%
    mutate(hpca_sub2 = 
        case_when(
            hpca_sub_temp %in% c("T_cell:CD8+_Central_memory", "T_cell:CD8+", "T_cell:CD8+_effector_memory_RA", "T_cell:CD8+_effector_memory", "T_cell:CD8+_naive") ~ "T_cell:CD8+", 
            hpca_sub_temp %in% c("T_cell:CD4+", "T_cell:CD4+_central_memory", "T_cell:CD4+_effector_memory", "T_cell:CD4+_Naive") ~ "T_cell:CD4+", hpca_sub == "T_cells" ~ "T_cell:Other", 
            TRUE ~ hpca_sub_temp)
        ) %>%
    mutate(hpca_sub2 = factor(hpca_sub2, levels = names(colors.hpca_sub2)))
rownames(newmetadata) <- colnames(filtered)
filtered@meta.data <- newmetadata
filtered$hpca_sub_temp <- NULL

filtered$group <- case_when(filtered$sample %in% FT_samplenames ~ "FT", TRUE ~ "Tumor")
filtered$sample <- recode_factor(filtered$sample, !!!samplelabels)

filtered <- filtered %>% FindClusters(graph.name = "wsnn", resolution = 0.1, verbose = FALSE)
# pdf(file=paste0(runid, "_umaps_celltyped.pdf"), width=30, height=10)
# umaps_wrapped(filtered)
# dev.off()

filtered <- filtered %>% FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), k.nn = 5)
filtered <- filtered %>% FindClusters(graph.name = "wsnn", resolution = 0.1, verbose = FALSE)
# pdf(file=paste0(runid, "_umaps_celltyped2.pdf"), width=30, height=10)
# umaps_wrapped(filtered, leftcolors = "hpca_sub2", leftcolorpal = colors.hpca_sub2)
# dev.off()

saveRDS(pred.GEX.hpca, paste0(runid, "_hpca.rds"))
saveRDS(pred.GEX.hpca_sub, paste0(runid, "_hpca_sub.rds"))
saveRDS(pred.GEX.hpca_sub_fine, paste0(runid, "_hpca_sub_fine.rds"))
saveRDS(filtered, paste0(runid, "_celltyped.rds"))

########## Buffer region for adding information to the objects ########################################

filtered <- readRDS(paste0(runid, "_celltyped.rds"))

filtered$subgroup2 <- data.frame(sample = filtered$sample) %>% left_join(infotable) %>% pull(subgroup2)
filtered$subgroup3 <- data.frame(sample = filtered$sample) %>% left_join(infotable) %>% pull(subgroup3)

saveRDS(filtered, paste0(runid, "_celltyped.rds"))

################################## Custom Cell Typing ##########################################
filtered <- readRDS(paste0(runid, "_celltyped.rds"))

filtered$compartments <- factor(case_when(filtered$seurat_clusters %in% c(1, 4, 9, 13, 14) ~ "Stromal", TRUE ~ "Epithelial"), levels = c("Epithelial", "Stromal"))

filtered$customtype2 <- case_when(
    filtered$compartments == "Epithelial" ~ "Epithelial_cells",
    TRUE ~ filtered$hpca_sub2
    )

filtered$customtype2 <- factor(filtered$customtype2, levels = names(colors.customtype2))
filtered$customtype <- factor(case_when((filtered$customtype2 == "Epithelial_cells" & filtered$group == "Tumor") ~ "Tumor", TRUE ~ filtered$customtype2), levels = names(colors.customtype))

pdf(file=paste0(runid, "_umaps_customtype2.pdf"), width=30, height=10)
umaps_wrapped(filtered, leftcolors = "customtype2", leftcolorpal = colors.customtype2, rightcolors = "subgroup2", rightcolorpal = colors.subgroup2)
dev.off()

celltype.table1 <- data.frame(celltype = filtered$customtype2, sample = filtered$sample) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.customtype2))) %>%
    left_join(infotable %>% mutate(n = -row_number())) %>%
    mutate(sample = reorder_within(sample, n, subgroup3))
saveRDS(celltype.table1, paste0(runid, "_props_customtype2.rds"))

p1 <- ggplot(celltype.table1) + 
    geom_bar(aes(fill=celltype, x=count, y=sample), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_reordered(expand = c(0,0)) +
    scale_fill_manual(values = colors.customtype2) +
    facet_grid(rows = vars(subgroup2), scales = "free_y", space = "free_y") +
    theme_cowplot() + 
    theme(legend.position = "right", axis.title = element_blank())

pdf(file=paste0(runid, "_props_customtype2.pdf"), width=18, height=4)
plot(p1)
dev.off()

celltype.table1 <- data.frame(celltype = filtered$customtype2, sample = filtered$sample) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.customtype2))) %>%
    left_join(infotable)

celltype.table2 <- celltype.table1 %>% filter(celltype %in% c("Monocyte", "Macrophage")) %>%
    dplyr::select(celltype, sample, count) %>%
    pivot_wider(names_from = "celltype", values_from = "count") %>%
    mutate(log2MacMonoRat = log2((Macrophage+1) / (Monocyte+1))) %>%
    mutate(group = case_when(str_starts(sample, "FT") ~ "FT", TRUE ~ "Tumor")) %>%
    left_join(infotable) %>%
    mutate(subgroup3 = factor(subgroup3, levels = c("BRCAwt", "BRCAmt", "Sensitive", "Resistant", "NACT"))) %>%
	mutate(subgroup = factor(case_when(subgroup3 %in% c("BRCAwt", "BRCAmt") ~ "FT", TRUE ~ subgroup3), levels = c("FT", "Sensitive", "Resistant", "NACT")))
saveRDS(celltype.table2, paste0(runid, "_MacMonoRatio.rds"))

p1 <- ggplot(celltype.table2, aes(x = subgroup, y = log2MacMonoRat)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = subgroup3), width = 0.2, size = 3, stroke=1, shape = 17) +
    scale_color_manual(values = colors.subgroup3) +
    stat_compare_means(label.y = 5.5) +
    theme_cowplot()

pdf(file=paste0(runid, "_MacMonoRatio.pdf"), width=5, height=4)
plot(p1)
dev.off()

saveRDS(filtered, paste0(runid, "_celltyped.rds"))

################################## Epithelial Sub Typing ##########################################
Epithelial <- subset(filtered, subset = customtype2 == "Epithelial_cells")

# pdf(file=paste0(runid, "_umaps_harmony_convergence.pdf"), width=8, height=6)
# set.seed(1)
# hm_epithelial <- Epithelial %>% RunHarmony(group.by.vars = 'sample', assay.use = 'RNA', reduction = 'pca', dims.use = 1:50, reduction.save = 'harmony_pca', project.dim = FALSE, plot_convergence = TRUE, lambda=0.8)
# dev.off()

set.seed(1)
hm_epithelial <- Epithelial %>% RunHarmony(group.by.vars = 'sample', assay.use = 'RNA', reduction = 'pca', dims.use = 1:50, reduction.save = 'harmony_pca', project.dim = FALSE, plot_convergence = FALSE, lambda=0.8)
DefaultAssay(hm_epithelial) <- "RNA"
hm_epithelial %<>% RunUMAP(reduction = "harmony_pca", dims = 1:20, n.neighbors = 20, min.dist = 0.2, reduction.name = "hmrnaumap", reduction.key = "hmrnaumap_", metric = "euclidean") 
hm_epithelial %<>% FindNeighbors(reduction = "harmony_pca", dims = 1:20, annoy.metric = "euclidean", k.param = 20)
hm_epithelial %<>% FindClusters(resolution = 0.18, verbose = FALSE)

pdf(file=paste0(runid, "_umaps_harmony.pdf"), width=18, height=6)
umaps_wrapped(hm_epithelial, rnadimreduc = "hmrnaumap", rightcolors = "seurat_clusters")
dev.off()

# Save the object, then run marker identification, cell tagging (below), then decide the cluster names, then follow below for plotting. 
# saveRDS(hm_epithelial, paste0(runid, "_epithelial.rds"))

epicelltypes <- c("Ciliated", "Secretory", "Cluster1", "Cluster2", "Cluster3", "Cycling")
hm_epithelial$epicluster <- factor(case_when(
    hm_epithelial$seurat_clusters == 0 ~ "Cluster1", 
    hm_epithelial$seurat_clusters == 1 ~ "Cluster2", 
    hm_epithelial$seurat_clusters == 2 ~ "Secretory", 
    hm_epithelial$seurat_clusters == 3 ~ "Ciliated", 
    hm_epithelial$seurat_clusters == 4 ~ "Cluster3",
    hm_epithelial$seurat_clusters == 5 ~ "Cycling"), levels = epicelltypes)

pdf(file=paste0(runid, "_umaps_harmony0.pdf"), width=18, height=6)
umaps_wrapped(hm_epithelial, rightcolors = "epicluster", rightcolorpal = epicolors, , leftcolors = "subgroup2", leftcolorpal = colors.subgroup2)
dev.off()
pdf(file=paste0(runid, "_umaps_harmony1a.pdf"), width=18, height=6)
umaps_wrapped(hm_epithelial, rnadimreduc = "hmrnaumap", rightcolors = "epicluster", rightcolorpal = epicolors, , leftcolors = "subgroup2", leftcolorpal = colors.subgroup2)
dev.off()

Tumor_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & group == "Tumor"))
Tumor_Epi <- AddModuleScore(
  object = Tumor_Epi,
  features = list(Double_sig),
  assay = "remap",
  name = 'sig_remap',
  ctrl = 10,
  seed = 123
)

hm_epithelial_malig <- subset(hm_epithelial, subset = group == "Tumor")
DefaultAssay(hm_epithelial_malig) <- "RNA"
hm_epithelial_malig %<>% RunUMAP(reduction = "harmony_pca", dims = 1:20, n.neighbors = 20, min.dist = 0.2, reduction.name = "hmrnaumap", reduction.key = "hmrnaumap_", metric = "euclidean") 

pdf(file=paste0(runid, "_umaps_harmony1b.pdf"), width=18, height=6)
umaps_wrapped(hm_epithelial_malig, rnadimreduc = "hmrnaumap", rightcolors = "epicluster", rightcolorpal = epicolors, , leftcolors = "subgroup2", leftcolorpal = colors.subgroup2)
dev.off()

hm_epithelial_malig$sig_remap1 <- data.frame(cell = colnames(hm_epithelial_malig)) %>% 
	left_join(data.frame(cell = colnames(Tumor_Epi), sig_remap1 = Tumor_Epi$sig_remap1)) %>%
	pull(sig_remap1)
hm_epithelial_malig$sig_remap1[is.na(hm_epithelial_malig$sig_remap1)] <- 0

pdf(file=paste0(runid, '_harmony_sig_malig.pdf'), height = 8, width = 8)
FeaturePlot(object = hm_epithelial_malig, features = "sig_remap1", min.cutoff = "q5", max.cutoff = "q95", reduction = "hmrnaumap", pt.size = 0.5, order = TRUE)
dev.off()

DefaultAssay(hm_epithelial) <- "RNA"
hm_epithelial <- CellCycleScoring(hm_epithelial, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

Idents(hm_epithelial) <- "epicluster"
pdf(file=paste0(runid, "_cellcycle_RidgePlot1.pdf"), width=6, height=6)
RidgePlot(hm_epithelial, features = c("S.Score", "G2M.Score"), ncol = 1)
dev.off()

########### cell tagging using markers ###################
# Rscript findAllMarkers_generic.R "/working" "Run10_epithelial.rds" "seurat_clusters" "RNA" "Run10_epithelial" 1> 1 2> 2 &

cancercellmarkers <- read.table(file = "cancercellmarkers.csv", sep = ",", header = TRUE) %>% dplyr::select(Cancer.cell.1:OV.081.specific) %>% slice_head(n = 50)
temp <- read.table(file = "FT_markers.csv", sep = ",", header = TRUE) %>% 
    filter(p_val_adj < 0.05) %>% filter(cluster %in% c(1, 2)) %>% dplyr::select(cluster, gene)
FTcellmarkers <- data.frame(CE = temp %>% filter(cluster == 1) %>% slice_head(n=50) %>% pull(gene), SE = temp %>% filter(cluster == 2) %>% slice_head(n=50) %>% pull(gene))

curated_markers <- c("Secretory_ref1", "Cycling1_ref2", "Cycling2_ref2", "Ciliated1_ref1", "Ciliated2_ref1", "Ciliated_ref2")
cancerFTmarkers <- bind_cols(cancercellmarkers, FTcellmarkers) %>% dplyr::select(c("SE", "Cycling.cancer.cell.1", "Cycling.cancer.cell.2", "Ciliated.cell.1", "Ciliated.cell.2", "CE")) %>% set_colnames(curated_markers)

epi.markers <- readRDS(paste0(runid, "_epithelial_markers.rds"))
top.markers.sub <- epi.markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% slice_head(n = 50)
epicellmarkers <- data.frame(sapply(as.character(0:5), function(x){top.markers.sub %>% ungroup() %>% filter(cluster == x) %>% dplyr::select(gene) %>% unlist() %>% unname()}))

marker_overlaps <- sapply(cancerFTmarkers, function(x){
    sapply(epicellmarkers, function(y){
        intersect(x, y)
    })
})

temp <- sapply(cancerFTmarkers, function(x){
    sapply(epicellmarkers, function(y){
        length(intersect(x, y))
    })
})
temp <- data.frame(temp) %>% 
    rownames_to_column("cluster") %>%
    pivot_longer(!cluster, names_to = "curated", values_to = "overlap") %>%
    mutate(cluster = factor(cluster, levels = paste0("X", 0:5))) %>%
    # mutate(cluster = factor(cluster, levels = names(epicolors))) %>%
    mutate(curated = factor(curated, levels = curated_markers))
saveRDS(temp, paste0(runid, "_cell_tagging.rds"))

pdf(file=paste0(runid, "_cell_tagging.pdf"), width=5, height=5)
p1 <- ggplot(temp, aes(x = cluster, y = curated)) +
    geom_tile(aes(fill = overlap)) +
    coord_equal() +
    scale_fill_gradient2(low = "white", high = "red") +
    guides(fill=guide_colorbar(ticks.colour = NA)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), strip.background = element_blank())
p1
dev.off()

table1 <- data.frame(cluster = hm_epithelial$epicluster, sample = hm_epithelial$sample) %>% table() %>% data.frame() %>%
    rename(count = Freq) %>%
    mutate(sample = factor(sample, levels = samplelabels)) %>%
    group_by(sample) %>%
    mutate(freq = count / sum(count)) %>%
	left_join(filtered@meta.data %>% dplyr::select(sample, subgroup2) %>% distinct())
saveRDS(table1, paste0(runid, "_props_epithelial.rds"))

p1 <- ggplot(table1 %>% mutate(sample = factor(sample, levels = rev(samplelabels)))) + 
    geom_bar(aes(fill=cluster, y=count, x=sample), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    # scale_fill_manual(values = jdb_palette("brewer_spectra")[]) +
    scale_fill_manual(values = epicolors) +
	facet_grid(rows = vars(subgroup2), scales = "free", space = "free") +
    coord_flip() + 
    theme_cowplot() + 
    theme(legend.position = "right", axis.title = element_blank())
p2 <- ggplot(table1 %>% mutate(cluster = factor(cluster))) + 
    geom_bar(aes(fill=sample, y=count, x=cluster), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_flip() + 
    theme_cowplot() + 
    theme(legend.position = "right", strip.background = element_blank(), strip.text = element_blank(), axis.title = element_blank())

pdf(file=paste0(runid, "_props_epithelial.pdf"), width=10, height=7)
plot(p1 / p2)
dev.off()

################## SingleR scores Heatmaps #####################

hpca <- celldex::HumanPrimaryCellAtlasData()
DefaultAssay(hm_epithelial) <- "RNA"
GEX_se <- as.SingleCellExperiment(DietSeurat(hm_epithelial, assays = "RNA"))
pred_hpca <- SingleR(test = GEX_se, ref = hpca, labels = hpca$label.main, BPPARAM=MulticoreParam(8))

epicolors <- c(`Ciliated` = "#E41A1C", `Secretory` = "#449B75", `Cluster1` = "#C66764", `Cluster2` = "#999999", `Cluster3` = "#AC5782", `Cycling` = "#FFE528")
pdf(file=paste0(runid, '_epi_analysis.pdf'), width=12, height=8)
plotScoreHeatmap(pred_hpca, clusters = hm_epithelial$epicluster, order.by = "clusters", show.labels = FALSE, use_raster = TRUE, annotation_colors=list(Clusters = epicolors))
dev.off()

hpca_sub <- hpca[, hpca$label.main %in% c("Hepatocytes", "Gametocytes", "Keratinocytes", "iPS_cells", "Embryonic_stem_cells", "Astrocyte", "Neurons", "Neuroepithelial_cell", "Epithelial_cells", "Smooth_muscle_cells", "Fibroblasts", "Osteoblasts", "Tissue_stem_cells", "Chondrocytes", "Endothelial_cells", "MSC")]
pred_hpca <- SingleR(test = GEX_se, ref = hpca_sub, labels = hpca_sub$label.main, BPPARAM=MulticoreParam(8))

epicolors <- c(`Ciliated` = "#E41A1C", `Secretory` = "#449B75", `Cluster1` = "#C66764", `Cluster2` = "#999999", `Cluster3` = "#AC5782", `Cycling` = "#FFE528")
pdf(file=paste0(runid, '_epi_analysis_1.pdf'), width=12, height=4)
plotScoreHeatmap(pred_hpca, clusters = hm_epithelial$epicluster, order.by = "clusters", show.labels = FALSE, use_raster = TRUE, annotation_colors=list(Clusters = epicolors))
dev.off()

hm_PDX_Epi <- readRDS(paste0(runid, "_epithelial_PDX.rds"))
DefaultAssay(hm_PDX_Epi) <- "RNA"
GEX_se <- as.SingleCellExperiment(DietSeurat(hm_PDX_Epi, assays = "RNA"))
pred_hpca <- SingleR(test = GEX_se, ref = hpca_sub, labels = hpca_sub$label.main, BPPARAM=MulticoreParam(8))

epicolors <- c(`Ciliated` = "#E41A1C", `PDXCluster1` = "#0072B2", `PDXCluster2` = "#D55E00", `PDXCluster3` = "#A7A6DA", `PDXCluster4` = "#97D540", `Cycling` = "#FFE528")
pdf(file=paste0(runid, '_epi_analysis_1_PDX.pdf'), width=12, height=4)
plotScoreHeatmap(pred_hpca, clusters = hm_PDX_Epi$epicluster, order.by = "clusters", show.labels = FALSE, use_raster = TRUE, annotation_colors=list(Clusters = epicolors))
dev.off()

##################### DotPlots #####################

markers <- readRDS(paste0(runid, "_celltyped_markers.rds")) %>% filter(p_val_adj < 0.05)

protein_coding_genes <- read.table(file ="protein_coding_genes.tsv", header=FALSE, stringsAsFactors=FALSE) %>% pull(V1)
markers %<>% dplyr::filter(gene %in% protein_coding_genes)

topmarkers <- markers %>% group_by(cluster) %>% slice_head(n = 10)
plotmarkers <- c(unique(topmarkers$gene))
pdf(file=paste0(runid, '_dotplot_filtered_customtype2_markersTop10.pdf'), width=35, height=6)
Idents(filtered) <- "customtype2"
DotPlot(filtered, assay = "RNA", features = plotmarkers, cluster.idents = FALSE) +
    scale_color_gradientn(colors = jdb_palette("solar_rojos"), limits = c(-2.5, 2.5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic", size=20), axis.title = element_blank(), axis.text.y = element_text(size=20))
dev.off()

epi_markers <- readRDS(paste0(runid, "_epithelial_markers.rds")) %>% filter(p_val_adj < 0.05)
protein_coding_genes <- read.table(file ="protein_coding_genes.tsv", header=FALSE, stringsAsFactors=FALSE) %>% pull(V1)
epi_markers %<>% dplyr::filter(gene %in% protein_coding_genes)
epicelltypes <- c("Ciliated", "Secretory", "Cluster1", "Cluster2", "Cluster3", "Cycling")
epi_markers$epicluster <- factor(case_when(
    epi_markers$cluster == 0 ~ "Cluster1", 
    epi_markers$cluster == 1 ~ "Cluster2", 
    epi_markers$cluster == 2 ~ "Secretory", 
    epi_markers$cluster == 3 ~ "Ciliated", 
    epi_markers$cluster == 4 ~ "Cluster3",
    epi_markers$cluster == 5 ~ "Cycling"), levels = epicelltypes)

topmarkers <- epi_markers %>% group_by(epicluster) %>% slice_head(n = 10)
plotmarkers <- c(unique(topmarkers$gene))
pdf(file=paste0(runid, '_dotplot_epithelial_markersTop10.pdf'), width=16, height=4)
Idents(hm_epithelial) <- "epicluster"
DotPlot(hm_epithelial, assay = "RNA", features = plotmarkers, cluster.idents = FALSE) +
    scale_color_gradientn(colors = jdb_palette("solar_rojos"), limits = c(-2.5, 2.5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic", size=20), axis.title = element_blank(), axis.text.y = element_text(size=20))
dev.off()

######### Save subtype objects ##############
dietfiltered <- DietSeurat(filtered, assays = c("RNA", "mATAC"))
Tumor <- subset(dietfiltered, subset = group == "Tumor")
FT <- subset(dietfiltered, subset = group == "FT")
Naive <- subset(dietfiltered, subset = subgroup2 == "Naive")
NACT <- subset(dietfiltered, subset = subgroup2 == "NACT")
Tumor <- recalculate(Tumor, rnaassay = "RNA", atacassay = "mATAC")
FT <- recalculate(FT, rnaassay = "RNA", atacassay = "mATAC")
Naive <- recalculate(Naive, rnaassay = "RNA", atacassay = "mATAC")
NACT <- recalculate(NACT, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(Tumor, paste0(runid, "_Tumor.rds"))
saveRDS(FT, paste0(runid, "_FT.rds"))
saveRDS(Naive, paste0(runid, "_Naive.rds"))
saveRDS(NACT, paste0(runid, "_NACT.rds"))

Epithelial <- subset(filtered, subset = customtype2 == "Epithelial_cells")
Epithelial <- recalculate(Epithelial, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(Epithelial, paste0(runid, "_epithelial_raw.rds"))

Epithelial <- subset(dietfiltered, subset = customtype2 == "Epithelial_cells")
Epithelial <- recalculate(Epithelial, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(Epithelial, paste0(runid, "_Epith.rds"))

subobj <- subset(dietfiltered, subset = customtype2 == "T_cell:CD8+")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_T_cell_CD8.rds"))
subobj <- subset(dietfiltered, subset = customtype2 == "NK_cell")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_NK_cell.rds"))
subobj <- subset(dietfiltered, subset = customtype2 == "Monocyte")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Monocyte.rds"))
subobj <- subset(dietfiltered, subset = customtype2 == "Macrophage")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Macrophage.rds"))
subobj <- subset(dietfiltered, subset = customtype2 == "Endothelial_cells")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Endothelial_cells.rds"))
subobj <- subset(dietfiltered, subset = customtype2 == "Fibroblasts")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Fibroblasts.rds"))

Epithelial <- subset(dietfiltered, subset = customtype2 == "Epithelial_cells")
Epithelial_Tumor <- subset(Epithelial, subset = group == "Tumor")
Epithelial_Tumor <- recalculate(Epithelial_Tumor, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(Epithelial_Tumor, paste0(runid, "_Epith_Tumor.rds"))

FT <- readRDS(paste0(runid, "_FT.rds"))
Naive <- readRDS(paste0(runid, "_Naive.rds"))
NACT <- readRDS(paste0(runid, "_NACT.rds"))

eFT <- subset(FT, subset = customtype2 == "Epithelial_cells")
eFT <- recalculate(eFT, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(eFT, paste0(runid, "_Epith_FT.rds"))
eNaive <- subset(Naive, subset = customtype2 == "Epithelial_cells")
eNaive <- recalculate(eNaive, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(eNaive, paste0(runid, "_Epith_Naive.rds"))
eNACT <- subset(NACT, subset = customtype2 == "Epithelial_cells")
eNACT <- recalculate(eNACT, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(eNACT, paste0(runid, "_Epith_NACT.rds"))

dietobj <- DietSeurat(hm_epithelial, assays = c("RNA", "mATAC"))
subobj <- subset(dietobj, subset = epicluster == "Ciliated")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Ciliated.rds"))
subobj <- subset(dietobj, subset = epicluster == "Secretory")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Secretory.rds"))
subobj <- subset(dietobj, subset = epicluster == "Cluster1")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Cluster1.rds"))
subobj <- subset(dietobj, subset = epicluster == "Cluster2")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Cluster2.rds"))
subobj <- subset(dietobj, subset = epicluster == "Cluster3")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Cluster3.rds"))
subobj <- subset(dietobj, subset = epicluster == "Cycling")
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Cycling.rds"))

subobj <- subset(dietobj, subset = (epicluster == "Ciliated" & group == "FT"))
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Ciliated_FT.rds"))
subobj <- subset(dietobj, subset = (epicluster == "Secretory" & group == "FT"))
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Secretory_FT.rds"))
subobj <- subset(dietobj, subset = (epicluster == "Ciliated" & group == "Tumor"))
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Ciliated_Tumor.rds"))
subobj <- subset(dietobj, subset = (epicluster == "Secretory" & group == "Tumor"))
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Secretory_Tumor.rds"))
subobj <- subset(dietobj, subset = (epicluster == "Cluster1" & group == "Tumor"))
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Cluster1_Tumor.rds"))
subobj <- subset(dietobj, subset = (epicluster == "Cluster2" & group == "Tumor"))
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Cluster2_Tumor.rds"))
subobj <- subset(dietobj, subset = (epicluster == "Cluster3" & group == "Tumor"))
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Cluster3_Tumor.rds"))
subobj <- subset(dietobj, subset = (epicluster == "Cycling" & group == "Tumor"))
subobj <- recalculate(subobj, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(subobj, paste0(runid, "_Cycling_Tumor.rds"))

PDX <- readRDS("RunPDX1_filtered.rds")
M_cells <- read.table("MouseCells.tsv") %>% pull(V1)
H_cells <- setdiff(colnames(PDX), M_cells)
PDX$CellName <- colnames(PDX)
PDX <- subset(PDX, subset = CellName %in% H_cells)
PDX <- recalculate(PDX, rnaassay = "RNA", atacassay = "mATAC")
saveRDS(PDX, paste0(runid, "_Epith_PDX.rds"))

######### Independent runs of time consuming steps (run via slurm - temporary slurm scripts) ###########

# # Markers
# Rscript findAllMarkers_generic.R "/working" "Run10_celltyped.rds" "customtype2" "RNA" "Run10_celltyped" 1> 1 2> 2 &

# # ChromVAR (Saved in same object)
# Rscript ChromVAR_generic.R /working "Run10_celltyped.rds" 1> 1 2> 2 &
# Rscript ChromVAR_generic.R /working "Run10_epithelial.rds" 1> 11 2> 22 &

# # DEGs
# Rscript degs_generic.R /working Run10_celltyped.rds Tumor FT customtype2 group Run10 1> 1 2> 2 &
# Rscript degs_generic.R /working Run10_celltyped.rds BRCAmt BRCAwt customtype2 subgroup3 Run10 1> 1 2> 2 &
# Rscript degs_generic.R /working Run10_celltyped.rds NACT Naive customtype2 subgroup2 Run10 1> 1 2> 2 &
# Rscript degs_generic.R /working Run10_celltyped.rds Naive FT customtype2 subgroup2 Run10 1> 1 2> 2 &
# Rscript degs_generic.R /working Run10_celltyped.rds NACT FT customtype2 subgroup2 Run10 1> 1 2> 2 &

# Rscript degs_generic.R /working Run10_celltyped.rds Resistant Sensitive customtype2 subgroup3 Run10 1> 1 2> 2 &
# Rscript degs_generic.R /working Run10_celltyped.rds NACT Sensitive customtype2 subgroup3 Run10 1> 1 2> 2 &
# Rscript degs_generic.R /working Run10_celltyped.rds NACT Resistant customtype2 subgroup3 Run10 1> 1 2> 2 &

# # Cicero
# Rscript Cicero_generic.R /working Run10_T_cell_CD8.rds seurat_clusters WNNUMAP Run10_T_cell_CD8 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_NK_cell.rds seurat_clusters WNNUMAP Run10_NK_cell 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Monocyte.rds seurat_clusters WNNUMAP Run10_Monocyte 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Macrophage.rds seurat_clusters WNNUMAP Run10_Macrophage 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Endothelial_cells.rds seurat_clusters WNNUMAP Run10_Endothelial_cells 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Fibroblasts.rds seurat_clusters WNNUMAP Run10_Fibroblasts 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Epith_Tumor.rds seurat_clusters WNNUMAP Run10_Epith_Tumor 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Epith_FT.rds seurat_clusters WNNUMAP Run10_Epith_FT 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Epith.rds seurat_clusters WNNUMAP Run10_Epith 1> 1 2> 2 &

# Rscript Cicero_generic.R /working Run10_Ciliated.rds seurat_clusters WNNUMAP Run10_Ciliated 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Secretory.rds seurat_clusters WNNUMAP Run10_Secretory 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Cluster1.rds seurat_clusters WNNUMAP Run10_Cluster1 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Cluster2.rds seurat_clusters WNNUMAP Run10_Cluster2 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Cluster3.rds seurat_clusters WNNUMAP Run10_Cluster3 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Cycling.rds seurat_clusters WNNUMAP Run10_Cycling 1> 1 2> 2 &

# Rscript Cicero_generic.R /working Run10_Ciliated_FT.rds seurat_clusters WNNUMAP Run10_Ciliated_FT 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Secretory_FT.rds seurat_clusters WNNUMAP Run10_Secretory_FT 1> 1 2> 2 &

# Rscript Cicero_generic.R /working Run10_Ciliated_Tumor.rds seurat_clusters WNNUMAP Run10_Ciliated_Tumor 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Secretory_Tumor.rds seurat_clusters WNNUMAP Run10_Secretory_Tumor 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Cluster1_Tumor.rds seurat_clusters WNNUMAP Run10_Cluster1_Tumor 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Cluster2_Tumor.rds seurat_clusters WNNUMAP Run10_Cluster2_Tumor 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Cluster3_Tumor.rds seurat_clusters WNNUMAP Run10_Cluster3_Tumor 1> 1 2> 2 &
# Rscript Cicero_generic.R /working Run10_Cycling_Tumor.rds seurat_clusters WNNUMAP Run10_Cycling_Tumor 1> 1 2> 2 &

# # FigR
# Rscript FigR_generic.R /working Run10_epithelial_raw.rds Run10_Epith mATAC hg38 1
# Rscript FigR_generic.R /working Run10_epithelial_raw.rds Run10_Epith mATAC hg38 2
# Rscript FigR_generic.R /working Run10_epithelial_raw.rds Run10_Epith mATAC hg38 3

# Rscript FigR_generic.R /working Run10_Epith_FT.rds Run10_Epith_FT mATAC hg38 1
# Rscript FigR_generic.R /working Run10_Epith_FT.rds Run10_Epith_FT mATAC hg38 2
# Rscript FigR_generic.R /working Run10_Epith_FT.rds Run10_Epith_FT mATAC hg38 3

# Rscript FigR_generic.R /working Run10_Epith_NACT.rds Run10_Epith_NACT mATAC hg38 1
# Rscript FigR_generic.R /working Run10_Epith_NACT.rds Run10_Epith_NACT mATAC hg38 2
# Rscript FigR_generic.R /working Run10_Epith_NACT.rds Run10_Epith_NACT mATAC hg38 3

# Rscript FigR_generic.R /working Run10_Epith_Naive.rds Run10_Epith_Naive mATAC hg38 1
# Rscript FigR_generic.R /working Run10_Epith_Naive.rds Run10_Epith_Naive mATAC hg38 2
# Rscript FigR_generic.R /working Run10_Epith_Naive.rds Run10_Epith_Naive mATAC hg38 3

# Rscript FigR_generic.R /working Run10_Epith_PDX.rds Run10_Epith_PDX mATAC hg38 1
# Rscript FigR_generic.R /working Run10_Epith_PDX.rds Run10_Epith_PDX mATAC hg38 2
# Rscript FigR_generic.R /working Run10_Epith_PDX.rds Run10_Epith_PDX mATAC hg38 3

# Rscript FigR_generic.R /working Run10_Epith_Tumor.rds Run10_Epith_Tumor mATAC hg38 1
# Rscript FigR_generic.R /working Run10_Epith_Tumor.rds Run10_Epith_Tumor mATAC hg38 2
# Rscript FigR_generic.R /working Run10_Epith_Tumor.rds Run10_Epith_Tumor mATAC hg38 3

############################ Differential expression volcano and differential enrichment of DBFs #############################

plot_differential_expression <- function(prefix, leftgroup, rightgroup, celltypes, celltype_colors, unitlabels){
    degsr <- read.table(file = paste0(prefix, "_degsr_", rightgroup, "_", leftgroup, ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))
    
    temp1 <- degsr %>% 
        filter(p_val_adj < 0.05) %>% 
        mutate(celltype = factor(celltype, levels = names(celltype_colors))) %>% 
        mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
        mutate(neg_log10_adj_pval = -log10(p_val_adj)) %>%
        filter(celltype %in% celltypes)
    
    forlabels <- temp1 %>% mutate(gene = case_when((gene %in% unitlabels) ~ gene, TRUE ~ ""))
    
    p1 <- ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
        ggrastr::rasterise(geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
        geom_hline(yintercept=-log10(0.05), linetype="dashed") +
        geom_label_repel(data = forlabels %>% filter(avg_log2FC > 0), aes(label = gene, segment.color=celltype), size = 3, max.overlaps = Inf, hjust = 0, nudge_x = 8, direction = "y", segment.size = 0.5, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 3, fontface = "italic") +
        geom_label_repel(data = forlabels %>% filter(avg_log2FC < 0), aes(label = gene, segment.color=celltype), size = 3, max.overlaps = Inf, hjust = 1, nudge_x = -8, direction = "y", segment.size = 0.5, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 3, fontface = "italic") +
        # coord_cartesian(xlim = c(-4, 4)) +
        scale_color_manual(values = celltype_colors, aesthetics = c("color", "segment.color")) +
        labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() +
        theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title=element_blank())
    p2 <- table((temp1 %>% filter(avg_log2FC < 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(celltype_colors))), direction = "DOWN") %>%
        bind_rows(table((temp1 %>% filter(avg_log2FC > 0))$celltype) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(celltype_colors))), direction = "UP")) %>%
        mutate(direction = factor(direction, levels = c("UP", "DOWN"))) %>% 
        ggplot(aes(fill=Var1, y=Freq, x=direction)) +
        geom_bar(stat="identity", width = 0.6) + 
        scale_fill_manual(values = celltype_colors) + 
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_flip() +
        theme_cowplot() +
        theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), plot.margin = margin(0.5, 0.1, 0.5, 0.1, "cm"))

    pdf(file = paste0(prefix, "_degsr_", rightgroup, "_", leftgroup, ".pdf"), width = 9, height = 9)
    plot(p1 + p2 + plot_layout(heights = c(15, 1)) + plot_annotation(theme = theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))))
    dev.off()
}

plot_differential_expression("Run10", "FT", "Naive", names(colors.customtype2), colors.customtype2, "")

differential_enrichment <- function(obj, assay, celltype, group, leftgroup, rightgroup, prefix){
    obj$celltype.group <- paste(obj@meta.data %>% pull(!!celltype), obj@meta.data %>% pull(!!group), sep = ".")
    Idents(obj) <- "celltype.group"
    DefaultAssay(obj) <- assay

    types <- unique(obj@meta.data %>% pull(!!celltype) %>% as.character())

    deus <- lapply(
        types, 
        FUN = function(type){
            curridents <- unique(Idents(obj))
            if (!(paste0(type, ".", leftgroup) %in% curridents) | !(paste0(type, ".", rightgroup) %in% curridents)) return(data.frame())
            if (length(WhichCells(obj, idents = paste0(type, ".", leftgroup))) < 3 | length(WhichCells(obj, idents = paste0(type, ".", rightgroup))) < 3) return(data.frame())
            temp1 <- FindMarkers(obj, ident.1 = paste0(type, ".", rightgroup), ident.2 = paste0(type, ".", leftgroup), mean.fxn = rowMeans, fc.name = "avg_diff") %>% rownames_to_column("unit")
            return(temp1)
        }
    )
    names(deus) <- types
    deus <- bind_rows(deus, .id = "celltype") 
    write.table(deus, file = paste0(prefix, "_differential_", assay, "_", rightgroup, "_", leftgroup, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    obj <- DietSeurat(obj, assays = assay)
    obj[[assay]] <- as(obj[[assay]], "Assay5")

    samples <- obj@meta.data %>% pull(sample) %>% as.character() %>% unique()
    oneouts <- lapply(samples, function(sam){subset(obj, subset = sample != sam)})

    deus_robust <- lapply(
        oneouts, 
        FUN = function(oneout){
            currsamples <- oneout@meta.data %>% pull(sample) %>% as.character() %>% unique()
            print(setdiff(samples, currsamples))
            curridents <- unique(Idents(oneout))
            ctsp_deus_list <- lapply(
                types,
                FUN = function(type){
                    print(paste0("    ", type))
                    if (!(paste0(type, ".", leftgroup) %in% curridents) | !(paste0(type, ".", rightgroup) %in% curridents)) return(c())
                    if (length(WhichCells(oneout, idents = paste0(type, ".", leftgroup))) < 3 | length(WhichCells(oneout, idents = paste0(type, ".", rightgroup))) < 3) return(c())
                    temp1 <- FindMarkers(oneout, ident.1 = paste0(type, ".", rightgroup), ident.2 = paste0(type, ".", leftgroup), mean.fxn = rowMeans, fc.name = "avg_diff") %>% rownames_to_column("unit")
                    return((temp1 %>% filter(p_val_adj < 0.05))$unit)
                }
            )
            names(ctsp_deus_list) <- types
            return(ctsp_deus_list)
        }
    )
    names(deus_robust) <- samples
    deus_robust <- lapply(deus_robust %>% purrr::transpose(), function(x){Reduce(intersect, x)})
    deusr <- bind_rows(lapply(deus_robust, function(unit){data.frame(unit)}), .id = "celltype") %>% left_join(deus)
    write.table(deusr, file = paste0(prefix, "_differentialr_", assay, "_", rightgroup, "_", leftgroup, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

differential_enrichment(filtered, "chromvar", "customtype2", "group", "FT", "Tumor", runid)
differential_enrichment(filtered, "remap", "customtype2", "group", "FT", "Tumor", runid)

FT <- subset(filtered, subset = group == "FT")
Tumor <- subset(filtered, subset = group == "Tumor")

differential_enrichment(FT, "chromvar", "customtype2", "subgroup3", "BRCAwt", "BRCAmt", runid)
differential_enrichment(FT, "remap", "customtype2", "subgroup3", "BRCAwt", "BRCAmt", runid)
differential_enrichment(Tumor, "chromvar", "customtype2", "subgroup2", "Naive", "NACT", runid)
differential_enrichment(Tumor, "remap", "customtype2", "subgroup2", "Naive", "NACT", runid)

Naive <- subset(Tumor, subset = subgroup2 == "Naive")

differential_enrichment(Naive, "chromvar", "customtype2", "subgroup3", "Sensitive", "Resistant", runid)
differential_enrichment(Naive, "remap", "customtype2", "subgroup3", "Sensitive", "Resistant", runid)

temp <- subset(filtered, subset = subgroup2 != "NACT")

differential_enrichment(temp, "chromvar", "customtype2", "subgroup2", "FT", "Naive", runid)
differential_enrichment(temp, "remap", "customtype2", "subgroup2", "FT", "Naive", runid)

####################### Figure 4 ################################

plot_volcano1 <- function(infile, outfile, unitlabels, nudge_val = 6){
    tbl <- read.table(file = infile, sep = "\t", skip = 1) %>%
		set_colnames(c("celltype", "unit", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj"))
    
    temp1 <- tbl %>% 
        filter(p_val_adj < 0.05) %>% 
        mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
        mutate(neg_log10_adj_pval = -log10(p_val_adj)) %>%
        filter(celltype %in% "Epithelial_cells")
    
    forlabels <- temp1 %>% mutate(unit = case_when((unit %in% unitlabels) ~ unit, TRUE ~ ""))
	
    p1 <- ggplot(temp1, aes(x = effect_size, y = neg_log10_adj_pval)) +
        ggrastr::rasterise(geom_point(data = temp1 %>% filter(!(unit %in% unitlabels)), color = "#e78ac3", size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
		ggrastr::rasterise(geom_point(data = temp1 %>% filter(unit %in% unitlabels), color = "black", fill = colors.customtype2["Epithelial_cells"], size = 1.5, stroke=0.3, shape = 21), dpi = 400) +
        geom_label_repel(data = forlabels %>% filter(effect_size > 0), aes(label = unit), size = 4, max.overlaps = Inf, hjust = 0, nudge_x = nudge_val, direction = "y", segment.size = 0.1, point.padding = 0.1, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 3, fontface = "italic") +
        geom_label_repel(data = forlabels %>% filter(effect_size < 0), aes(label = unit), size = 4, max.overlaps = Inf, hjust = 1, nudge_x = -nudge_val, direction = "y", segment.size = 0.1, point.padding = 0.1, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 3, fontface = "italic") +
        # coord_cartesian(xlim = c(-4, 4)) +
        # scale_color_manual(values = color_values, aesthetics = c("color")) +
        theme_cowplot() +
        theme(axis.title = element_blank())
    
    pdf(file = outfile, width = 5, height = 6)
    plot(p1)
    dev.off()
}

plot_volcano2 <- function(infile, outfile, unitlabels, nudge_val = 6){
    tbl <- read.table(file = infile, sep = "\t", skip = 1) %>%
		set_colnames(c("celltype", "unit", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj"))
    
    temp1 <- tbl %>% 
        filter(p_val_adj < 0.05) %>% 
        mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
        mutate(neg_log10_adj_pval = -log10(p_val_adj)) %>%
        filter(celltype %in% "Epithelial_cells")
    
    # forlabels <- temp1 %>% mutate(unit = case_when((unit %in% unitlabels) ~ unit, TRUE ~ ""))
	forlabels <- temp1 %>% left_join(unitlabels) %>% drop_na()
	color_values = unitlabels %>% dplyr::select(group, color) %>% distinct() %>% pull(color)
	names(color_values) = unitlabels %>% dplyr::select(group, color) %>% distinct() %>% pull(group)
    
    p1 <- ggplot(temp1, aes(x = effect_size, y = neg_log10_adj_pval)) +
        ggrastr::rasterise(geom_point(data = temp1 %>% filter(!(unit %in% unitlabels$unit)), color = "#e78ac3", size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
		ggrastr::rasterise(geom_point(data = temp1 %>% filter(unit %in% unitlabels$unit), color = "black", fill = colors.customtype2["Epithelial_cells"], size = 1.5, stroke=0.3, shape = 21), dpi = 400) +
        geom_label_repel(data = forlabels %>% filter(effect_size > 0), aes(label = unit, color = group), size = 4, max.overlaps = Inf, hjust = 0, nudge_x = nudge_val, direction = "y", segment.size = 0.1, point.padding = 0.1, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 3, fontface = "italic") +
        geom_label_repel(data = forlabels %>% filter(effect_size < 0), aes(label = unit, color = group), size = 4, max.overlaps = Inf, hjust = 1, nudge_x = -nudge_val, direction = "y", segment.size = 0.1, point.padding = 0.1, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 3, fontface = "italic") +
        # coord_cartesian(xlim = c(-4, 4)) +
        scale_color_manual(values = color_values, aesthetics = c("color")) +
        theme_cowplot() +
        theme(axis.title = element_blank())
    
    pdf(file = outfile, width = 5, height = 6)
    plot(p1)
    dev.off()
}

NN_labels = data.frame(unit = c(base::intersect(Matondo_Chemoresponse, NACT), base::intersect(stress_associated, NACT)), 
	color = c(rep("#377EB8", length(base::intersect(Matondo_Chemoresponse, NACT))), rep("#e01c1a", length(base::intersect(stress_associated, NACT)))), 
	group = c(rep("Matondo_Chemoresponse", length(base::intersect(Matondo_Chemoresponse, NACT))), rep("stress_associated", length(base::intersect(stress_associated, NACT)))))

plot_volcano2("Run10_degsr_NACT_Naive.tsv", "Run10_degsr_NACT_Naive.pdf", NN_labels, nudge_val = 0)

SR_labels = data.frame(unit = c(base::intersect(MultiOv, Resistant), base::intersect(stress_associated, Resistant)), 
	color = c(rep("#008018", length(base::intersect(MultiOv, Resistant))), rep("#e01c1a", length(base::intersect(stress_associated, Resistant)))), 
	group = c(rep("MultiOv", length(base::intersect(MultiOv, Resistant))), rep("stress_associated", length(base::intersect(stress_associated, Resistant)))))

plot_volcano2("Run10_degsr_Resistant_Sensitive.tsv", "Run10_degsr_Resistant_Sensitive.pdf", SR_labels, nudge_val = 0)

plot_volcano1("Run10_daps_Naive_FT.tsv", "Run10_daps_Naive_FT.pdf", c(), nudge_val = 0)
plot_volcano1("Run10_daps_NACT_Naive.tsv", "Run10_daps_NACT_Naive.pdf", c(), nudge_val = 0)
plot_volcano1("Run10_daps_Resistant_Sensitive.tsv", "Run10_daps_Resistant_Sensitive.pdf", c(), nudge_val = 0)

plot_volcano1("Run10_differential_remap_NACT_Naive.tsv", "Run10_differential_remap_NACT_Naive.pdf", OvSp, nudge_val = 2)
plot_volcano1("Run10_differential_remap_Resistant_Sensitive.tsv", "Run10_differential_remap_Resistant_Sensitive.pdf", OvSp, nudge_val = 2)

library(eulerr)

set1 <- read.table(file = "Run10_degsr_NACT_Naive.tsv", sep = "\t", skip = 1) %>% 
	set_colnames(c("celltype", "unit", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & effect_size > 0 & celltype == "Epithelial_cells") %>% pull(unit)
set2 <- read.table(file = "Run10_degsr_Resistant_Sensitive.tsv", sep = "\t", skip = 1) %>% 
	set_colnames(c("celltype", "unit", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & effect_size > 0 & celltype == "Epithelial_cells") %>% pull(unit)

length(set1)
length(set2)
length(base::intersect(set1, set2))

pdf(file = "Venn1.pdf", width = 3, height = 3)
fit <- euler(c(A = length(setdiff(set1, set2)), B = length(setdiff(set2, set1)), "A&B" = length(base::intersect(set1, set2))))
plot(fit, quantities = TRUE, fills = c("#326c7b", "#952320"))
dev.off()

set1 <- read.table(file = "Run10_daps_NACT_Naive.tsv", sep = "\t", skip = 1) %>% 
	set_colnames(c("celltype", "unit", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & effect_size > 0 & celltype == "Epithelial_cells") %>% pull(unit)
set2 <- read.table(file = "Run10_daps_Resistant_Sensitive.tsv", sep = "\t", skip = 1) %>% 
	set_colnames(c("celltype", "unit", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & effect_size > 0 & celltype == "Epithelial_cells") %>% pull(unit)

length(set1)
length(set2)
length(base::intersect(set1, set2))

pdf(file = "Venn2.pdf", width = 3, height = 3)
fit <- euler(c(A = length(setdiff(set1, set2)), B = length(setdiff(set2, set1)), "A&B" = length(base::intersect(set1, set2))))
plot(fit, quantities = TRUE, fills = c("#326c7b", "#952320"))
dev.off()

set1 <- read.table(file = "Run10_differential_remap_NACT_Naive.tsv", sep = "\t", skip = 1) %>% 
	set_colnames(c("celltype", "unit", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & effect_size > 0 & celltype == "Epithelial_cells") %>% pull(unit)
set2 <- read.table(file = "Run10_differential_remap_Resistant_Sensitive.tsv", sep = "\t", skip = 1) %>% 
	set_colnames(c("celltype", "unit", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & effect_size > 0 & celltype == "Epithelial_cells") %>% pull(unit)

length(set1)
length(set2)
length(base::intersect(set1, set2))

pdf(file = "Venn3.pdf", width = 3, height = 3)
fit <- euler(c(A = length(setdiff(set1, set2)), B = length(setdiff(set2, set1)), "A&B" = length(base::intersect(set1, set2))))
plot(fit, quantities = TRUE, fills = c("#326c7b", "#952320"))
dev.off()

################################# Define Persister Signature #########################################

Tumor_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & group == "Tumor"))

# # DAPS
# Rscript daps_generic.R /Run10_v2 Run10_celltyped_v2.rds Naive FT customtype2 subgroup2 Run10 1> 1 2> 2 &
# Rscript daps_generic.R /Run10_v2 Run10_celltyped_v2.rds NACT Naive customtype2 subgroup2 Run10 1> 1 2> 2 &
# Rscript daps_generic.R /Run10_v2 Run10_celltyped_v2.rds Resistant Sensitive customtype2 subgroup3 Run10 1> 11 2> 22 &

NaivNACT_P <- read.table("Run10_daps_NACT_Naive.tsv", sep = "\t", header = TRUE) %>% drop_na() %>% dplyr::filter(celltype =="Epithelial_cells" & p_val_adj < 0.05 & avg_log2FC > 0) %>% pull(gene)
SensResi_P <- read.table("Run10_daps_Resistant_Sensitive.tsv", sep = "\t", header = TRUE) %>% drop_na() %>% dplyr::filter(celltype =="Epithelial_cells" & p_val_adj < 0.05 & avg_log2FC > 0) %>% pull(gene)
length(setdiff(NaivNACT_P, SensResi_P))
length(base::intersect(NaivNACT_P, SensResi_P))
length(setdiff(SensResi_P, NaivNACT_P))
patient_EPS_peaks <- base::intersect(NaivNACT_P, SensResi_P)
write.table(data.frame(peak = patient_EPS_peaks) %>% separate(peak, c("chr", "start", "end"), sep = "-") %>% arrange(chr, as.numeric(start)), "patient_EPS_peaks.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

NaivNACT_P <- read.table("Run10_daps_NACT_Naive.tsv", sep = "\t", header = TRUE) %>% drop_na() %>% dplyr::filter(celltype =="Epithelial_cells" & p_val_adj < 0.05)
SensResi_P <- read.table("Run10_daps_Resistant_Sensitive.tsv", sep = "\t", header = TRUE) %>% drop_na() %>% dplyr::filter(celltype =="Epithelial_cells" & p_val_adj < 0.05)
temp1 <- NaivNACT_P %>% 
	filter(p_val_adj < 0.05) %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))
temp2 <- SensResi_P %>% 
	filter(p_val_adj < 0.05) %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))

pdf(file="Run10_DAPs.pdf", width = 12, height = 6)
p1 <- ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(size = 0.5, stroke=0.1, shape = 16, color = colors.customtype2["Epithelial_cells"]), dpi = 400) +
	theme_cowplot() + ggtitle("Naive vs. NACT")
p2 <- ggplot(temp2, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(size = 0.5, stroke=0.1, shape = 16, color = colors.customtype2["Epithelial_cells"]), dpi = 400) +
	theme_cowplot() + ggtitle("Sensitive vs. Resistant")
p1 | p2
dev.off()

library(ReMapEnrich)
remapCatalog <- bedToGranges("myremap_v2.bed")
set.seed(1)
en <- enrichment(patient_EPS_peaks %>% StringToGRanges(), remapCatalog, chromSizes = loadChromSizes("hg38"), byChrom = FALSE)
temp <- data.frame(en) %>% drop_na() %>% arrange(-`q.significance`) %>% filter(category %in% rownames(filtered$RNA))
plotdata <- temp %>% mutate(category = factor(category, levels = temp$category)) %>% mutate(rn = row_number())
patient_EPS_factors <- plotdata %>% pull(category) %>% head(n = 120) %>% as.character() %>% sort()
top120 <- plotdata %>% pull(category) %>% head(n = 120) %>% as.character() %>% sort()

temp <- data.frame(en) %>% drop_na() %>% arrange(`q.significance`)
plotdata2 <- temp %>% mutate(category = factor(category, levels = temp$category)) %>% mutate(issig = category %in% patient_EPS_factors)
pdf(file=paste0(runid, "_EPSv2_step1.pdf"), width=30, height=12)
p1 <- ggplot(plotdata2 %>% filter(category %in% patient_EPS_factors), aes(x = category, y = `q.significance`)) +
	geom_col() + 
	theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- ggplot(data = plotdata2 %>% arrange(issig), mapping = aes(x = category, y = `q.significance`, color = issig)) +
    ggrastr::rasterise(geom_point(size=2, stroke=0.1, shape=16), dpi=400) +
	scale_color_manual(values = c("grey70", "red")) +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p1 / p2
dev.off()

NeoExprSignature <- Tumor_Epi$RNA$counts[patient_EPS_factors,] %>% t() %>% data.frame() %>% 
	mutate(group = Tumor_Epi$subgroup2) %>%
	filter(group == "NACT") %>%
	dplyr::select(-group) %>%
	pivot_longer(everything(), names_to = "gene", values_to = "expr") %>%
	mutate(expressed = case_when(expr > 0 ~ 1, TRUE ~ 0)) %>%
	group_by(gene) %>%
	summarize(percexpr = sum(expressed)*100/n())

patient_EPS_factors <- NeoExprSignature %>% arrange(-percexpr) %>% head(n = 100) %>% pull(gene) %>% sort()
write.table(data.frame(factor = patient_EPS_factors), "patient_EPS_factors.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

temp <- data.frame(en) %>% drop_na() %>% arrange(`q.significance`) %>% filter(category %in% rownames(filtered$RNA))
plotdata2 <- temp %>% mutate(category = factor(category, levels = temp$category)) %>% mutate(issig = category %in% patient_EPS_factors)
# labs = plotdata2 %>% filter(category %in% patient_EPS_factors) %>% pull(category) %>% tail(n = 15)
labs = c("CREBBP", "EP300", "BRD2", "BRD4", "HDAC1", "HDAC2", "SMARCB1", "ARID1A", "SMARCA4", "SMARCA2", "SMARCC1", "SMARCE1", "MAML1", "RBPJ", "SOX2")
pdf(file=paste0(runid, "_EPSv2_step2.pdf"), width=6, height=6)
ggplot(data = plotdata2 %>% arrange(issig), mapping = aes(x = category, y = `q.significance`, color = issig)) +
    ggrastr::rasterise(geom_point(size=2, stroke=0.1, shape=16), dpi=400) +
	ggrastr::rasterise(geom_point(data = plotdata2 %>% filter(category %in% labs), color = "black", fill = "red", size = 2, stroke=0.3, shape = 21), dpi=400) +
	geom_label_repel(data = plotdata2 %>% filter(category %in% labs), aes(label = category), size = 6, max.overlaps = Inf, force = 3, segment.size = 0.1) +
	scale_color_manual(values = c("grey70", "red")) +
	scale_x_discrete(expand=c(0.1,0.1)) +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
dev.off()

temp <- data.frame(en) %>% drop_na() %>% arrange(`q.significance`) %>% filter(category %in% rownames(filtered$RNA))
plotdata2 <- temp %>% mutate(category = factor(category, levels = temp$category)) %>% mutate(issig = category %in% patient_EPS_factors) %>%
	mutate(issig = case_when(!(category %in% top120) ~ NA, TRUE ~ issig))
pdf(file=paste0(runid, "_EPSv2_step2_v2.pdf"), width=6, height=6)
ggplot(data = plotdata2 %>% arrange(issig), mapping = aes(x = category, y = `q.significance`, color = issig)) +
    ggrastr::rasterise(geom_point(size=2, stroke=0.1, shape=16), dpi=400) +
	scale_color_manual(values = c("grey70", "red"), na.value = "black") +
	scale_x_discrete(expand=c(0.1,0.1)) +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
dev.off()

labs = c("CREBBP", "EP300", "BRD2", "BRD4", "HDAC1", "HDAC2", "SMARCB1", "ARID1A", "SMARCA4", "SMARCA2", "SMARCC1", "SMARCE1", "MAML1", "RBPJ", "SOX2")
pdf(file=paste0(runid, "_EPSv2_step2_small.pdf"), width=3.5, height=4)
ggplot(data = plotdata2 %>% arrange(issig), mapping = aes(x = category, y = `q.significance`, color = issig)) +
    ggrastr::rasterise(geom_point(size=2, stroke=0.1, shape=16), dpi=400) +
	ggrastr::rasterise(geom_point(data = plotdata2 %>% filter(category %in% labs), color = "black", fill = "red", size = 2, stroke=0.3, shape = 21), dpi=400) +
	geom_label_repel(data = plotdata2 %>% filter(category %in% labs), aes(label = category), size = 6, max.overlaps = Inf, force = 3, segment.size = 0.1) +
	scale_color_manual(values = c("grey70", "red"), na.value = "black") +
	scale_x_discrete(expand=c(0.1,0.1)) +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
dev.off()

data.frame(en) %>% drop_na() %>% arrange(-`q.significance`) %>% filter(category %in% patient_EPS_factors) %>% pull(category) %>% cat()

write.table(data.frame(en) %>% drop_na() %>% arrange(-`q.significance`) %>% filter(category %in% rownames(filtered$RNA)) %>% left_join(NeoExprSignature, by = c("category" = "gene")) %>% mutate(rank = row_number()), "Run10_flat_EPS_RemapEnrich.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(data.frame(en) %>% drop_na() %>% arrange(-`q.significance`) %>% filter(category %in% patient_EPS_factors) %>% mutate(rank = row_number()) %>% dplyr::select(rank, category), "Run10_flat_EPS_final.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########################### Plot Persister Scores FeaturePlots and ViolinPlots #########################################

hm_epithelial <- readRDS(paste0(runid, "_epithelial.rds"))
hm_epithelial_malig <- subset(hm_epithelial, subset = group == "Tumor")

Tumor_Epi$epicluster <- data.frame(cell = colnames(Tumor_Epi)) %>% left_join(hm_epithelial_malig@meta.data %>% mutate(cell = colnames(hm_epithelial_malig))) %>% pull(epicluster)
Tumor_Epi[["harmony_pca"]] <- hm_epithelial_malig[["harmony_pca"]]

DefaultAssay(Tumor_Epi) <- "RNA"
Tumor_Epi %<>% RunUMAP(reduction = "harmony_pca", dims = 1:20, n.neighbors = 20, min.dist = 0.2, reduction.name = "hmrnaumap", reduction.key = "hmrnaumap_", metric = "euclidean")

Tumor_Epi <- AddModuleScore(
  object = Tumor_Epi,
  features = list(patient_EPS_factors),
  assay = "remap",
  name = 'EPSv2_factors',
  ctrl = 10,
  seed = 123
)

pdf(file=paste0(runid, '_harmony_EPSv2_FeaturePlot_patients.pdf'), height = 8, width = 8)
FeaturePlot(object = Tumor_Epi, features = "EPSv2_factors1", min.cutoff = "q5", max.cutoff = "q95", reduction = "hmrnaumap", pt.size = 0.5, order = TRUE)
dev.off()

Tumor_Epi$persist <- case_when(Tumor_Epi$EPSv2_factors1 > quantile(Tumor_Epi$EPSv2_factors1)[4] ~ "High", Tumor_Epi$EPSv2_factors1 < quantile(Tumor_Epi$EPSv2_factors1)[2] ~ "Low", TRUE ~ "Mid")

plotdata <- Tumor_Epi@meta.data %>% 
	dplyr::select(persist, sample, subgroup3) %>%
	group_by(persist, sample, subgroup3) %>%
	mutate(persist = factor(persist, levels = c("Low", "Mid", "High"))) %>%
	summarize(n = n())
pdf(file=paste0(runid, "_EPSv2_High_Low_props_patients.pdf"), width=6, height=7)
p1 <- ggplot(plotdata %>% group_by(persist) %>% summarize(n = sum(n)), aes(fill=persist, y=n, x="ident")) +
    geom_bar(position="fill", stat="identity", color = "black") + 
	scale_fill_manual(values = c(Low = "#fdb515", Mid = "#cccccb", High = "#3953a4")) +
	theme_cowplot() + theme(axis.text.x = element_blank(), legend.position = "none", axis.title = element_blank())
p2 <- ggplot(plotdata, aes(fill=persist, y=n, x=sample)) +
    geom_bar(position="fill", stat="identity", color = "black") + 
	facet_grid(cols = vars(subgroup3), scales = "free", space = "free") +
	scale_fill_manual(values = c(Low = "#fdb515", Mid = "#cccccb", High = "#3953a4")) +
	theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank())
p1 + p2 + plot_layout(widths = c(1, 4))
dev.off()

saveRDS(Tumor_Epi, "Tumor_Epi.rds")

table1 <- Tumor_Epi@meta.data %>% dplyr::select(sample, subgroup3, EPSv2_factors1, epicluster) %>%
	filter(EPSv2_factors1 != Inf & EPSv2_factors1 != -Inf)

pdf(file=paste0(runid, "_EPSv2_VlnPlots_patients.pdf"), width=14, height=6)
p1 <- ggplot(table1, aes(x=subgroup3, y=EPSv2_factors1)) + 
    geom_violin(colour = "black", fill = "grey80") +
	stat_compare_means(aes(label = ..p.signif..), comparisons = list(c("Sensitive", "Resistant"), c("Resistant", "NACT"))) +
	# ggrastr::rasterise(geom_jitter(size=0.001, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
	stat_summary(fun=mean, geom="point", size=1, color="black") +
    # coord_cartesian(ylim = c(-4, 12)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- ggplot(table1, aes(x=sample, y=EPSv2_factors1)) + 
    geom_violin(colour = "black", fill = "grey80") +
	# ggrastr::rasterise(geom_jitter(size=0.1, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
	stat_summary(fun=mean, geom="point", size=1, color="black") +
    facet_grid(cols = vars(subgroup3), scales = "free", space = "free") +
	# coord_cartesian(ylim = c(-4, 11)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 <- ggplot(table1, aes(x=epicluster, y=EPSv2_factors1)) + 
    geom_violin(colour = "black", fill = "grey80") +
	# ggrastr::rasterise(geom_jitter(size=0.001, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
	stat_summary(fun=mean, geom="point", size=1, color="black") +
    # coord_cartesian(ylim = c(-5, 11)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1 + p2 + p3 + plot_layout(widths = c(2, 4, 3))
dev.off()

############################ Correlation analysis of DBF enrichments (Cooperativity) ##################################

DBFs <- sort(rownames(filtered$remap))
cleanlist <- sort(base::intersect(DBFs, rownames(filtered$RNA)))
filtered$remap$data[is.na(filtered$remap$data)] <- 0
filtered$remap$data[filtered$remap$data == Inf] <- 0

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
smooth_scores <- function(obj){
    raw_scores <- obj$remap$data[cleanlist,]
    cellkNN <- FNN::get.knn(obj[["pca"]][[]],k=50)$nn.index
    rownames(cellkNN) <- colnames(raw_scores)
    sc <- FigR::smoothScoresNN(NNmat=cellkNN,mat=raw_scores,nCores=16)
    return(sc %>% t() %>% data.frame())
}

# FT_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & subgroup2 == "FT"))
# Naive_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & subgroup2 == "Naive"))
# NACT_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & subgroup2 == "NACT"))
# Sens_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & subgroup3 == "Sensitive"))
# Resi_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & subgroup3 == "Resistant"))

# forcor_FT_Epi <- smooth_scores(FT_Epi)
# forcor_Naive_Epi <- smooth_scores(Naive_Epi)
# forcor_NACT_Epi <- smooth_scores(NACT_Epi)
# forcor_Sens_Epi <- smooth_scores(Sens_Epi)
# forcor_Resi_Epi <- smooth_scores(Resi_Epi)
# save(forcor_FT_Epi, forcor_Naive_Epi, forcor_NACT_Epi, forcor_Sens_Epi, forcor_Resi_Epi, file="Run10_remap_smooth_groups_ALL.RData")
load("../Run10/Run10_remap_smooth_groups_ALL.RData")

# pearsons_FT_Epi <- forcor_FT_Epi %>% do(data.frame(t(cor(., .))))
pearsons_Naive_Epi <- forcor_Naive_Epi %>% do(data.frame(t(cor(., .))))
pearsons_NACT_Epi <- forcor_NACT_Epi %>% do(data.frame(t(cor(., .))))
pearsons_Sens_Epi <- forcor_Sens_Epi %>% do(data.frame(t(cor(., .))))
pearsons_Resi_Epi <- forcor_Resi_Epi %>% do(data.frame(t(cor(., .))))

labs <- c()

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
plot_corr_heatmap <- function(corrs, label){
	dists <- as.dist(1 - corrs)
	tree <- hclust(dists, method="complete")
	dend <- as.dendrogram(tree)
	orderedcorrs <- corrs[order.dendrogram(dend), order.dendrogram(dend)]
    color.scheme <- rev(brewer.pal(8,"RdBu"))
	ha = HeatmapAnnotation(isSig = (colnames(orderedcorrs) %in% patient_EPS_factors), col = list(isSig = c(`TRUE` = "black", `FALSE` = "white")), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
	# ra = rowAnnotation(isSig = (colnames(orderedcorrs) %in% patient_EPS_factors), col = list(isSig = c(`TRUE` = "black", `FALSE` = "white")), simple_anno_size = unit(5, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
	# For test
	# ra = rowAnnotation(link = anno_mark(at = match(labs, colnames(orderedcorrs)), labels = labs, which = "row", side = "left"))
	pdf(file=paste0("Run10_remap_Pearson_", label, ".pdf"), width=4.5, height=4)
		# For Figure
		plot(Heatmap(matrix=orderedcorrs, show_row_dend=FALSE, show_column_dend=FALSE, cluster_rows=FALSE, cluster_columns=FALSE, col=col_fun, show_row_names=FALSE, show_column_names=FALSE, show_heatmap_legend=FALSE, left_annotation=rowAnnotation(link = anno_mark(at = match(labs, colnames(orderedcorrs)), labels = labs, which = "row", side = "left", labels_gp = gpar(fontsize = 60))), top_annotation=ha, use_raster = TRUE, raster_by_magick = TRUE))
		# For test
		# plot(Heatmap(matrix=orderedcorrs, show_row_dend=FALSE, show_column_dend=FALSE, cluster_rows=FALSE, cluster_columns=FALSE, col = col_fun, show_row_names=FALSE, show_column_names=FALSE, show_heatmap_legend=FALSE, top_annotation = ha, left_annotation = ra, use_raster = TRUE))
	dev.off()
}

plot_corr_heatmap(pearsons_Naive_Epi, "ALL_Naive_Epi")
plot_corr_heatmap(pearsons_NACT_Epi, "ALL_NACT_Epi")
plot_corr_heatmap(pearsons_Sens_Epi, "ALL_Sens_Epi")
plot_corr_heatmap(pearsons_Resi_Epi, "ALL_Resi_Epi")

# PDX_Epi <- readRDS(paste0(runid, "_epithelial_PDX.rds"))
# PDX_Epi$remap$data[is.na(PDX_Epi$remap$data)] <- 0
# PDX_Epi$remap$data[PDX_Epi$remap$data == Inf] <- 0
# Control_Epi <- subset(PDX_Epi, group == "PDX_C")
# Treated_Epi <- subset(PDX_Epi, group == "PDX_T")
# PH27C_Epi <- subset(PDX_Epi, (patient == "PH27" & group == "PDX_C"))
# PH626C_Epi <- subset(PDX_Epi, (patient == "PH626" & group == "PDX_C"))
# PH27T_Epi <- subset(PDX_Epi, (patient == "PH27" & group == "PDX_T"))
# PH626T_Epi <- subset(PDX_Epi, (patient == "PH626" & group == "PDX_T"))

# forcor_Control_Epi <- smooth_scores(Control_Epi)
# forcor_Treated_Epi <- smooth_scores(Treated_Epi)
# forcor_PH27C_Epi <- smooth_scores(PH27C_Epi)
# forcor_PH626C_Epi <- smooth_scores(PH626C_Epi)
# forcor_PH27T_Epi <- smooth_scores(PH27T_Epi)
# forcor_PH626T_Epi <- smooth_scores(PH626T_Epi)
# save(forcor_Control_Epi, forcor_Treated_Epi, forcor_PH27C_Epi, forcor_PH626C_Epi, forcor_PH27T_Epi, forcor_PH626T_Epi, file="Run10_remap_smooth_groups_PDX_ALL.RData")
load("Run10_remap_smooth_groups_PDX_ALL.RData")

pearsons_PH27C_Epi <- forcor_PH27C_Epi %>% do(data.frame(t(cor(., .))))
pearsons_PH626C_Epi <- forcor_PH626C_Epi %>% do(data.frame(t(cor(., .))))
pearsons_PH27T_Epi <- forcor_PH27T_Epi %>% do(data.frame(t(cor(., .))))
pearsons_PH626T_Epi <- forcor_PH626T_Epi %>% do(data.frame(t(cor(., .))))

plot_corr_heatmap(pearsons_PH27C_Epi, "ALL_PH27C_Epi_cc")
plot_corr_heatmap(pearsons_PH626C_Epi, "ALL_PH626C_Epi_cc")
plot_corr_heatmap(pearsons_PH27T_Epi, "ALL_PH27T_Epi_cc")
plot_corr_heatmap(pearsons_PH626T_Epi, "ALL_PH626T_Epi_cc")

############# Correlation for FT vs. Naive ################

factors <- sort(c("ASCL1", "BCL11B", "BNC2", "CEBPA", "CEBPD", "CREB5", "CTBP2", "CTCFL", "DNMT1", "E2F1", "E2F2", "E2F3", "E2F7", "E2F8", "EBF2", "EGR1", "EGR2", "EGR3", "ETS2", "ETV4", "FLI1", "FOS", "FOSB", "FOSL1", "FOXM1", "FOXP2", "GLIS3", "GRHL1", "GRHL2", "HCFC1R1", "HES1", "HMGA1", "HMGA2", "HNF1A", "HNF1B", "HOXB13", "HOXC5", "HOXC6", "IRF1", "IRF4", "IRF8", "IRX3", "JUN", "JUNB", "KLF2", "KLF4", "KMT2C", "LHX1", "LHX4", "LIN28B", "MAFF", "MAML1", "MBD3", "MCM3", "MECOM", "MEF2C", "MEOX1", "MSX2", "MTF2", "MYBL2", "MYCN", "MYOCD", "NCOA2", "NCOA3", "NFATC1", "NFIA", "NFIB", "NR1H2", "NR4A1", "NR4A3", "NR5A2", "NRIP1", "ONECUT2", "PAX2", "PAX7", "PAX8", "PCGF2", "PDX1", "PGR", "PHF19", "PHOX2A", "POU3F1", "POU6F2", "PRDM1", "PRRX2", "SALL4", "SIX2", "SMARCD3", "SNAI3", "SOX11", "SOX17", "SOX2", "SOX21", "SOX3", "SP6", "STAT4", "TAL1", "TCF7", "TCF7L1", "TCFL5", "TFAP2A", "THRB", "TLE3", "TP73", "TRPS1", "VAX2", "WDHD1", "WT1", "YY1AP1", "ZEB2", "ZFP42", "ZFP57", "ZFY", "ZNF124", "ZNF16", "ZNF19", "ZNF300", "ZNF385D", "ZNF439", "ZNF462", "ZNF492", "ZNF519", "ZNF595", "ZNF675", "ZNF680", "ZNF684", "ZNF695", "ZNF707", "ZNF730", "ZNF93", "ZXDC"))
factors <- intersect(factors, rownames(filtered$remap))

smooth_scores <- function(obj){
    raw_scores <- obj$remap$data[factors,]
    cellkNN <- FNN::get.knn(obj[["pca"]][[]],k=50)$nn.index
    rownames(cellkNN) <- colnames(raw_scores)
    sc <- FigR::smoothScoresNN(NNmat=cellkNN,mat=raw_scores,nCores=16)
    return(sc %>% t() %>% data.frame())
}

# FT_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & subgroup2 == "FT"))
# Naive_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & subgroup2 == "Naive"))

# forcor_FT_Epi <- smooth_scores(FT_Epi)
# forcor_Naive_Epi <- smooth_scores(Naive_Epi)
# save(forcor_FT_Epi, forcor_Naive_Epi, file="Run10_remap_smooth_fig2.RData")
load("Run10_remap_smooth_fig2.RData")

pearsons_FT_Epi <- forcor_FT_Epi %>% do(data.frame(t(cor(., .))))
pearsons_Naive_Epi <- forcor_Naive_Epi %>% do(data.frame(t(cor(., .))))
to_remove <- union(rownames(pearsons_FT_Epi)[which(is.na(pearsons_FT_Epi[1,]))], rownames(pearsons_Naive_Epi)[which(is.na(pearsons_Naive_Epi[1,]))])
forcor_FT_Epi <- forcor_FT_Epi[,setdiff(colnames(forcor_FT_Epi), to_remove)]
forcor_Naive_Epi <- forcor_Naive_Epi[,setdiff(colnames(forcor_Naive_Epi), to_remove)]
pearsons_FT_Epi <- forcor_FT_Epi %>% do(data.frame(t(cor(., .))))
pearsons_Naive_Epi <- forcor_Naive_Epi %>% do(data.frame(t(cor(., .))))

labs <- c("E2F1", "E2F3", "E2F7", "E2F8", "FOXM1", "MYCN", "PGR", "NFIA", "TP73", "JUN", "FOS", "FOSL1", "FOSB", "SOX17", "SOX21", "SOX2", "KMT2C", "DNMT1", "PAX8", "WT1", "MECOM", "ZNF462", "SALL4")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
plot_corr_heatmap_factors_wLabs <- function(pearsons, label){
    mtx <- pearsons %>% ungroup() %>% as.matrix()
    rownames(mtx) <- colnames(mtx)
    corrs <- mtx
    dists <- as.dist(1 - corrs)
    tree <- hclust(dists, method="complete")
    dend <- as.dendrogram(tree)
    pdf(file=paste0("Run10_remap_Pearson_", label, "_factors.pdf"), width=8, height=8)
        plot(Heatmap(matrix=corrs, show_row_dend=FALSE, show_column_dend=TRUE, cluster_rows=dend, cluster_columns=dend, col = col_fun(seq(-1, 1)), show_row_names=FALSE, show_column_names=FALSE, show_heatmap_legend=FALSE, left_annotation=rowAnnotation(link = anno_mark(at = match(labs, colnames(corrs)), labels = labs, which = "row", side = "left")), row_dend_reorder = FALSE, column_dend_reorder = FALSE))
    dev.off()
}
plot_corr_heatmap_factors_wLabs(pearsons_FT_Epi, "FT_Epi")
plot_corr_heatmap_factors_wLabs(pearsons_Naive_Epi, "Naive_Epi")

##################################### Croft et al. data analysis - Figure 5 #######################################

filtered_c <- readRDS("Chromatin_PapercombinedPeaks_share.rds")

filtered_chromvar_pre <- subset(filtered_c, subset = chemo %in% "pre")
filtered_epithelial_pre <- subset(filtered_chromvar_pre, subset = highLevelType %in% "Epithelial")
filtered_epithelial_no_A <- subset(
  filtered_epithelial_pre,
  subset = sample_label != "pre_A")
filtered_epithelial_no_A@meta.data$sample_label %>% table()

DefaultAssay(filtered_epithelial_no_A) <- "remap"
Idents(filtered_epithelial_no_A) <- "sample_label"

filtered_epithelial_no_A <- AddModuleScore(
  object = filtered_epithelial_no_A,
  features = list(patient_EPS_factors),
  assay = "remap",
  name = 'EPSv2_factors',
  ctrl = 10,
  seed = 123
)

table1 <- filtered_epithelial_no_A@meta.data %>% dplyr::select(EPSv2_factors1, sample_label) %>%
  dplyr::filter(EPSv2_factors1 != Inf & EPSv2_factors1 != -Inf)

table1$sample_label <- factor(table1$sample_label, levels = c("pre_C", "pre_E", "pre_B", "pre_D"))

pdf(file=paste0(runid, "_EPSv2_VlnPlots_Croft.pdf"), width=5, height=6)
ggplot(table1, aes(x=sample_label, y=EPSv2_factors1)) + 
  geom_violin(colour = "black") +
  ggrastr::rasterise(geom_jitter(size=0.01, stroke=0.1, shape=16), dpi=400) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  stat_summary(fun=mean, geom="point", size=0.7, color="#F61A23") +
  theme_cowplot() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

library(harmony)
set.seed(1)
hm_epithelial <- filtered_epithelial_no_A %>% RunHarmony(group.by.vars = 'sample_label', assay.use = 'RNA', reduction = 'harmony', dims.use = 1:50, reduction.save = 'harmony_pca', project.dim = FALSE, plot_convergence = FALSE, lambda=0.8)
DefaultAssay(hm_epithelial) <- "RNA"
hm_epithelial %<>% RunUMAP(reduction = "harmony_pca", dims = 1:20, n.neighbors = 20, min.dist = 0.2, reduction.name = "hmrnaumap", reduction.key = "hmrnaumap_", metric = "euclidean") 
hm_epithelial %<>% FindNeighbors(reduction = "harmony_pca", dims = 1:20, annoy.metric = "euclidean", k.param = 20)
hm_epithelial %<>% FindClusters(resolution = 0.18, verbose = FALSE)

pdf(file=paste0('Run10_harmony_EPSv2_FeaturePlot_Croft.pdf'), height = 4, width = 5)
FeaturePlot(object = hm_epithelial, features = "EPSv2_factors1", min.cutoff = "q5", max.cutoff = "q95", reduction = "hmrnaumap", pt.size = 0.5, order = TRUE) + theme(plot.title = element_blank())
dev.off()

pdf(file=paste0('Run10_harmony_EPSv2_samples_Croft.pdf'), height = 4, width = 5)
DimPlot(object = hm_epithelial, reduction = "hmrnaumap", pt.size = 0.5)
dev.off()

################################# Explore targets of signature (FigR analysis) - Figure 6 #########################

signature <- patient_EPS_factors
Epith_Tumor_figR <- read.table("Run10_Epith_Tumor_figR_motif_relaxed_sub.tsv", header = TRUE) %>% 
	set_colnames(c("DORC", "Motif", "Score")) %>%
	mutate(Score = as.numeric(Score))
degsr_NN <- read.table(file = "Run10_degsr_NACT_Naive.tsv", sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells")
degsr_SR <- read.table(file = "Run10_degsr_Resistant_Sensitive.tsv", sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells")

Regulated <- Epith_Tumor_figR %>% 
	filter(Motif %in% signature) %>% 
	mutate(Reg = case_when(Score > 0 ~ "Act", TRUE ~ "Rep")) %>%
	group_by(DORC, Reg) %>%
	summarise(n = n()) %>%
	pivot_wider(names_from = "Reg", values_from = "n") %>%
	replace_na(list(Act = 0, Rep = 0)) %>%
	mutate(PropAct = Act/(Act + Rep))

degsr_sig <- degsr_NN %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj)) %>%
	left_join(Regulated, join_by(gene == DORC))
temp1 <- degsr_sig %>% 
	mutate(Expr = case_when(avg_log2FC > 0 ~ "UP", TRUE ~ "DOWN")) %>% 
	mutate(Reg = case_when(is.na(PropAct) ~ "NoInfo", PropAct == 0.5 ~ "Neutral", PropAct > 0.5 ~ "MoreAct", PropAct < 0.5 ~ "MoreRep"))
DOWN_pie_degsr <- temp1 %>% filter(Expr == "DOWN") %>%
	group_by(Reg) %>%
	summarise(n = n()) %>%
	ggplot(aes(x="", y=n, fill=Reg)) +
	geom_bar(stat="identity", width=1, color="black") +
	coord_polar("y", start=0) +
	scale_fill_manual(values = c(`MoreAct` = "#1a9850", `MoreRep` = "#d73027", `Neutral` = "#FEAE2F", `NoInfo` = "grey80")) +
	theme_void() + 
	theme(legend.position="none")
UP_pie_degsr <- temp1 %>% filter(Expr == "UP") %>%
	group_by(Reg) %>%
	summarise(n = n()) %>%
	ggplot(aes(x="", y=n, fill=Reg)) +
	geom_bar(stat="identity", width=1, color="black") +
	coord_polar("y", start=0) +
	scale_fill_manual(values = c(`MoreAct` = "#1a9850", `MoreRep` = "#d73027", `Neutral` = "#FEAE2F", `NoInfo` = "grey80")) +
	theme_void() + 
	theme(legend.position="none")

pdf(file="Run10_degsr_NACT_Naive_targets.pdf", width=8, height=8)
p1 <- ggplot(degsr_sig, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(is.na(PropAct)), size = 1.5, stroke=0.1, shape = 16, color = "grey80"), dpi = 400) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(!is.na(PropAct)), mapping = aes(color = PropAct), size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	# scale_color_distiller(palette = "RdYlGn", direction = 1) +
	scale_color_gradient2(low = "#d73027", mid = "#FEAE2F", high = "#1a9850", midpoint = 0.5) +
	theme_cowplot()
print((p1 + annotation_custom(ggplotGrob(DOWN_pie_degsr), xmin = -10, xmax = -4, ymin = 190, ymax = 250) + annotation_custom(ggplotGrob(UP_pie_degsr), xmin = 4, xmax = 10, ymin = 190, ymax = 250)), newpage = FALSE)
dev.off()

pdf(file="Run10_degsr_NACT_Naive_targets1.pdf", width=3, height=3)
p1 <- ggplot(degsr_sig, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(is.na(PropAct)), size = 2, stroke=0.1, shape = 16, color = "grey80"), dpi = 400) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(!is.na(PropAct)), mapping = aes(color = PropAct), size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
	scale_color_gradient2(low = "#d73027", mid = "#FEAE2F", high = "#1a9850", midpoint = 0.5) +
	theme_cowplot() + theme(legend.position = "none")
p1
dev.off()

degsr_sig <- degsr_SR %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj)) %>%
	left_join(Regulated, join_by(gene == DORC))
temp1 <- degsr_sig %>% 
	mutate(Expr = case_when(avg_log2FC > 0 ~ "UP", TRUE ~ "DOWN")) %>% 
	mutate(Reg = case_when(is.na(PropAct) ~ "NoInfo", PropAct == 0.5 ~ "Neutral", PropAct > 0.5 ~ "MoreAct", PropAct < 0.5 ~ "MoreRep"))
DOWN_pie_degsr <- temp1 %>% filter(Expr == "DOWN") %>%
	group_by(Reg) %>%
	summarise(n = n()) %>%
	ggplot(aes(x="", y=n, fill=Reg)) +
	geom_bar(stat="identity", width=1, color="black") +
	coord_polar("y", start=0) +
	scale_fill_manual(values = c(`MoreAct` = "#1a9850", `MoreRep` = "#d73027", `Neutral` = "#FEAE2F", `NoInfo` = "grey80")) +
	theme_void() + 
	theme(legend.position="none")
UP_pie_degsr <- temp1 %>% filter(Expr == "UP") %>%
	group_by(Reg) %>%
	summarise(n = n()) %>%
	ggplot(aes(x="", y=n, fill=Reg)) +
	geom_bar(stat="identity", width=1, color="black") +
	coord_polar("y", start=0) +
	scale_fill_manual(values = c(`MoreAct` = "#1a9850", `MoreRep` = "#d73027", `Neutral` = "#FEAE2F", `NoInfo` = "grey80")) +
	theme_void() + 
	theme(legend.position="none")

pdf(file="Run10_degsr_Resistant_Sensitive_targets.pdf", width=8, height=8)
p1 <- ggplot(degsr_sig, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(is.na(PropAct)), size = 1.5, stroke=0.1, shape = 16, color = "grey80"), dpi = 400) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(!is.na(PropAct)), mapping = aes(color = PropAct), size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	# scale_color_distiller(palette = "RdYlGn", direction = 1) +
	scale_color_gradient2(low = "#d73027", mid = "#FEAE2F", high = "#1a9850", midpoint = 0.5) +
	theme_cowplot()
print((p1 + annotation_custom(ggplotGrob(DOWN_pie_degsr), xmin = -10, xmax = -4, ymin = 190, ymax = 250) + annotation_custom(ggplotGrob(UP_pie_degsr), xmin = 4, xmax = 10, ymin = 190, ymax = 250)), newpage = FALSE)
dev.off()

pdf(file="Run10_degsr_Resistant_Sensitive_targets1.pdf", width=3, height=3)
p1 <- ggplot(degsr_sig, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(is.na(PropAct)), size = 2, stroke=0.1, shape = 16, color = "grey80"), dpi = 400) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(!is.na(PropAct)), mapping = aes(color = PropAct), size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
	scale_color_gradient2(low = "#d73027", mid = "#FEAE2F", high = "#1a9850", midpoint = 0.5) +
	theme_cowplot() + theme(legend.position = "none")
p1
dev.off()

degsr_sig <- degsr_NN %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj)) %>%
	left_join(Regulated, join_by(gene == DORC))
temp1 <- degsr_sig %>% 
	mutate(Expr = case_when(avg_log2FC > 0 ~ "UP", TRUE ~ "DOWN")) %>% 
	mutate(Reg = case_when(is.na(PropAct) ~ "NoInfo", PropAct == 0.5 ~ "Neutral", PropAct > 0.5 ~ "MoreAct", PropAct < 0.5 ~ "MoreRep"))

rightgreengenes <- temp1 %>% filter(Expr == "UP" & Reg == "MoreAct") %>% pull(gene)
leftredgenes <- temp1 %>% filter(Expr == "DOWN" & Reg == "MoreRep") %>% pull(gene)

gostres1 <- gost(query = leftredgenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))
gostres2 <- gost(query = rightgreengenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))

# highlight_terms <- c("response to stimulus", "developmental process", "cell communication", "Transcriptional misregulation in cancer", "Pathways in cancer", "Interleukin-4 and Interleukin-13 signaling", "cell division", "DNA repair", "Cell Cycle")
highlight_terms <- unique(c("regulation of cell communication", "cell communication", "regulation of cell adhesion", "cell migration", "regulation of cell migration", "cellular response to oxygen", "response to oxygen containing compound", "response to stress", "regulation of response to stress", "regulation of cell cycle", "vasculature development", "blood vessel development", "response to hormone", "RHO GTPase cycle", "Signaling by Rho GTPases", "Signaling by Receptor Tyrosine Kinases", "Endocytosis", "Tight junction", "response to stimulus", "developmental process", "cell communication", "Transcriptional misregulation in cancer", "Pathways in cancer", "Interleukin-4 and Interleukin-13 signaling", "cell division", "DNA repair", "Cell Cycle"))
highlight_terms <- c("response to stress", "response to stimulus", "Transcriptional misregulation in cancer", "Interleukin−4 and Interleukin−13 signaling", "Interferon alpha/beta signaling", "cell cycle", "Cell Cycle", "Transcriptional Regulation by TP53", "Regulation of TP53 Activity", "Polycomb repressive complex", "regulation of primary metabolic process", "Cellular responses to stress", "Cellular responses to stimuli")

# highlight_terms <- gostres1$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres1$result <- gostres1$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))
# highlight_terms <- gostres2$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres2$result <- gostres2$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))

pdf(file="Run10_degsr_NACT_Naive_targets_GOST.pdf", width=8, height=8)
p1 <- gostplot_MOD(gostres1, capped = FALSE, interactive = FALSE)
p2 <- gostplot_MOD(gostres2, capped = FALSE, interactive = FALSE)
p1 | p2
dev.off()

degsr_sig <- degsr_SR %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
    mutate(neg_log10_adj_pval = -log10(p_val_adj)) %>%
	left_join(Regulated, join_by(gene == DORC))
temp1 <- degsr_sig %>% 
	mutate(Expr = case_when(avg_log2FC > 0 ~ "UP", TRUE ~ "DOWN")) %>% 
	mutate(Reg = case_when(is.na(PropAct) ~ "NoInfo", PropAct == 0.5 ~ "Neutral", PropAct > 0.5 ~ "MoreAct", PropAct < 0.5 ~ "MoreRep"))

rightgreengenes <- temp1 %>% filter(Expr == "UP" & Reg == "MoreAct") %>% pull(gene)
leftredgenes <- temp1 %>% filter(Expr == "DOWN" & Reg == "MoreRep") %>% pull(gene)

gostres1 <- gost(query = leftredgenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))
gostres2 <- gost(query = rightgreengenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))

# highlight_terms <- c("response to stimulus", "developmental process", "cell communication", "Transcriptional misregulation in cancer", "Pathways in cancer", "Interleukin-4 and Interleukin-13 signaling", "cell division", "DNA repair", "Cell Cycle")
# highlight_terms <- gostres1$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres1$result <- gostres1$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))
# highlight_terms <- gostres2$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres2$result <- gostres2$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))

pdf(file="Run10_degsr_Resistant_Sensitive_targets_GOST.pdf", width=8, height=8)
p1 <- gostplot_MOD(gostres1, capped = FALSE, interactive = FALSE)
p2 <- gostplot_MOD(gostres2, capped = FALSE, interactive = FALSE)
p1 | p2
dev.off()

######################### Cutsite Dist (Supplementary) ################################

temp2 <- data.frame(cutsites = colSums(filtered[["mATAC"]]@data), celltype = filtered$customtype2, status1 = filtered$subgroup2, status2 = filtered$subgroup3, status3 = filtered$subgroup5) %>%
    mutate(status1 = factor(status1, levels=c("FT", "Naive", "NACT")), status2 = factor(status2, levels=c("BRCAwt", "BRCAmt", "Sensitive", "Resistant", "NACT")), status3 = factor(status3))
temp3 <- temp2 %>% group_by(celltype, status1) %>% summarise(n = n())
temp4 <- temp2 %>% group_by(celltype, status2) %>% summarise(n = n())
temp5 <- temp2 %>% group_by(celltype, status3) %>% summarise(n = n())

pdf(file='Run10_cutsite-dist.pdf', width=16, height=9)
ggplot(temp2, aes(x = status1, y = cutsites, color = celltype)) +
    ggrastr::rasterise(geom_jitter(size=0.5, stroke=0.1, shape=16), dpi=400) +
    geom_violin() +
    stat_summary(fun=mean, geom="point", size=3, color="black", shape=3) +
    geom_text(data = temp3, mapping = aes(x = status1, label = paste0(n)), y = -1000, angle = 45, color = "black") +
    facet_grid(cols = vars(celltype), scales="free") +
    scale_color_manual(values = colors.customtype2) +
    ylab("Number of cutsites") +
    # scale_y_break(c(50000, 350000), ticklabels=c(350000, 360000), expand = expansion(mult = c(0.2, 0))) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.text.x = element_text(size = 10), axis.title.x = element_blank())
dev.off()

temp6 <- temp2 %>% group_by(celltype, status2) %>% summarise(mu = mean(cutsites)) %>% filter(celltype == "Epithelial_cells")
pdf(file='Run10_cutsite-dist_groups.pdf', width=5, height=6)
ggplot(temp2 %>% filter(celltype == "Epithelial_cells" & !(status1 == "FT")), aes(x = status2, y = cutsites, color = celltype)) +
    ggrastr::rasterise(geom_jitter(size=0.5, stroke=0.1, shape=16), dpi=400) +
    geom_violin() +
    stat_summary(fun=mean, geom="point", size=3, color="black", shape=3) +
    geom_text(data = temp4 %>% filter(celltype == "Epithelial_cells" & status2 %in% c("Sensitive", "Resistant", "NACT")), mapping = aes(x = status2, label = paste0(n)), y = -1000, angle = 45, color = "black") +
    geom_text(data = temp6 %>% filter(celltype == "Epithelial_cells" & status2 %in% c("Sensitive", "Resistant", "NACT")), mapping = aes(x = status2, y = mu, label = paste0(round(mu))), nudge_y = 10000, color = "black") +
    facet_grid(cols = vars(celltype), scales="free") +
    scale_color_manual(values = colors.customtype2) +
    ylab("Number of cutsites") +
    # scale_y_break(c(50000, 350000), ticklabels=c(350000, 360000), expand = expansion(mult = c(0.2, 0))) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.text.x = element_text(size = 10), axis.title.x = element_blank())
dev.off()

########################### FigR plots - DEGRN - #############################

# cd /working
# cat Run10_Epith_FT_figR_motif_relaxed.tsv | awk '{ if ($10 >= 1 || $10 <= -1) print $1,$2,$10 }' OFS='\t' > Run10_Epith_FT_figR_motif_sub.tsv
# cat Run10_Epith_Naive_figR_motif_relaxed.tsv | awk '{ if ($10 >= 1 || $10 <= -1) print $1,$2,$10 }' OFS='\t' > Run10_Epith_Naive_figR_motif_sub.tsv
# cat Run10_Epith_NACT_figR_motif_relaxed.tsv | awk '{ if ($10 >= 1 || $10 <= -1) print $1,$2,$10 }' OFS='\t' > Run10_Epith_NACT_figR_motif_sub.tsv
# cat Run10_Epith_PDX_figR_motif_relaxed.tsv | awk '{ if ($10 >= 1 || $10 <= -1) print $1,$2,$10 }' OFS='\t' > Run10_Epith_PDX_figR_motif_sub.tsv

figR_motif_Epith_FT <- read.table("Run10_Epith_FT_figR_motif_sub.tsv", header = TRUE) %>% mutate(Score = as.numeric(Score))
figR_motif_Epith_Naive <- read.table("Run10_Epith_Naive_figR_motif_sub.tsv", header = TRUE) %>% mutate(Score = as.numeric(Score))
figR_motif_Epith_NACT <- read.table("Run10_Epith_NACT_figR_motif_sub.tsv", header = TRUE) %>% mutate(Score = as.numeric(Score))

degsr_Naive_FT <- read.table(file = "Run10_degsr_Naive_FT.tsv", sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% filter(p_val_adj < 0.05)
degsr_NACT_Naive <- read.table(file = "Run10_degsr_NACT_Naive.tsv", sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% filter(p_val_adj < 0.05)

plotDEGRN <- function(ct, reg_tbl, degsr, thresh = 2, label = "", vertex_mult = 1, motifsub = NULL){
	if(length(motifsub) > 0){ reg_tbl <- reg_tbl %>% filter(Motif %in% motifsub) }
    ct_degsr <- degsr %>% filter(celltype == ct)
    if(nrow(ct_degsr) == 0) {return(ggplot() + theme_void())}
    # selectdegs <- ct_degsr %>% pull(gene) %>% head(n = 1000)
    net_tbl_data <- reg_tbl %>% dplyr::select(DORC, Motif, Score) %>%
        filter((abs(Score) > thresh) & (DORC %in% ct_degsr$gene) & (Motif %in% ct_degsr$gene)) %>%
        # group_by(DORC) %>% slice_max(abs(Score), n = 5) %>%
        left_join(ct_degsr %>% dplyr::select(gene, avg_log2FC), by = join_by(DORC == gene)) %>% rename(DORC_FC = avg_log2FC) %>%
        left_join(ct_degsr %>% dplyr::select(gene, avg_log2FC), by = join_by(Motif == gene)) %>% rename(Motif_FC = avg_log2FC) %>%
        filter((Score > 0 & (sign(DORC_FC) * sign(Motif_FC) == 1)) | (Score < 0 & (sign(DORC_FC) * sign(Motif_FC) == -1)))
    net_tbl <- net_tbl_data %>% dplyr::select(Motif, DORC, Score) %>%
        set_colnames(c("from", "to", "Score")) %>%
        filter(from != to)
    if(nrow(net_tbl) == 0) {return(ggplot() + theme_void())}
    net <- network(net_tbl, directed = TRUE)
    net %v% "avg_log2FC" <- (data.frame(gene = (net %v% "vertex.names")) %>% left_join(ct_degsr))$avg_log2FC * vertex_mult
	net %v% "vertex.names" <- case_when((net %v% "vertex.names") %in% motifsub ~ (net %v% "vertex.names"), TRUE ~ "")
    net <- ggnetwork(net, layout = "fruchtermanreingold", niter = 1000, arrow.gap = 0.015)
    # print(paste0("vertices = ", network.size(net)))
    # print(paste0("edges = ", network.edgecount(net)))
    if (label == "") {tag <- ct} else {tag <- label}
    p <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(data = net, aes(color = Score), arrow = arrow(length = unit(4, "pt"), type = "closed")) +
        scale_colour_gradientn(colours = c("darkred",muted("red"),"white","white","white",muted("green"),"darkgreen"), values = rescale(c(-5, -0.81, -0.8, 0, 0.8, 0.81, 5)), limits = c(-5, 5), name = "", oob=squish) +
        new_scale_color() +
        geom_nodes(data = net, aes(color = avg_log2FC), size = 6) +
        geom_nodetext(data = net, aes(label = vertex.names), color = "black", size = 2) +
        scale_color_gradient2(low = "darkred", mid = "white", high = "darkgreen", limits = c(-5, 5), name = "", oob=squish) +
        ggtitle(tag) +
        theme_blank()
    return(p)
}

full_list <- c("KLF6", "JUN", "JUNB", "FOS", "FOSB", "PAX8", "WT1", "FOXM1", "ESR1", "PGR", "RUNX3", "BARX2", "ATF3", "TBX2", "TP73", "RFX2", "STAT4", "RELA", "AR", "ZNF496", "FOXP2")

pdf(file=paste0(runid, "_DEGRNs_FT_Naive_NACT_Epith_degsr_sub.pdf"), width=30, height=20)
plot1 <- plotDEGRN("Epithelial_cells", reg_tbl = figR_motif_Epith_FT, degsr = degsr_Naive_FT, thresh = 1, label = "FT", vertex_mult = -1, motifsub = full_list)
plot2 <- plotDEGRN("Epithelial_cells", reg_tbl = figR_motif_Epith_Naive, degsr = degsr_Naive_FT, thresh = 1, label = "Naive", vertex_mult = 1, motifsub = full_list)
plot3 <- plotDEGRN("Epithelial_cells", reg_tbl = figR_motif_Epith_Naive, degsr = degsr_NACT_Naive, thresh = 1, label = "Naive", vertex_mult = -1, motifsub = full_list)
plot4 <- plotDEGRN("Epithelial_cells", reg_tbl = figR_motif_Epith_NACT, degsr = degsr_NACT_Naive, thresh = 1, label = "NACT", vertex_mult = 1, motifsub = full_list)
(plot1 | plot2 | (ggplot() + theme_void())) / ((ggplot() + theme_void()) | plot3 | plot4)
dev.off()

make_RegPie_plots <- function(figr1, figr2, degsr, outfilename, w=25, h=12){
    degs_list <- unique(degsr$gene)

    figR_motif_1_sub <- figr1 %>% filter(abs(Score) > 1) %>% dplyr::select(DORC, Motif, Score) %>%
        filter(DORC %in% degs_list & Motif %in% degs_list) %>% rename(score1 = Score)
    figR_motif_2_sub <- figr2 %>% filter(abs(Score) > 1) %>% dplyr::select(DORC, Motif, Score) %>%
        filter(DORC %in% degs_list & Motif %in% degs_list) %>% rename(score2 = Score)

    ct_degsr <- degsr %>% filter(celltype == "Epithelial_cells")
    net_tbl_data_1 <- figR_motif_1_sub %>% set_colnames(c("DORC", "Motif", "Score")) %>%
        filter((DORC %in% ct_degsr$gene) & (Motif %in% ct_degsr$gene) & (abs(Score) > 1)) %>%
        left_join(ct_degsr %>% dplyr::select(gene, avg_log2FC), by = join_by(DORC == gene)) %>% rename(DORC_FC = avg_log2FC) %>%
        left_join(ct_degsr %>% dplyr::select(gene, avg_log2FC), by = join_by(Motif == gene)) %>% rename(Motif_FC = avg_log2FC) %>%
        mutate(gexpr_pattern = case_when((Score > 0 & (sign(DORC_FC) * sign(Motif_FC) == 1)) | (Score < 0 & (sign(DORC_FC) * sign(Motif_FC) == -1)) ~ "consistent", TRUE ~ "inconsistent")) %>%
        filter(DORC != Motif)
    summ1 <- net_tbl_data_1 %>% 
        mutate(reg = case_when((Score > 0 & gexpr_pattern == "consistent") ~ "Act", (Score < 0 & gexpr_pattern == "consistent") ~ "Rep", TRUE ~ "inconsistent")) %>%
        group_by(Motif, reg) %>% 
        summarize(n = n()) %>%
        pivot_wider(names_from = "reg", values_from = "n") %>% 
        replace_na(list(Act = 0, Rep = 0)) %>% 
        mutate(total = Act + Rep) %>%
        arrange(-total) %>%
        # filter(total > 10) %>%
        dplyr::select(-total) %>% 
        pivot_longer(!Motif, names_to = "Reg", values_to = "count")
    net_tbl_1 <- net_tbl_data_1 %>% dplyr::select(Motif, DORC, Score) %>%
        set_colnames(c("from", "to", "Score")) %>%
        filter(from != to)
    net_tbl_data_2 <- figR_motif_2_sub %>% set_colnames(c("DORC", "Motif", "Score")) %>%
        filter((DORC %in% ct_degsr$gene) & (Motif %in% ct_degsr$gene) & (abs(Score) > 1)) %>%
        left_join(ct_degsr %>% dplyr::select(gene, avg_log2FC), by = join_by(DORC == gene)) %>% rename(DORC_FC = avg_log2FC) %>%
        left_join(ct_degsr %>% dplyr::select(gene, avg_log2FC), by = join_by(Motif == gene)) %>% rename(Motif_FC = avg_log2FC) %>%
        mutate(gexpr_pattern = case_when((Score > 0 & (sign(DORC_FC) * sign(Motif_FC) == 1)) | (Score < 0 & (sign(DORC_FC) * sign(Motif_FC) == -1)) ~ "consistent", TRUE ~ "inconsistent")) %>%
        filter(DORC != Motif)
    summ2 <- net_tbl_data_2 %>% 
        mutate(reg = case_when((Score > 0 & gexpr_pattern == "consistent") ~ "Act", (Score < 0 & gexpr_pattern == "consistent") ~ "Rep", TRUE ~ "inconsistent")) %>%
        group_by(Motif, reg) %>% 
        summarize(n = n()) %>%
        pivot_wider(names_from = "reg", values_from = "n") %>% 
        replace_na(list(Act = 0, Rep = 0)) %>% 
        mutate(total = Act + Rep) %>%
        arrange(-total) %>%
        # filter(total > 10) %>%
        dplyr::select(-total) %>% 
        pivot_longer(!Motif, names_to = "Reg", values_to = "count")
    net_tbl_2 <- net_tbl_data_2 %>% dplyr::select(Motif, DORC, Score) %>%
        set_colnames(c("from", "to", "Score")) %>%
        filter(from != to)

    TFs_in_both <- intersect(summ1$Motif, summ2$Motif)
    temp <- summ1 %>% mutate(group = "1") %>% bind_rows(summ2 %>% mutate(group = "2")) %>% replace_na(list(count = 0)) %>% filter(Motif %in% TFs_in_both)
    ord3 <- temp %>% filter(group == "2") %>% dplyr::select(-group) %>% group_by(Motif) %>% mutate(per =  100 *count/sum(count)) %>% ungroup() %>% filter(Reg == "Act") %>% arrange(-per, -count) %>% pull(Motif)
    temp1 <- temp %>% filter(Motif %in% ord3) %>% mutate(Motif = factor(Motif, levels = ord3))

    pdf(file=outfilename, width=w, height=h)
    p1 <- ggplot(temp1, aes(fill=Reg, y=count, x=Motif)) + 
        geom_bar(position="stack", stat="identity", color = "black") + 
        facet_grid(rows = vars(group)) +
        scale_fill_manual(values = c("darkgreen", "red", "grey50")) +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    p2 <- ggplot(temp1, aes(fill=Reg, y=count, x=Motif)) + 
        geom_bar(position="fill", stat="identity", color = "black") + 
        facet_grid(rows = vars(group)) +
        scale_fill_manual(values = c("darkgreen", "red", "grey50")) +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    plot(p1 / p2)
    dev.off()
}
make_RegPie_plots(figR_motif_Epith_FT, figR_motif_Epith_Naive, degsr_Naive_FT, "Run10_RegPie_Bar_FT_Naive_Epith.pdf", w=34, h=12)

######################### NUMBAT / CNV ##################################

numbatouts <- paste0(datadir, "numbat_runs/", samplenames[6:14], "/out")
nbs <- lapply(numbatouts, function(x){Numbat$new(out_dir = x)})

posteriors <- mapply(function(x, y){
	data.frame(cell = colnames(filtered), sample = filtered$sample, celltype = filtered$customtype) %>% 
		filter(sample == x) %>% separate(col = "cell", into = c(NA, "cell"), sep = "_") %>% 
		left_join(y$clone_post %>% dplyr::select(cell, p_cnv_x, p_cnv_y, p_cnv, compartment_opt, clone_opt, GT_opt)) %>% drop_na()
}, samplelabels[6:13], nbs[1:8], SIMPLIFY = FALSE)
temp <- data.frame(cell = colnames(filtered), sample = filtered$sample, celltype = filtered$customtype) %>% 
	filter(sample == "PH1303")%>% mutate(cell = str_sub(cell, start = 8L, end = -1L)) %>% 
	left_join(nbs[[9]]$clone_post %>% dplyr::select(cell, p_cnv_x, p_cnv_y, p_cnv, compartment_opt, clone_opt, GT_opt)) %>% drop_na()
posteriors <- bind_rows(bind_rows(posteriors), temp) %>% mutate(sample = factor(sample, levels = samplelabels[6:14])) %>%
	mutate(CNVstatus = case_when(p_cnv > 0.55 ~ "high", p_cnv < 0.45 ~ "low", TRUE ~ "ambiguous")) %>%
	left_join(data.frame(sample = samplelabels, label = names(samplelabels))) %>%
	mutate(cell = paste0(label, "_", cell))

# pdf(file=paste0('Run10_posteriors.pdf'), height = 6, width = 6)
# ggplot(posteriors, aes(p_cnv)) +
	# geom_histogram(binwidth = 0.01)
# dev.off()

plotdata <- filtered@meta.data %>% dplyr::select(customtype2, seurat_clusters, sample) %>% rownames_to_column("cell") %>%
	left_join(posteriors %>% dplyr::select(cell, CNVstatus)) %>%
	replace_na(list(CNVstatus = "NotTested")) %>%
	mutate(CNVstatus = factor(CNVstatus, levels = c("high", "ambiguous", "low", "NotTested")))
plotdata1 <- plotdata %>% group_by(customtype2, sample, CNVstatus) %>% 
	summarize(n = n())

saveRDS(plotdata1, 'Run10_CNVs.rds')
pdf(file=paste0('Run10_CNVs.pdf'), height = 25, width = 5)
ggplot(plotdata1) +
	geom_bar(aes(fill=CNVstatus, x=n, y=sample), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
	facet_grid(rows = vars(customtype2), scales = "free", space = "free") +
	theme_cowplot()
dev.off()

########################### MGATK / Mitochondrial Clonotyping ###############################################

process_mito <- function(samplename, samplelabel){
    variant_table <- read.table(paste0(datadir, samplename, "/outs/mtgt/final/", samplename, ".variant_stats.tsv.gz"), sep="\t", header=TRUE)
    variant_table <- subset(variant_table, subset = n_cells_conf_detected >= 5 & strand_correlation >= 0.65)
    mgatk <- ReadMGATK(dir = paste0(datadir, samplename, "/outs/mtgt/final/"))
    mito <- CreateSeuratObject(counts = mgatk$counts, assay = "mito")
    cells <- filtered@meta.data %>% filter(sample == samplelabel) %>% pull(gex_barcode)
    mito <- subset(mito, cells = cells)
    colnames(mito) <- data.frame(gex_barcode = colnames(mito)) %>% left_join(filtered@meta.data %>% filter(sample == samplelabel) %>% rownames_to_column("bc")) %>% pull(bc)
    mito <- AlleleFreq(mito, variants = variant_table$variant, assay = "mito")
    DefaultAssay(mito) <- "alleles"
    mito <- FindClonotypes(mito)
    return(mito)
}

mitobjs <- mapply(process_mito, samplenames[6:length(samplenames)], samplelabels[6:length(samplelabels)])
names(mitobjs) <- samplelabels[6:length(samplelabels)]
merged_mitobjs <- merge(x = mitobjs[[1]], y = mitobjs[2:length(mitobjs)])
mitobjs <- lapply(mitobjs, FindClonotypes)
mitobjs <- mapply(function(x, sample){
    x$mitocluster <- paste0(sample, "_", x$seurat_clusters)
    return(x)
}, mitobjs, names(mitobjs))
merged_mitobjs <- merge(x = mitobjs[[1]], y = mitobjs[2:length(mitobjs)])

# plots <- lapply(mitobjs, function(obj){
    # return(DoHeatmap(obj, features = VariableFeatures(obj), slot = "data", disp.max = 1) + scale_fill_viridis_c())
# })
# pdf(file=paste0(runid, "_mito_heatmaps1.pdf"), width=15, height=18)
# wrap_plots(plots, ncol=2)
# dev.off()

Tumor <- subset(filtered, subset = group == "Tumor")
Tumor[["alleles"]] <- merged_mitobjs[["alleles"]]
Idents(Tumor) <- "customtype2"
DefaultAssay(Tumor) <- "alleles"
Tumor$mitocluster <- data.frame(cell = colnames(Tumor)) %>% left_join(merged_mitobjs@meta.data %>% rownames_to_column("cell") %>% dplyr::select(cell, mitocluster)) %>% pull(mitocluster)
Tumor$sample <- factor(Tumor$sample, levels = samplelabels[6:length(samplelabels)])
DefaultAssay(Tumor) <- "alleles"
TumorSlim <- DietSeurat(Tumor, assays = "alleles")
TumorSlim@meta.data <- TumorSlim@meta.data %>% dplyr::select("customtype", "sample")

saveRDS(TumorSlim, paste0(runid, "_mito_heatmaps2.rds"))
pdf(file=paste0(runid, "_mito_heatmaps2.pdf"), width=40, height=20)
DoHeatmap(TumorSlim, features = rownames(TumorSlim), group.by = c("customtype", "sample"), slot = "data", disp.max = 1) & scale_fill_viridis_c()
dev.off()

process_mito <- function(samplename, samplelabel){
    variant_table <- read.table(paste0(datadir, samplename, "/outs/mtgt/final/", samplename, ".variant_stats.tsv.gz"), sep="\t", header=TRUE)
    variant_table <- subset(variant_table, subset = n_cells_conf_detected >= 5 & strand_correlation >= 0.65)
    mgatk <- ReadMGATK(dir = paste0(datadir, samplename, "/outs/mtgt/final/"))
    mito <- CreateSeuratObject(counts = mgatk$counts, assay = "mito")
    cells <- filtered_PDX@meta.data %>% filter(sample == samplelabel) %>% pull(gex_barcode)
    mito <- subset(mito, cells = cells)
    colnames(mito) <- data.frame(gex_barcode = colnames(mito)) %>% left_join(filtered_PDX@meta.data %>% filter(sample == samplelabel) %>% rownames_to_column("bc")) %>% pull(bc)
    mito <- AlleleFreq(mito, variants = variant_table$variant, assay = "mito")
    DefaultAssay(mito) <- "alleles"
    mito <- FindClonotypes(mito)
    return(mito)
}

PDX_samplenames <- c("PH27_C", "PH27_T", "PH235_C", "PH235_T", "PH626_C", "PH626_T")
mitobjs <- mapply(process_mito, PDX_samplenames, PDX_samplenames)
names(mitobjs) <- PDX_samplenames
merged_mitobjs <- merge(x = mitobjs[[1]], y = mitobjs[2:length(mitobjs)])
mitobjs <- lapply(mitobjs, FindClonotypes)
mitobjs <- mapply(function(x, sample){
    x$mitocluster <- paste0(sample, "_", x$seurat_clusters)
    return(x)
}, mitobjs, names(mitobjs))
merged_mitobjs <- merge(x = mitobjs[[1]], y = mitobjs[2:length(mitobjs)])

# plots <- lapply(mitobjs, function(obj){
    # return(DoHeatmap(obj, features = VariableFeatures(obj), slot = "data", disp.max = 1) + scale_fill_viridis_c())
# })
# pdf(file=paste0(runid, "_mito_heatmaps1_PDX.pdf"), width=15, height=18)
# wrap_plots(plots, ncol=2)
# dev.off()

filtered_PDX[["alleles"]] <- merged_mitobjs[["alleles"]]
Idents(filtered_PDX) <- "epicluster"
DefaultAssay(filtered_PDX) <- "alleles"
filtered_PDX$mitocluster <- data.frame(cell = colnames(filtered_PDX)) %>% left_join(merged_mitobjs@meta.data %>% rownames_to_column("cell") %>% dplyr::select(cell, mitocluster)) %>% pull(mitocluster)
filtered_PDX$sample <- factor(filtered_PDX$sample, levels = PDX_samplenames)
DefaultAssay(filtered_PDX) <- "alleles"
filtered_PDXSlim <- DietSeurat(filtered_PDX, assays = "alleles")
filtered_PDXSlim@meta.data <- filtered_PDXSlim@meta.data %>% dplyr::select("epicluster", "sample")

saveRDS(filtered_PDXSlim, paste0(runid, "_mito_heatmaps2_PDX.rds"))
pdf(file=paste0(runid, "_mito_heatmaps2_PDX.pdf"), width=40, height=20)
DoHeatmap(filtered_PDXSlim, features = rownames(filtered_PDXSlim), group.by = c("epicluster", "sample"), slot = "data", disp.max = 1) & scale_fill_viridis_c()
dev.off()

########################### CICERO / coaccessibility ###############################
gene_bodies <- read.table("features.tsv", sep = "\t") %>% 
    dplyr::select(V2, V4, V5, V6) %>% set_colnames(c("gene", "chr", "start", "end")) %>%
    tidyr::unite("grange", chr:end, remove = TRUE, sep = "-")
gene_bodies <- gene_bodies %>% filter(!str_starts(gene, "chr")) %>% filter(str_starts(grange, "chr"))
gene_body_granges <- StringToGRanges(gene_bodies$grange)
names(gene_body_granges) <- as.character(gene_bodies$gene)
gene_body_granges$gene_name <- as.character(gene_bodies$gene)

connsfiles <- c("Run10_Ciliated_conns.rds", "Run10_Secretory_conns.rds", "Run10_Cluster1_conns.rds", "Run10_Cluster2_conns.rds", "Run10_Cluster3_conns.rds", "Run10_Cycling_conns.rds")
conns.list <- lapply(connsfiles, function(x){readRDS(x) %>% filter(coaccess > 0.1)})
connsfiles_FT <- c("Run10_Ciliated_FT_conns.rds", "Run10_Secretory_FT_conns.rds")
conns.list_FT <- lapply(connsfiles_FT, function(x){readRDS(x) %>% filter(coaccess > 0.1)})
connsfiles_Tumor <- c("Run10_Ciliated_Tumor_conns.rds", "Run10_Secretory_Tumor_conns.rds", "Run10_Cluster1_Tumor_conns.rds", "Run10_Cluster2_Tumor_conns.rds", "Run10_Cluster3_Tumor_conns.rds", "Run10_Cycling_Tumor_conns.rds")
conns.list_Tumor <- lapply(connsfiles_Tumor, function(x){readRDS(x) %>% filter(coaccess > 0.1)})
cleanup_conns <- function(conns_var){
    conns_temp <- conns_var %>% 
        separate_wider_delim(Peak1, "-", names = c("chr1", "start1", "end1")) %>% 
        separate_wider_delim(Peak2, "-", names = c("chr2", "start2", "end2")) %>% 
        mutate(start1 = as.integer(start1), end1 = as.integer(end1), start2 = as.integer(start2), end2 = as.integer(end2))
    conns_temp <- conns_temp %>% filter(start2 > start1) %>%
        bind_rows(conns_temp %>% filter(start2 < start1) %>%
        mutate(chrtemp = chr1, starttemp = start1, endtemp = end1, chr1 = chr2, start1 = start2, end1 = end2) %>%
        mutate(chr2 = chrtemp, start2 = starttemp, end2 = endtemp) %>% dplyr::select(-c(chrtemp, starttemp, endtemp))) %>%
        tidyr::unite("Peak1", chr1:end1, remove = TRUE, sep = "-") %>% tidyr::unite("Peak2", chr2:end2, remove = TRUE, sep = "-") %>%
        group_by(Peak1, Peak2) %>% summarise(coaccess = mean(coaccess)) %>% ungroup()
    return(conns_temp)
}
conns.list <- lapply(conns.list, cleanup_conns)
conns.list_FT <- lapply(conns.list, cleanup_conns)
conns.list_Tumor <- lapply(conns.list, cleanup_conns)

myConnsPlot <- function(obj, reg, genes, ids, ident_colors, conns_list){
	conns_list <- lapply(conns_list, function(x){x %>% mutate(coaccess = 1)})
	link_plots <- lapply(conns_list, function(x){
		Links(obj) <- ConnectionsToLinks(x)
		link_plot <- LinkPlot(obj, region = reg) + theme(legend.position = "none")
		return(link_plot)
	})
	cov_plot <- CoveragePlot(obj, region = reg, features = genes, annotation = FALSE, peaks = FALSE, links = FALSE, idents = ids, scale.factor = 1e7, ymax = 200) & scale_fill_manual(values = ident_colors)
    expr_plot <- ExpressionPlot(obj, features = genes, assay = "RNA", idents = ids) & scale_fill_manual(values = ident_colors)
    gene_plot <- AnnotationPlot(obj, region = reg)
    return(CombineTracks(plotlist = list(cov_plot, link_plots[[1]], link_plots[[2]], link_plots[[3]], link_plots[[4]], link_plots[[5]], link_plots[[6]], gene_plot), expression.plot = expr_plot, heights = c(12, 2, 2, 2, 2, 2, 2, 1), widths = c(12, 3)))
}
myConnsPlot_short <- function(obj, reg, genes, ids, ident_colors, conns_list){
	conns_list <- lapply(conns_list, function(x){x %>% mutate(coaccess = 1)})
	link_plots <- lapply(conns_list, function(x){
		Links(obj) <- ConnectionsToLinks(x)
		link_plot <- LinkPlot(obj, region = reg) + theme(legend.position = "none")
		return(link_plot)
	})
	cov_plot <- CoveragePlot(obj, region = reg, features = genes, annotation = FALSE, peaks = FALSE, links = FALSE, idents = ids, scale.factor = 1e7, ymax = 200) & scale_fill_manual(values = ident_colors)
    expr_plot <- ExpressionPlot(obj, features = genes, assay = "RNA", idents = ids) & scale_fill_manual(values = ident_colors)
    gene_plot <- AnnotationPlot(obj, region = reg)
    return(CombineTracks(plotlist = list(cov_plot, link_plots[[1]], link_plots[[2]], gene_plot), expression.plot = expr_plot, heights = c(4, 2, 2, 1), widths = c(4, 1)))
}

filterGeneTSSConns <- function(conns, gene){
    TSSg <- FigR::hg38TSSRanges
    names(TSSg) <- as.character(TSSg$gene_name)
    allPeaks <- unique(c(as.character(conns$Peak1), as.character(conns$Peak2)))
    OV <- findOverlaps(query = TSSg, subject = StringToGRanges(allPeaks))
    OVd <- OV %>% as.data.frame() %>% dplyr::rename("Gene"="queryHits","Peak"="subjectHits")
    OVd$Gene <- as.character(TSSg$gene_name)[OVd$Gene]
    OVd$Peak <- as.character(allPeaks)[OVd$Peak]
    TSSpeaks <- (OVd %>% filter(Gene == gene))$Peak
    return(conns %>% filter(Peak1 %in% TSSpeaks | Peak2 %in% TSSpeaks))
}
filterGeneBodyConns <- function(conns, gene){
    allPeaks <- unique(c(as.character(conns$Peak1), as.character(conns$Peak2)))
    OV <- findOverlaps(query = gene_body_granges, subject = StringToGRanges(allPeaks))
    OVd <- OV %>% as.data.frame() %>% dplyr::rename("Gene"="queryHits","Peak"="subjectHits")
    OVd$Gene <- as.character(gene_body_granges$gene_name)[OVd$Gene]
    OVd$Peak <- as.character(allPeaks)[OVd$Peak]
    TSSpeaks <- (OVd %>% filter(Gene == gene))$Peak
    return(conns %>% filter(Peak1 %in% TSSpeaks | Peak2 %in% TSSpeaks))
}
conns_list <- lapply(conns.list, filterGeneBodyConns, gene = "PAX8")

ids <- names(epicolors)
ident_colors <- epicolors
Idents(hm_epithelial) <- "epicluster"
connsPlotAroundGene <- function(obj, conns_list, g, thresh = 0.1){
    grange <- gene_windows %>% filter(gene == g) %>% dplyr::select(grange) %>% unlist() %>% unname()
	conns_list <- lapply(conns_list, function(x){x %>% filter(coaccess > thresh)})
	conns_list <- lapply(conns_list, filterGeneTSSConns, gene = g)
    p1 <- myConnsPlot(obj, grange, c(g), ids, ident_colors, conns_list)
    plot(p1)
}

connsPlotAroundGene_short <- function(obj, conns_list, g, thresh = 0.1){
    grange <- gene_windows %>% filter(gene == g) %>% dplyr::select(grange) %>% unlist() %>% unname()
	conns_list <- lapply(conns_list, function(x){x %>% filter(coaccess > thresh)})
	conns_list <- lapply(conns_list, filterGeneTSSConns, gene = g)
    p1 <- myConnsPlot_short(obj, grange, c(g), ids[1:2], ident_colors[1:2], conns_list)
    plot(p1)
}

gene_windows <- read.table("features.tsv", sep = "\t") %>% 
    dplyr::select(V2, V4, V5, V6) %>% set_colnames(c("gene", "chr", "start", "end")) %>%
    mutate(start = case_when(start - 300000 < 0 ~ 1, TRUE ~ start - 300000), end = end + 300000) %>% 
    tidyr::unite("grange", chr:end, remove = TRUE, sep = "-")

# pdf("Run10_conns_around_genes_MSLN.pdf", width=12, height=12)
# connsPlotAroundGene(hm_epithelial, conns.list, "MSLN", thresh = 0.1)
# dev.off()
# pdf("Run10_conns_around_genes_FOXJ1.pdf", width=12, height=12)
# connsPlotAroundGene(hm_epithelial, conns.list, "FOXJ1", thresh = 0.1)
# dev.off()
# pdf("Run10_conns_around_genes_TOP2A.pdf", width=12, height=12)
# connsPlotAroundGene(hm_epithelial, conns.list, "TOP2A", thresh = 0.1)
# dev.off()

DefaultAssay(hm_epithelial) <- "mATAC"
hm_epithelial_FT <- subset(hm_epithelial, subset = group == "FT")
hm_epithelial_Tumor <- subset(hm_epithelial, subset = group == "Tumor")

pdf("Run10_conns_around_genes_Tumor_MSLN.pdf", width=12, height=12)
connsPlotAroundGene(hm_epithelial_Tumor, conns.list_Tumor, "MSLN", thresh = 0.25)
dev.off()
pdf("Run10_conns_around_genes_Tumor_FOXJ1.pdf", width=12, height=12)
connsPlotAroundGene(hm_epithelial_Tumor, conns.list_Tumor, "FOXJ1", thresh = 0.25)
dev.off()
pdf("Run10_conns_around_genes_Tumor_TOP2A.pdf", width=12, height=12)
connsPlotAroundGene(hm_epithelial_Tumor, conns.list_Tumor, "TOP2A", thresh = 0.1)
dev.off()
pdf("Run10_conns_around_genes_Tumor_CDK1.pdf", width=12, height=12)
connsPlotAroundGene(hm_epithelial_Tumor, conns.list_Tumor, "CDK1", thresh = 0.1)
dev.off()

pdf("Run10_conns_around_genes_FT_MSLN.pdf", width=12, height=4)
connsPlotAroundGene_short(hm_epithelial_FT, conns.list_FT, "MSLN", thresh = 0.25)
dev.off()
pdf("Run10_conns_around_genes_FT_FOXJ1.pdf", width=12, height=4)
connsPlotAroundGene_short(hm_epithelial_FT, conns.list_FT, "FOXJ1", thresh = 0.25)
dev.off()
pdf("Run10_conns_around_genes_FT_TOP2A.pdf", width=12, height=4)
connsPlotAroundGene_short(hm_epithelial_FT, conns.list_FT, "TOP2A", thresh = 0.1)
dev.off()
pdf("Run10_conns_around_genes_FT_CDK1.pdf", width=12, height=4)
connsPlotAroundGene_short(hm_epithelial_FT, conns.list_FT, "CDK1", thresh = 0.1)
dev.off()

pdf("Run10_conns_around_genes_FT_PAX8_GeneExpr.pdf", width=6, height=6)
hm_epithelial_FT <- subset(hm_epithelial, subset = group == "FT")
DefaultAssay(hm_epithelial_FT) <- "RNA"
Idents(hm_epithelial_FT) <- "epicluster"
VlnPlot(hm_epithelial_FT, features = c("PAX8"))
dev.off()

####################### Mutation data ##################################

sampleord = c("PH1199", "PH1230", "PH1302", "PH1239", "PH1243", "PH1303", "PH1214", "PH1238", "PH1252")

clpat = read.csv("clinical_patientdata.csv", header = TRUE) %>%
	mutate(sample = factor(sample, levels = sampleord)) %>%
	mutate(stage = factor(stage), CRS = factor(CRS), HR = factor(HR)) %>%
	mutate(group = factor(group, levels = c("Sensitive", "Resistant", "NACT")))

mtpat = read.csv("mutations_patientdata.csv", header = FALSE) %>%
	mutate(sample = case_match(V1, "Ph1238" ~ "PH1238", "Ph1239" ~ "PH1239", "Ph1243" ~ "PH1243", "Ph1302" ~ "PH1302", .default = V1)) %>%
	mutate(status = case_match(V2, "Likely Benign" ~ "Benign", "Benign" ~ "Benign", "Likely pathogenic/VUS" ~ "Pathogenic", "VUS" ~ "VUS", "Pathogenic" ~ "Pathogenic")) %>%
	mutate(status = factor(status, levels = c("Normal", "Benign", "VUS", "Pathogenic"), ordered = TRUE)) %>%
	rename(gene = V3) %>% dplyr::select(sample, gene, status) %>%
	group_by(sample, gene) %>% summarise(status = max(status)) %>%
	mutate(sample = factor(sample, levels = sampleord)) %>%
	tibble() %>%
	complete(sample, gene, fill = list(status = "Normal"))

mtpdx = read.csv("mutations_PDX.csv", header = FALSE) %>%
	mutate(sample = factor(V1, levels = c("PH27", "PH235", "PH626"))) %>%
	mutate(status = case_match(V2, "Likely_benign" ~ "Benign", "Benign" ~ "Benign", "Benign/Likely_benign" ~ "Benign", "Uncertain_significance" ~ "VUS", "Pathogenic/Likely_pathogenic" ~ "Pathogenic", "Pathogenic" ~ "Pathogenic")) %>%
	mutate(status = factor(status, levels = c("Benign", "VUS", "Pathogenic"), ordered = TRUE)) %>%
	rename(gene = V3) %>% dplyr::select(sample, gene, status) %>%
	group_by(sample, gene) %>% summarise(status = max(status))

hmtheme1 = theme_cowplot() + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.margin = unit(c(0, 0, 0, 0), "cm"), axis.ticks.length = unit(0, "cm"))
hmtheme2 = theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks = element_blank(), axis.line = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.margin = unit(c(0, 0, 0, 0), "cm"), axis.ticks.length = unit(0, "cm"))

pdf(file = "Oncoplot_patientdata.pdf", width = 4, height = 8)
p1 = ggplot(mtpat, aes(sample, gene, fill = status)) + 
    geom_tile(color = "black") + 
	scale_fill_manual(values = c(`Benign` = "#FFF000", `VUS` = "#FF8400", `Pathogenic` = "#E64300", `Normal` = "white")) +
	scale_x_discrete(expand=c(0,0)) +
	scale_y_discrete(expand=c(0,0)) +
	coord_fixed() +
	hmtheme1
p2 = ggplot(clpat, aes(sample, "age", fill = age)) + 
    geom_tile(color = "black") + 
	scale_fill_steps(low = "#56B1F7", high = "#132B43") +
	scale_x_discrete(expand=c(0,0)) +
	scale_y_discrete(expand=c(0,0)) +
	coord_fixed() +
	hmtheme1
p3 = ggplot(clpat, aes(sample, "stage", fill = stage)) + 
    geom_tile(color = "black") + 
	scale_fill_manual(values = c("#A847DB", "#973FC4", "#6F2F91", "#582572")) +
	scale_x_discrete(expand=c(0,0)) +
	scale_y_discrete(expand=c(0,0)) +
	coord_fixed() +
	hmtheme1
p4 = ggplot(clpat, aes(sample, "CRS", fill = CRS)) + 
    geom_tile(color = "black") + 
	scale_fill_manual(values = c("#ABC803", "#5D9D00"), na.value = "grey80") +
	scale_x_discrete(expand=c(0,0)) +
	scale_y_discrete(expand=c(0,0)) +
	coord_fixed() +
	hmtheme1
p5 = ggplot(clpat, aes(sample, "HR", fill = HR)) + 
    geom_tile(color = "black") + 
	scale_fill_manual(values = c("#E90026", "#579BE9"), na.value = "grey80") +
	scale_x_discrete(expand=c(0,0)) +
	scale_y_discrete(expand=c(0,0)) +
	coord_fixed() +
	hmtheme1
p6 = ggplot(clpat, aes(sample, "group", fill = group)) + 
    geom_tile(color = "black") + 
	scale_fill_manual(values = c(`Sensitive` = "#489fa7", `Resistant` = "#c45757", `NACT` = "#306c7b")) +
	scale_x_discrete(expand=c(0,0)) +
	scale_y_discrete(expand=c(0,0)) +
	coord_fixed() +
	hmtheme2
(p1 / p2 / p3 / p4 / p5 / p6) + plot_layout(guides = 'collect')
dev.off()

#  + xlab(NULL) + ylab(NULL)

pdf(file = "Oncoplot_PDX.pdf", width = 3, height = 3)
ggplot(mtpdx, aes(sample, gene, fill = status)) + 
    geom_tile() + 
	scale_fill_manual(values = c(`Benign` = "#FFF000", `VUS` = "#FF8400", `Pathogenic` = "#E64300")) +
    theme_cowplot() +
    xlab(NULL) +
    ylab(NULL) + 
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##################################### Figure 6c #######################################

signature <- patient_EPS_factors

Epith_Tumor_figR <- read.table("Run10_Epith_Tumor_figR_motif_relaxed_sub.tsv", header = TRUE) %>% 
  set_colnames(c("DORC", "Motif", "Score")) %>%
  mutate(Score = as.numeric(Score))

get_gene_summary <- function(gene_list,
                              data = Epith_Tumor_figR,
                              motif_filter = signature,
                              fill_colors = c("MoreAct" = "#1a9850",
                                              "MoreRep" = "#d73027",
                                              "Neutral" = "#FEAE2F",
                                              "NoInfo" = "grey80")) {
  
  total_genes <- length(gene_list)
  
  # Filter and process the data
  temp <- data %>% 
    dplyr::filter(Motif %in% motif_filter) %>% 
    dplyr::filter(DORC %in% gene_list) %>% 
    dplyr::mutate(Reg = factor(case_when(
      Score > 0 ~ "Act",
      TRUE ~ "Rep"
    ), levels = c("Act", "Rep"))) %>%
    group_by(DORC, Reg) %>%
    summarise(n = n(), .groups = "drop") %>%
	complete(DORC, Reg, fill = list(n = 0)) %>%
    pivot_wider(names_from = "Reg", values_from = "n") %>%
    replace_na(list(Act = 0, Rep = 0)) %>%
    dplyr::mutate(PropAct = Act / (Act + Rep)) %>% 
    dplyr::mutate(Reg = case_when(
      is.na(PropAct)       ~ "NoInfo",
      PropAct == 0.5       ~ "Neutral",
      PropAct > 0.5        ~ "MoreAct",
      PropAct < 0.5        ~ "MoreRep"
    ))
  
  # Summarise by regulatory group and add a "NoInfo" row if needed
  summary_table <- temp %>% 
    group_by(Reg) %>% 
    summarise(n = n(), .groups = "drop") %>% 
    { 
      sum_n <- sum(.$n)
      bind_rows(., tibble(Reg = "NoInfo", n = total_genes - sum_n))
    } %>% 
    dplyr::mutate(PropSig = n / total_genes) %>% 
    dplyr::mutate(Reg = factor(Reg, levels = c("NoInfo", "MoreRep", "Neutral", "MoreAct")))
  
  return(summary_table)
}

get_MoreAct_genes <- function(gene_list,
                              data = Epith_Tumor_figR,
                              motif_filter = signature,
                              fill_colors = c("MoreAct" = "#1a9850",
                                              "MoreRep" = "#d73027",
                                              "Neutral" = "#FEAE2F",
                                              "NoInfo" = "grey80")) {
  
  total_genes <- length(gene_list)
  
  # Filter and process the data
  temp <- data %>% 
    dplyr::filter(Motif %in% motif_filter) %>% 
    dplyr::filter(DORC %in% gene_list) %>% 
    dplyr::mutate(Reg = case_when(
      Score > 0 ~ "Act",
      TRUE ~ "Rep"
    )) %>%
    group_by(DORC, Reg) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = "Reg", values_from = "n") %>%
    replace_na(list(Act = 0, Rep = 0)) %>%
    dplyr::mutate(PropAct = Act / (Act + Rep)) %>% 
    dplyr::mutate(Reg = case_when(
      is.na(PropAct)       ~ "NoInfo",
      PropAct == 0.5       ~ "Neutral",
      PropAct > 0.5        ~ "MoreAct",
      PropAct < 0.5        ~ "MoreRep"
    ))
  
  MoreAct_genes <- temp %>% filter(Reg == "MoreAct") %>% pull(DORC)
  
  return(MoreAct_genes)
}

gene_sets = list(cytokine_apoptosis = cytokine_apoptosis, stress_associated = stress_associated, Interferon_signaling = Interferon_signaling, differentiated = differentiated, EMT_associated = EMT_associated, AntigenPresentation = AntigenPresentation, Proteasomal_degradation = Proteasomal_degradation, TCA_cycle = TCA_cycle, Quiescence_associated = Quiescence_associated, Proliferative_DNA_repair = Proliferative_DNA_repair, RNA_processing = RNA_processing, PI3K = PI3K, Notch = Notch, commongenes = commongenes)

name_map = data.frame(gene_set = names(gene_sets), name = c("Cytokine & Apoptosis", "Stress Associated", "Interferon Signaling", "Glycosylation", "EMT Associated", "Antigen Presentation", "Proteasomal Degradation", "TCA cycle", "Quiescence Associated", "Proliferative DNA repair", "RNA processing", "PI3K", "Notch", "NACT-Resistant overlap")) %>% mutate(size = sapply(gene_sets, length)) %>% mutate(name = paste0(name, " [", size, "]"))

summaries = lapply(gene_sets, get_gene_summary) %>% bind_rows(.id = "gene_set") %>% mutate(gene_set = factor(gene_set, levels = rev(names(gene_sets)))) %>% left_join(name_map)
plotdata = summaries %>% mutate(name = factor(name, levels = summaries %>% filter(Reg == "MoreAct") %>% arrange(PropSig) %>% pull(name)))

pdf(file=paste0('Run10_pathway_specific_target_proportions_1.pdf'), height = 4.5, width = 6)
ggplot(plotdata, aes(x = n, y = name, fill = Reg)) +
    geom_bar(stat = "identity", position="fill", color = "black") +
	geom_text(plotdata %>% filter(Reg == "MoreAct") %>% mutate(n = as.character(n)), mapping = aes(label = n, x = PropSig), color = "white", nudge_x = -0.02, hjust = 1) +
	scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c("MoreAct" = "#1a9850", "MoreRep" = "#d73027", "Neutral" = "#FEAE2F", "NoInfo" = "grey80")) +
    theme_cowplot() +
    theme(legend.position = "none", axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank())
dev.off()

write.table(sapply(MoreAct_genes, function(x){return(paste(x, collapse = " "))}) %>% data.frame() %>% rownames_to_column("gene_set"), "Fig6_green_genes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

stress_associated_targets <- get_MoreAct_genes(stress_associated)
# ATF3 BCL6 CDKN1A CEBPD CREB5 DUSP1 EGR1 FOS GADD45B HBEGF HLA-G JUN MCL1 MYC NR4A1 PLK3 SNAI2 SOCS3 TNF

EMT_associated_targets <- get_MoreAct_genes(EMT_associated)
# BCAR1 BIRC2 COL1A2 DST E2F3 EPHB2 FADD ITCH ITGA2 ITGA6 MMP1 MMP12 MX1 PML SFN SMAD3

Quiescence_associated_targets <- get_MoreAct_genes(Quiescence_associated)
# ADAMTS6 AQP1 ARHGEF26 B2M BMPR1B CALD1 CD99 CDC42EP4 CXCL14 DTNA EEPD1 FAM107A GABARAPL2 GATM GLIS3 GLUL HSPB1 KAT2A MIR99AHG NDP NEAT1 NFIA NMB PDLIM3 PEA15 PIK3C2A PLPP3 PNRC1 RFX4 RGMA RSRP1 SARAF SLC6A9 SOX2 WASL

##################################### Figure 7 (PDX) #################################

filtered_PDX <- readRDS(paste0(runid, "_epithelial_PDX.rds"))

filtered_PDX <- AddModuleScore(
  object = filtered_PDX,
  features = list(patient_EPS_factors),
  assay = "remap",
  name = 'EPSv2_factors',
  ctrl = 10,
  seed = 123
)

pdf(file=paste0(runid, '_harmony_EPSv2_FeaturePlot_PDX.pdf'), height = 8, width = 8)
FeaturePlot(object = filtered_PDX, features = "EPSv2_factors1", min.cutoff = "q5", max.cutoff = "q95", reduction = "hmrnaumap", pt.size = 0.5, order = TRUE)
dev.off()

hm_PDX_Epi_2 <- filtered_PDX
scores <- hm_PDX_Epi_2@meta.data$EPSv2_factors1
q_low  <- quantile(scores, .25)
q_high <- quantile(scores, .75)

hm_PDX_Epi_2$quartile <- ifelse(scores <= q_low,  "Low", ifelse(scores >= q_high, "High", "Mid"))
hm_PDX_Epi_2$quartile <- factor(hm_PDX_Epi_2$quartile, levels = c("High", "Low", "Mid"))

p1 <- DimPlot(
    object    = hm_PDX_Epi_2,
    reduction = "hmrnaumap",
    group.by  = "quartile",
    pt.size   = 0.5,
    cols      = c(
      Low  = "#fdb515",
      Mid  = "grey80",
      High = "blue"
    )
  ) + 
ggtitle("Persister High & Low") +
theme_cowplot() +
theme(legend.title=element_blank())
pdf(file='Run10_harmony_EPSv2_High_Low_PDX.pdf', height = 8, width = 8)
ggrastr::rasterise(p1, dpi = 400)
dev.off()

celltype.table1 <- data.frame(celltype = hm_PDX_Epi_2$quartile, sample = hm_PDX_Epi_2$sample) %>% table() %>% data.frame() %>%
  dplyr::mutate(celltype = factor(celltype)) %>%
  dplyr::mutate(sample = factor(sample, levels = c("PH27_C", "PH27_T", "PH235_C", "PH235_T", "PH626_C", "PH626_T"))) %>%
  group_by(sample) %>%
  dplyr::mutate(freq = Freq / sum(Freq))

celltype.table1$celltype <- factor(celltype.table1$celltype, levels = c("Low", "Mid", "High"))

df2 <- celltype.table1 %>%
  separate(sample, into = c("PDX","State"), sep = "_") %>%
  dplyr::mutate(
    State = recode(State, C = "Naive", T = "Residual")
  ) %>%
  dplyr::select(PDX, State, celltype, freq)

df2$PDX <- factor(df2$PDX, levels = c("PH27","PH235","PH626"))

is_alluvia_form(df2, axes = 1:2, silent = TRUE)
p4 <- ggplot(df2,
       aes(x        = State,
           stratum  = celltype,
           alluvium = celltype,
           y        = freq,
           fill     = celltype,
           label    = celltype)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum(alpha = 0.8) +
  facet_wrap(~ PDX, scales = "free_x") +
  scale_fill_manual(values = c(Low  = "#feb500",
                               Mid  = "grey80",
                               High = "blue")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "High vs Low Persister Cells",
    x     = NULL,
    y     = "Proportion"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

ggsave("Run10_EPSv2_High_Low_Sankey_PDX.pdf", p4, width  = 6, height = 6)

############################# PDX violin plots #################################

table1 <- filtered_PDX@meta.data %>% dplyr::select(sample, group, patient, epicluster, EPSv2_factors1) %>%
	filter(EPSv2_factors1 != Inf & EPSv2_factors1 != -Inf) %>%
	mutate(patient = factor(patient, levels = c("PH27", "PH235", "PH626")))

pdf(file=paste0(runid, "_EPSv2_VlnPlots_PDX.pdf"), width=10, height=6)
p1 <- ggplot(table1, aes(x=patient, y=EPSv2_factors1)) + 
    geom_violin(colour = "black", fill = "grey80") +
	stat_compare_means(aes(label = ..p.signif..), comparisons = list(c("PH27", "PH235"), c("PH27", "PH626"), c("PH235", "PH626"))) +
	# ggrastr::rasterise(geom_jitter(size=0.001, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
	stat_summary(fun=mean, geom="point", size=1, color="black") +
    # coord_cartesian(ylim = c(-4, 12)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- ggplot(table1, aes(x=group, y=EPSv2_factors1)) + 
    geom_violin(colour = "black", fill = "grey80") +
	stat_compare_means(aes(label = ..p.signif..)) +
	# ggrastr::rasterise(geom_jitter(size=0.1, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
	stat_summary(fun=mean, geom="point", size=1, color="black") +
    # coord_cartesian(ylim = c(-4, 11)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 <- ggplot(table1, aes(x=epicluster, y=EPSv2_factors1)) + 
    geom_violin(colour = "black", fill = "grey80") +
	# ggrastr::rasterise(geom_jitter(size=0.001, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
	stat_summary(fun=mean, geom="point", size=1, color="black") +
    # coord_cartesian(ylim = c(-5, 11)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1 + p2 + p3 + plot_layout(widths = c(2, 2, 3))
dev.off()

pdf(file=paste0(runid, "_EPSv2_VlnPlots_PDX_Untreated.pdf"), width=5, height=6)
ggplot(table1 %>% filter(group == "PDX_C"), aes(x=patient, y=EPSv2_factors1)) + 
    geom_violin(colour = "black", fill = "grey80") +
	stat_compare_means(aes(label = ..p.signif..), comparisons = list(c("PH27", "PH235"), c("PH27", "PH626"), c("PH235", "PH626"))) +
	# ggrastr::rasterise(geom_jitter(size=0.001, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
	stat_summary(fun=mean, geom="point", size=1, color="black") +
    # coord_cartesian(ylim = c(-4, 12)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

###########################################################################
############################### Primed genes ##############################
###########################################################################

filtered_PDX$persist <- case_when(filtered_PDX$EPSv2_factors1 > quantile(filtered_PDX$EPSv2_factors1)[4] ~ "High", filtered_PDX$EPSv2_factors1 < quantile(filtered_PDX$EPSv2_factors1)[2] ~ "Low", TRUE ~ "Mid")

filtered_PDX$sample_persist <- factor(paste0(filtered_PDX$sample, "_", filtered_PDX$persist), levels = c("PH27_C_Low", "PH27_C_High", "PH235_C_Low", "PH235_C_High", "PH626_C_Low", "PH626_C_High", "PH27_T_Low", "PH27_T_High", "PH235_T_Low", "PH235_T_High", "PH626_T_Low", "PH626_T_High"))
filtered_PDX_extremes <- subset(filtered_PDX, subset = (persist %in% c("High", "Low")))
filtered_PDX_extremes$group_persist <- factor(paste0(filtered_PDX_extremes$group, "_", filtered_PDX_extremes$persist), levels = c("PDX_C_Low", "PDX_C_High", "PDX_T_Low", "PDX_T_High"))

DefaultAssay(filtered_PDX_extremes) <- "mATAC"
Idents(filtered_PDX_extremes) <- "persist"
da_peaks1 <- FindMarkers(
  object = filtered_PDX_extremes,
  ident.1 = "High",
  ident.2 = "Low",
  test.use = 'wilcox'
) %>% filter(p_val_adj < 0.05)
da_peaks <- da_peaks1 %>% filter(avg_log2FC > 0) %>% rownames()
dd_peaks <- da_peaks1 %>% filter(avg_log2FC < 0) %>% rownames()

write.table(data.frame(peak = da_peaks) %>% separate(peak, c("chr", "start", "end"), sep = "-") %>% arrange(chr, as.numeric(start)), "high_peaks.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(data.frame(peak = dd_peaks) %>% separate(peak, c("chr", "start", "end"), sep = "-") %>% arrange(chr, as.numeric(start)), "low_peaks.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# wc -l low_peaks.bed # 21622 low_peaks.bed
# wc -l high_peaks.bed # 12495 high_peaks.bed
# wc -l patient_EPS_peaks.bed # 1762 patient_EPS_peaks.bed
# bedtools intersect -a high_peaks.bed -b patient_EPS_peaks.bed -u | wc -l # 644
# bedtools intersect -a low_peaks.bed -b patient_EPS_peaks.bed -u | wc -l # 22

######################### BigWigs - DeepTools heatmap #############################

DefaultAssay(filtered_PDX_extremes) <- "mATAC"
Idents(filtered_PDX_extremes) <- "group_persist"
ExportGroupBW(filtered_PDX_extremes,
    assay = 'mATAC',
    group.by = 'group_persist',
    idents = NULL,
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    cutoff = NULL,
    genome = genome,
    outdir = "Run10",
    verbose=TRUE
)

# cd Run10
# module load python/2.7.18
# computeMatrix reference-point -S PDX_C_Low-TileSize-100-normMethod-rc.bw PDX_T_Low-TileSize-100-normMethod-rc.bw PDX_C_High-TileSize-100-normMethod-rc.bw PDX_T_High-TileSize-100-normMethod-rc.bw \
                              # -R low_peaks.bed \
                              # --referencePoint center \
                              # -b 1000 -a 1000 -p 16 \
                              # --missingDataAsZero \
                              # --outFileNameMatrix low_peaks.clust \
                              # --outFileName low_peaks.clust.tab.gz \
                              # --sortRegions descend \
                              # --sortUsing mean &
# computeMatrix reference-point -S PDX_C_Low-TileSize-100-normMethod-rc.bw PDX_T_Low-TileSize-100-normMethod-rc.bw PDX_C_High-TileSize-100-normMethod-rc.bw PDX_T_High-TileSize-100-normMethod-rc.bw \
                              # -R high_peaks.bed \
                              # --referencePoint center \
                              # -b 1000 -a 1000 -p 16 \
                              # --missingDataAsZero \
                              # --outFileNameMatrix high_peaks.clust \
                              # --outFileName high_peaks.clust.tab.gz \
                              # --sortRegions descend \
                              # --sortUsing mean &

# plotHeatmap -m low_peaks.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out low_peaks.pdf --samplesLabel C_Low T_Low C_High T_High --heatmapHeight 8 &
# plotHeatmap -m low_peaks.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out low_peaks.png --samplesLabel C_Low T_Low C_High T_High &
# plotHeatmap -m high_peaks.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out high_peaks.pdf --samplesLabel C_Low T_Low C_High T_High --heatmapHeight 8 &
# plotHeatmap -m high_peaks.clust.tab.gz --colorList 'white,blue' --sortRegions keep -out high_peaks.png --samplesLabel C_Low T_Low C_High T_High &

############### Identifying primed genes #############

# Get targets
DefaultAssay(filtered_PDX_extremes) <- "mATAC"
genes <- ClosestFeature(
  object = filtered_PDX_extremes,
  regions = da_peaks
  # regions = PDX_EPS_peaks
) %>% filter(gene_biotype == "protein_coding" & distance < 150000)
targets_by_proximity <- genes$gene_name %>% unique() %>% sort()

Epith_Tumor_figR <- read.table("RunPDX_figR_motif_relaxed_sub.tsv", header = TRUE) %>% set_colnames(c("DORC", "Motif", "Score"))
targets_by_FigR <- Epith_Tumor_figR %>% filter(Motif %in% patient_EPS_factors) %>% pull(DORC) %>% unique() %>% sort()

targets <- base::intersect(targets_by_proximity, targets_by_FigR)
# targets <- targets_by_proximity

objs <- lapply(c("PH27", "PH235", "PH626"), function(x){subset(filtered_PDX_extremes, subset = patient == x)})

get_change_tbl <- function(obj){
	DefaultAssay(obj) <- "RNA"
	Idents(obj) <- "group_persist"
	degs1 <- FindMarkers(
	  object = obj,
	  ident.1 = "PDX_C_High",
	  ident.2 = "PDX_C_Low"
	) %>% filter(p_val_adj < 0.05)
	CHLdf <- degs1 %>% rownames_to_column("gene") %>% 
		mutate(CHL = case_when(avg_log2FC > 0 ~ "UP", TRUE ~ "DOWN")) %>% 
		dplyr::select(gene, CHL)

	DefaultAssay(obj) <- "RNA"
	Idents(obj) <- "group_persist"
	degs1 <- FindMarkers(
	  object = obj,
	  ident.1 = "PDX_T_High",
	  ident.2 = "PDX_T_Low"
	) %>% filter(p_val_adj < 0.05)
	THLdf <- degs1 %>% rownames_to_column("gene") %>% 
		mutate(THL = case_when(avg_log2FC > 0 ~ "UP", TRUE ~ "DOWN")) %>% 
		dplyr::select(gene, THL)

	DefaultAssay(obj) <- "RNA"
	Idents(obj) <- "group_persist"
	degs1 <- FindMarkers(
	  object = obj,
	  ident.1 = "PDX_T_Low",
	  ident.2 = "PDX_C_Low"
	) %>% filter(p_val_adj < 0.05)
	CTLdf <- degs1 %>% rownames_to_column("gene") %>% 
		mutate(CTL = case_when(avg_log2FC > 0 ~ "UP", TRUE ~ "DOWN")) %>% 
		dplyr::select(gene, CTL)

	DefaultAssay(obj) <- "RNA"
	Idents(obj) <- "group_persist"
	degs1 <- FindMarkers(
	  object = obj,
	  ident.1 = "PDX_T_High",
	  ident.2 = "PDX_C_High"
	) %>% filter(p_val_adj < 0.05)
	CTHdf <- degs1 %>% rownames_to_column("gene") %>% 
		mutate(CTH = case_when(avg_log2FC > 0 ~ "UP", TRUE ~ "DOWN")) %>% 
		dplyr::select(gene, CTH)

	target_changes <- CHLdf %>% full_join(THLdf) %>% full_join(CTLdf) %>% full_join(CTHdf) %>% filter(gene %in% targets)
	return(target_changes)
}

obj <- filtered_PDX_extremes
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "group_persist"
degs1 <- FindMarkers(
  object = obj,
  ident.1 = "PDX_C_High",
  ident.2 = "PDX_C_Low"
) %>% filter(p_val_adj < 0.05)
CHLdf <- degs1 %>% rownames_to_column("gene")

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "group_persist"
degs1 <- FindMarkers(
  object = obj,
  ident.1 = "PDX_T_High",
  ident.2 = "PDX_C_High"
) %>% filter(p_val_adj < 0.05)
CTHdf <- degs1 %>% rownames_to_column("gene")
CTHfc <- FoldChange(obj, ident.1 = "PDX_T_High", ident.2 = "PDX_C_High") %>% rownames_to_column("gene")

change_tbls_targets <- lapply(c(filtered_PDX_extremes, objs), get_change_tbl)
change_tbls <- change_tbls_targets
temp1 <- change_tbls[[1]] %>% left_join(change_tbls[[2]] %>% mutate(PH27 = 1)) %>% left_join(change_tbls[[3]] %>% mutate(PH235 = 1)) %>% left_join(change_tbls[[4]] %>% mutate(PH626 = 1))

# write.table(temp1 %>% replace_na(list(CHL = "NoChange", THL = "NoChange", CTL = "NoChange", CTH = "NoChange", PH27 = 0, PH235 = 0, PH626 = 0)), "Run10_PDX_primed_genes_analysis.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

temp1 %<>% replace_na(list(CHL = "NoChange", THL = "NoChange", CTL = "NoChange", CTH = "NoChange", PH27 = 0, PH235 = 0, PH626 = 0))

# Using both targets_by_FigR and targets_by_proximity filter in the beginning
temp1 %>% filter(CTH != "NoChange") %>% nrow() # 227
temp1 %>% filter(CTL != "NoChange") %>% nrow() # 370
temp1 %>% filter(CTH != "NoChange" & CTL != "NoChange") %>% nrow() # 110
temp1 %>% filter(CTH != "NoChange" & CTL == "NoChange") %>% nrow() # 117
temp1 %>% filter(CTH == "NoChange" & CTL != "NoChange") %>% nrow() # 260
temp1 %>% filter(CTH == "UP" & CTL == "NoChange") %>% nrow() # 92
temp1 %>% filter(CTH == "DOWN" & CTL == "NoChange") %>% nrow() # 25
temp1 %>% filter(CTH == "UP" & CTL == "NoChange" & CHL != "NoChange") %>% nrow() # 48
temp1 %>% filter(CTH == "UP" & CTL == "NoChange" & CHL == "NoChange") %>% nrow() # 44

genelist44 <- temp1 %>% filter(CTH == "UP" & CTL == "NoChange" & CHL == "NoChange") %>% pull(gene)

volcdata <- CTHdf %>% filter(gene %in% (temp1 %>% filter(CTH != "NoChange" & CTL == "NoChange") %>% pull(gene))) %>% left_join(CHLdf %>% dplyr::select(gene, avg_log2FC) %>% rename(CHL = avg_log2FC)) %>% replace_na(list(CHL = 0))
volcdata_hl <- temp1 %>% filter(CTH == "UP" & CTL == "NoChange" & CHL == "NoChange") %>% pull(gene)

volcdata <- volcdata %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))

pdf("Run10_primed_genes_volcano.pdf", width=5, height=5)
ggplot(volcdata, aes(x = avg_log2FC, y = neg_log10_adj_pval, color = CHL)) +
	# geom_point(size = 2, stroke=0.1, shape = 16) +
	geom_point(data = volcdata %>% filter(!(gene %in% volcdata_hl)), size = 2, stroke=0.1, shape = 16) +
	geom_point(data = volcdata %>% filter(gene %in% volcdata_hl), mapping = aes(fill = CHL), size = 2, stroke=0.5, shape = 21, color = "black") +
	# geom_label_repel(data = volcdata %>% filter(gene == "GIPC1"), aes(label = gene), size = 3, max.overlaps = Inf, hjust = 0, nudge_x = -1, direction = "y", segment.size = 0.5, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 3, fontface = "italic") +
	scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-5,5),breaks=scales::breaks_pretty(n=3)) + # ,oob = scales::squish
	scale_fill_gradientn(colours = jdb_palette("solar_extra"),limits=c(-5,5),breaks=scales::breaks_pretty(n=3)) + # ,oob = scales::squish
	# coord_cartesian(xlim = c(-4, 4)) +
	labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
	theme_cowplot()
dev.off()

DefaultAssay(filtered) <- "RNA"
Idents(filtered) <- "subgroup2"
degs1 <- FindMarkers(
  object = filtered,
  ident.1 = "NACT",
  ident.2 = "Naive"
) %>% filter(p_val_adj < 0.05)
Naive_vs_NACT <- degs1 %>% rownames_to_column("gene")
up_in_NACT_vs_Naive <- degs1 %>% filter(avg_log2FC > 0) %>% rownames() %>% unique() %>% sort()
down_in_NACT_vs_Naive <- degs1 %>% filter(avg_log2FC < 0) %>% rownames() %>% unique() %>% sort()

genelist33 <- setdiff(temp1 %>% filter(CTH == "UP" & CTL == "NoChange" & CHL == "NoChange") %>% pull(gene), down_in_NACT_vs_Naive)
# ELMO1 GIPC1 STK35 KIF16B SLC26A9 ELK3 NHSL2 LIF MEF2A POU2F2 NOS3 RHOB CYHR1 C4orf19 PLSCR2 RRM2B F11R CMPK1 INSL6 ATP6V0A1 GPNMB OGDH SLC2A12 NCK2 RALB PLEKHG5 SDHC TMC7 PTPRR ARL14 TMEM40 CDHR2 S1PR1

write.table(data.frame(gene = genelist33), "Run10_primed_genes.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

temp1 %>% filter(CTH == "UP") %>% nrow() # 177
temp1 %>% filter(CTH == "UP" & CTL == "NoChange" & CHL == "NoChange") %>% nrow() # 44
setdiff(temp1 %>% filter(CTH == "UP" & CTL == "NoChange" & CHL == "NoChange") %>% pull(gene), down_in_NACT_vs_Naive) %>% length() # 33

primed_pie <- data.frame(n = c(33, 11, 133), label = as.character(c(33, 11, 133)))

primed_pie <- primed_pie %>% 
	arrange(desc(label)) %>%
	mutate(prop = n / sum(primed_pie$n) *100) %>%
	mutate(ypos = cumsum(prop)- 0.5*prop )

pdf("Run10_primed_genes_pie.pdf", width=5, height=5)
ggplot(primed_pie, aes(x="", y=prop, fill=label)) +
	geom_bar(stat="identity", width=1, color="white") +
	coord_polar("y", start=0) +
	theme_void() + 
	theme(legend.position="none") +
	geom_text(aes(y = ypos, label = label), color = "white", size=6) +
	scale_fill_brewer(palette="Set1")
dev.off()
