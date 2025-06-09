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

sample_levels <- c("FT243", "FT244", "FT245", "FT281", "FT283", "PH1199", "PH1230", "PH1239", "PH1243", "PH1302", "PH1303", "PH1214", "PH1238", "PH1252")
subgroup2_levels <- c("FT", "Naive", "NACT")
subgroup3_levels <- c("BRCAwt", "BRCAmt", "Sensitive", "Resistant", "NACT")
infotable <- data.frame(
    sample = sample_levels,
	subgroup2 = factor(c("FT", "FT", "FT", "FT", "FT", "Naive", "Naive", "Naive", "Naive", "Naive", "Naive", "NACT", "NACT", "NACT"), levels = subgroup2_levels),
    subgroup3 = factor(c("BRCAwt", "BRCAwt", "BRCAwt", "BRCAmt", "BRCAmt", "Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Resistant", "NACT", "NACT", "NACT"), levels = subgroup3_levels)
)

Double_sig <- c("AFF4", "APC", "AR", "ARID1A", "ARNTL", "ATF3", "ATOH8", "ATXN7L3", "BICRA", "BRCA1", "BRD2", "BRD4", "BRD9", "CARM1", "CASZ1", "CBX1", "CBX3", "CDKN1B", "CEBPB", "CEBPG", "CHD2", "CHD4", "CREB5", "CREBBP", "DAXX", "DEK", "DPF2", "E2F1", "E2F7", "EHF", "ELF3", "ELL", "ELL2", "EP300", "EPAS1", "ESR1", "ESR2", "ESRRA", "ETV1", "FOS", "FOSB", "FOSL1", "FOSL2", "FOXA2", "FOXJ2", "FOXM1", "GATA6", "GATAD1", "GLI2", "GPS2", "GTF2A2", "HDAC2", "HIF1A", "HLF", "HMGB2", "ICE2", "JUN", "JUNB", "JUND", "KDM5B", "KLF3", "KLF4", "KLF5", "KLF6", "KMT2C", "KMT2D", "MAFB", "MAFF", "MAML1", "MAX", "MBD2", "MECOM", "MED1", "MED25", "MED26", "MGA", "MLXIP", "MYC", "NCAPH2", "NCOA2", "NCOA3", "NELFE", "NFE2L1", "NFE2L2", "NFIA", "NFIC", "NFIL3", "NFIX", "NFKB1", "NIPBL", "NR1H2", "NR2F2", "NR3C1", "NRF1", "NSD2", "ONECUT2", "PAF1", "PARP1", "PAX6", "PAX8", "PGR", "PHIP", "PPARG", "PRMT5", "RAD21", "RBAK", "RBPJ", "RCOR1", "RELA", "RFX1", "RFX2", "RFX5", "RUVBL2", "RXRA", "SFMBT1", "SIN3A", "SMAD3", "SMAD4", "SMARCA2", "SMARCA4", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCE1", "SMC1A", "SMC3", "SOX2", "SOX4", "SOX9", "SP1", "SPDEF", "SRC", "SREBF2", "SS18", "STAG1", "STAT3", "SUPT5H", "SUPT6H", "TCF7L2", "TEAD1", "TEAD2", "TEAD4", "TFAP2C", "TFAP4", "TFEB", "TP53", "TP63", "TSC22D4", "TSHZ1", "UBN1", "USF1", "USF2", "VDR", "WT1", "XBP1", "YY1AP1", "ZBTB24", "ZC3H8", "ZEB1", "ZHX1", "ZIC2", "ZMYM2", "ZMYND8", "ZNF114", "ZNF19", "ZNF26", "ZNF280D", "ZNF398", "ZNF548", "ZNF554", "ZNF563", "ZNF624", "ZNF667", "ZNF677", "ZNF750", "ZNF791", "ZNF90", "ZNF92")

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

sample_levels <- c("FT243", "FT244", "FT245", "FT281", "FT283", "PH1199", "PH1230", "PH1239", "PH1243", "PH1302", "PH1303", "PH1214", "PH1238", "PH1252")
filtered$sample <- factor(filtered$sample, levels = sample_levels)

subgroup2_levels <- c("FT", "Naive", "NACT")
subgroup3_levels <- c("BRCAwt", "BRCAmt", "Sensitive", "Resistant", "NACT")
infotable <- data.frame(
    sample = sample_levels,
	subgroup2 = factor(c("FT", "FT", "FT", "FT", "FT", "Naive", "Naive", "Naive", "Naive", "Naive", "Naive", "NACT", "NACT", "NACT"), levels = subgroup2_levels),
    subgroup3 = factor(c("BRCAwt", "BRCAwt", "BRCAwt", "BRCAmt", "BRCAmt", "Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Resistant", "NACT", "NACT", "NACT"), levels = subgroup3_levels)
)
filtered$subgroup2 <- data.frame(sample = filtered$sample) %>% left_join(infotable) %>% pull(subgroup2)
filtered$subgroup3 <- data.frame(sample = filtered$sample) %>% left_join(infotable) %>% pull(subgroup3)

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

########### cell tagging ###################
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

# Markers
Rscript findAllMarkers_generic.R "/working" "Run10_celltyped.rds" "customtype2" "RNA" "Run10_celltyped" 1> 1 2> 2 &

# ChromVAR (Saved in same object)
Rscript ChromVAR_generic.R /working "Run10_celltyped.rds" 1> 1 2> 2 &
Rscript ChromVAR_generic.R /working "Run10_epithelial.rds" 1> 11 2> 22 &

# DEGs
Rscript degs_generic.R /working Run10_celltyped.rds Tumor FT customtype2 group Run10 1> 1 2> 2 &
Rscript degs_generic.R /working Run10_celltyped.rds BRCAmt BRCAwt customtype2 subgroup3 Run10 1> 1 2> 2 &
Rscript degs_generic.R /working Run10_celltyped.rds NACT Naive customtype2 subgroup2 Run10 1> 1 2> 2 &
Rscript degs_generic.R /working Run10_celltyped.rds Naive FT customtype2 subgroup2 Run10 1> 1 2> 2 &
Rscript degs_generic.R /working Run10_celltyped.rds NACT FT customtype2 subgroup2 Run10 1> 1 2> 2 &

Rscript degs_generic.R /working Run10_celltyped.rds Resistant Sensitive customtype2 subgroup3 Run10 1> 1 2> 2 &
Rscript degs_generic.R /working Run10_celltyped.rds NACT Sensitive customtype2 subgroup3 Run10 1> 1 2> 2 &
Rscript degs_generic.R /working Run10_celltyped.rds NACT Resistant customtype2 subgroup3 Run10 1> 1 2> 2 &

# Rscript degs_generic.R /working Run10_celltyped.rds Sensitive FT customtype2 subgroup3 Run10 1> 1 2> 2 &
# Rscript degs_generic.R /working Run10_celltyped.rds Resistant FT customtype2 subgroup3 Run10 1> 1 2> 2 &

# Cicero
Rscript Cicero_generic.R /working Run10_T_cell_CD8.rds seurat_clusters WNNUMAP Run10_T_cell_CD8 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_NK_cell.rds seurat_clusters WNNUMAP Run10_NK_cell 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Monocyte.rds seurat_clusters WNNUMAP Run10_Monocyte 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Macrophage.rds seurat_clusters WNNUMAP Run10_Macrophage 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Endothelial_cells.rds seurat_clusters WNNUMAP Run10_Endothelial_cells 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Fibroblasts.rds seurat_clusters WNNUMAP Run10_Fibroblasts 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Epith_Tumor.rds seurat_clusters WNNUMAP Run10_Epith_Tumor 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Epith_FT.rds seurat_clusters WNNUMAP Run10_Epith_FT 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Epith.rds seurat_clusters WNNUMAP Run10_Epith 1> 1 2> 2 &

Rscript Cicero_generic.R /working Run10_Ciliated.rds seurat_clusters WNNUMAP Run10_Ciliated 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Secretory.rds seurat_clusters WNNUMAP Run10_Secretory 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Cluster1.rds seurat_clusters WNNUMAP Run10_Cluster1 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Cluster2.rds seurat_clusters WNNUMAP Run10_Cluster2 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Cluster3.rds seurat_clusters WNNUMAP Run10_Cluster3 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Cycling.rds seurat_clusters WNNUMAP Run10_Cycling 1> 1 2> 2 &

Rscript Cicero_generic.R /working Run10_Ciliated_FT.rds seurat_clusters WNNUMAP Run10_Ciliated_FT 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Secretory_FT.rds seurat_clusters WNNUMAP Run10_Secretory_FT 1> 1 2> 2 &

Rscript Cicero_generic.R /working Run10_Ciliated_Tumor.rds seurat_clusters WNNUMAP Run10_Ciliated_Tumor 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Secretory_Tumor.rds seurat_clusters WNNUMAP Run10_Secretory_Tumor 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Cluster1_Tumor.rds seurat_clusters WNNUMAP Run10_Cluster1_Tumor 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Cluster2_Tumor.rds seurat_clusters WNNUMAP Run10_Cluster2_Tumor 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Cluster3_Tumor.rds seurat_clusters WNNUMAP Run10_Cluster3_Tumor 1> 1 2> 2 &
Rscript Cicero_generic.R /working Run10_Cycling_Tumor.rds seurat_clusters WNNUMAP Run10_Cycling_Tumor 1> 1 2> 2 &

# FigR
Rscript FigR_generic.R /working Run10_epithelial_raw.rds Run10_Epith mATAC hg38 1
Rscript FigR_generic.R /working Run10_epithelial_raw.rds Run10_Epith mATAC hg38 2
Rscript FigR_generic.R /working Run10_epithelial_raw.rds Run10_Epith mATAC hg38 3

Rscript FigR_generic.R /working Run10_Epith_FT.rds Run10_Epith_FT mATAC hg38 1
Rscript FigR_generic.R /working Run10_Epith_FT.rds Run10_Epith_FT mATAC hg38 2
Rscript FigR_generic.R /working Run10_Epith_FT.rds Run10_Epith_FT mATAC hg38 3

Rscript FigR_generic.R /working Run10_Epith_NACT.rds Run10_Epith_NACT mATAC hg38 1
Rscript FigR_generic.R /working Run10_Epith_NACT.rds Run10_Epith_NACT mATAC hg38 2
Rscript FigR_generic.R /working Run10_Epith_NACT.rds Run10_Epith_NACT mATAC hg38 3

Rscript FigR_generic.R /working Run10_Epith_Naive.rds Run10_Epith_Naive mATAC hg38 1
Rscript FigR_generic.R /working Run10_Epith_Naive.rds Run10_Epith_Naive mATAC hg38 2
Rscript FigR_generic.R /working Run10_Epith_Naive.rds Run10_Epith_Naive mATAC hg38 3

Rscript FigR_generic.R /working Run10_Epith_PDX.rds Run10_Epith_PDX mATAC hg38 1
Rscript FigR_generic.R /working Run10_Epith_PDX.rds Run10_Epith_PDX mATAC hg38 2
Rscript FigR_generic.R /working Run10_Epith_PDX.rds Run10_Epith_PDX mATAC hg38 3

Rscript FigR_generic.R /working Run10_Epith_Tumor.rds Run10_Epith_Tumor mATAC hg38 1
Rscript FigR_generic.R /working Run10_Epith_Tumor.rds Run10_Epith_Tumor mATAC hg38 2
Rscript FigR_generic.R /working Run10_Epith_Tumor.rds Run10_Epith_Tumor mATAC hg38 3

############################ Plot - DEGS - DAMS - DARS #############################

degsr <- read.table(file = paste0("Run10", "_degsr_", "Persister", "_", "NonPersister", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05)
degsr_NACT_Naive <- read.table(file = paste0("Run10", "_degsr_", "NACT", "_", "Naive", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells") %>% dplyr::select(celltype, gene, avg_log2FC)
degsr_Resistant_Sensitive <- read.table(file = paste0("Run10", "_degsr_", "Resistant", "_", "Sensitive", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells") %>% dplyr::select(celltype, gene, avg_log2FC)
UPR <- c("ASNS", "ATF3", "CALR", "CEBPB", "DNAJB11", "DNAJC3", "ERN1", "GFPT1", "HSP90B1", "HSPA5", "KDELR3", "PDIA5", "TPP1", "ERO1A", "HSPA9", "IARS", "NOP56", "PSAT1", "TARS")
ATF4_targets <- c("MTHFD2", "PYCR1", "GARS", "ASNS", "FTH1", "RAB32", "PHGDH", "SLC1A5", "HMOX1")
degsr_NACT_Naive %>% filter(gene %in% UPR)
degsr_NACT_Naive %>% filter(gene %in% ATF4_targets)
degsr_Resistant_Sensitive %>% filter(gene %in% UPR)
degsr_Resistant_Sensitive %>% filter(gene %in% ATF4_targets)

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

unitlabels <- c("SOX17", "MECOM", "PAX8", "WT1", "FOXK2", "SOX2", "POU5F1", "NANOG", "ESR1", "PGR", "CD47", "MAML1")
plot_differential_expression("Run10", "FT", "Tumor", names(colors.customtype2), colors.customtype2, unitlabels)
# plot_differential_expression("Run10", "FT", "Naive", names(colors.customtype2), colors.customtype2, unitlabels)
unitlabels <- c("THRAP3", "ZNF462", "SOX17", "RPA2", "MPHOSPH8", "HCFC1R1", "HMGA1", "ZNF781", "ZFP57", "GRHL1", "WT1", "MECOM", "PAX8", "ZNF217", "ESR1", "NR2F6", "TGIF2", "MEIS1", "PBX1")
plot_differential_expression("Run10", "FT", "Naive", "Epithelial_cells", colors.customtype2["Epithelial_cells"], unitlabels)
plot_differential_expression("Run10", "BRCAwt", "BRCAmt", names(colors.customtype2), colors.customtype2, unitlabels)
plot_differential_expression("Run10", "Naive", "NACT", names(colors.customtype2), colors.customtype2, unitlabels)
# plot_differential_expression("Run10", "Naive", "NACT", "Epithelial_cells", colors.customtype2["Epithelial_cells"], "")
plot_differential_expression("Run10", "Sensitive", "Resistant", names(colors.customtype2), colors.customtype2, unitlabels)
# plot_differential_expression("Run10", "Sensitive", "Resistant", "Epithelial_cells", colors.customtype2["Epithelial_cells"], "")

colors.customtype2_alex <- c(`T_cell:CD8+` = "#A9E0AB", `T_cell:CD4+` = "#665299", `T_cell:Other` = "#31c4e9", `NK_cell` = "#E2E8A6", `B_cell` = "#D68645", `Monocyte` = "#B4AEE5", `Macrophage` = "#B81254", `DC` = "#67A5E5", `Neutrophils` = "#3B69B7", `Endothelial_cells` = "#ffd92f", `Fibroblasts` = "#B81223", `Epithelial_cells` = "#1E4A3A")
plot_differential_expression("Run10", "FT", "Naive", names(colors.customtype2), colors.customtype2_alex, c())

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

filtered2 <- subset(filtered, subset = nCount_mATAC < 5000)
Tumor2 <- subset(filtered2, subset = group == "Tumor")
differential_enrichment(Tumor2, "remap", "customtype2", "subgroup2", "Naive", "NACT", "Run10_filter2")
Naive2 <- subset(Tumor2, subset = subgroup2 == "Naive")
differential_enrichment(Naive2, "remap", "customtype2", "subgroup3", "Sensitive", "Resistant", "Run10_filter2")
filtered3 <- subset(filtered, subset = nCount_mATAC > 5000)
Tumor3 <- subset(filtered3, subset = group == "Tumor")
differential_enrichment(Tumor3, "remap", "customtype2", "subgroup2", "Naive", "NACT", "Run10_filter3")
Naive3 <- subset(Tumor3, subset = subgroup2 == "Naive")
differential_enrichment(Naive3, "remap", "customtype2", "subgroup3", "Sensitive", "Resistant", "Run10_filter3")


plot_differential_enrichment <- function(prefix, assay, leftgroup, rightgroup, celltypes, celltype_colors, unitlabels){
    deus <- read.table(file = paste0(prefix, "_differential_", assay, "_", rightgroup, "_", leftgroup, ".tsv"), sep = "\t", skip = 1) %>% set_colnames(c("celltype", "unit", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj"))
    deusr <- read.table(file = paste0(prefix, "_differentialr_", assay, "_", rightgroup, "_", leftgroup, ".tsv"), sep = "\t", skip = 1) %>% set_colnames(c("celltype", "unit", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj"))

    temp1 <- deus %>% 
        filter(p_val_adj < 0.05) %>%
        mutate(celltype = factor(celltype, levels = names(celltype_colors))) %>% 
        mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
        mutate(neg_log10_p_val_adj = -log10(p_val_adj)) %>%
        filter(celltype %in% celltypes)
    temp2 <- deusr %>% 
        filter(p_val_adj < 0.05) %>%
        mutate(celltype = factor(celltype, levels = names(celltype_colors))) %>% 
        mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
        mutate(neg_log10_p_val_adj = -log10(p_val_adj)) %>%
        filter(celltype %in% celltypes)
    temp2a <- temp2 %>% filter(p_val_adj != 1e-305)
	temp2b <- temp2 %>% filter(p_val_adj == 1e-305) %>% mutate(neg_log10_p_val_adj = jitter(neg_log10_p_val_adj, amount = 5))
	forlabels1 <- temp1 %>% mutate(unit = case_when((unit %in% unitlabels) ~ unit, TRUE ~ ""))
    forlabels2a <- temp2a %>% mutate(unit = case_when((unit %in% unitlabels) ~ unit, TRUE ~ ""))
	forlabels2b <- temp2b %>% mutate(unit = case_when((unit %in% unitlabels) ~ unit, TRUE ~ ""))
    p1 <- ggplot(temp1, aes(x = avg_diff, y = neg_log10_p_val_adj)) +
        geom_vline(xintercept = 0, linetype = "dashed") +
		geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
		geom_point(aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
        geom_label_repel(data = forlabels1 %>% filter(avg_diff > 0), aes(label = unit, segment.color=celltype), size = 3, max.overlaps = Inf, hjust = 0, nudge_x = 1, direction = "y", segment.size = 0.2, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 1.5) +
        geom_label_repel(data = forlabels1 %>% filter(avg_diff < 0), aes(label = unit, segment.color=celltype), size = 3, max.overlaps = Inf, hjust = 1, nudge_x = -1, direction = "y", segment.size = 0.2, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 1.5) +
        labs(x = "Average difference in assay score", y = "-log10(Adjusted P value)") +
        scale_color_manual(values = celltype_colors, aesthetics = c("color", "segment.color")) +
        coord_cartesian(ylim=c(0, 330)) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() +
        theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title = element_blank())
    p2 <- ggplot(mapping = aes(x = avg_diff, y = neg_log10_p_val_adj)) +
        geom_vline(xintercept = 0, linetype = "dashed") +
		geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
		geom_point(data = temp2a, aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
        geom_point(data = temp2b, aes(color = celltype), size = 1.5, stroke=0.1, shape = 16) +
        geom_label_repel(data = forlabels2a %>% filter(avg_diff > 0), aes(label = unit, segment.color=celltype), size = 3, max.overlaps = Inf, hjust = 0, nudge_x = 1, direction = "y", segment.size = 0.2, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 1.5) +
        geom_label_repel(data = forlabels2a %>% filter(avg_diff < 0), aes(label = unit, segment.color=celltype), size = 3, max.overlaps = Inf, hjust = 1, nudge_x = -1, direction = "y", segment.size = 0.2, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 1.5) +
        geom_label_repel(data = forlabels2b %>% filter(avg_diff > 0), aes(label = unit, segment.color=celltype), size = 3, max.overlaps = Inf, hjust = 0, nudge_x = 1, direction = "y", segment.size = 0.2, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 1.5) +
        geom_label_repel(data = forlabels2b %>% filter(avg_diff < 0), aes(label = unit, segment.color=celltype), size = 3, max.overlaps = Inf, hjust = 1, nudge_x = -1, direction = "y", segment.size = 0.2, show.legend = FALSE, arrow = arrow(length = unit(0.010, "npc")), force = 1.5) +
        labs(x = "Average difference in assay score", y = "-log10(Adjusted P value)") +
        scale_color_manual(values = celltype_colors, aesthetics = c("color", "segment.color")) +
        coord_cartesian(ylim=c(0, 330)) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() +
        theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.title = element_blank())
    pdf(file = paste0(prefix, "_differential_", assay, "_", rightgroup, "_", leftgroup, ".pdf"), width = 15, height = 6)
    plot(p1 | p2)
    dev.off()
}

unitlabels <- c("SOX17", "MECOM", "PAX8", "WT1", "FOXK2", "SOX2", "POU5F1", "NANOG", "ESR1", "PGR")
plot_differential_enrichment("Run10", "chromvar", "FT", "Tumor", names(colors.customtype2), colors.customtype2, unitlabels)
plot_differential_enrichment("Run10", "remap", "FT", "Tumor", names(colors.customtype2), colors.customtype2, unitlabels)
plot_differential_enrichment("Run10", "chromvar", "BRCAwt", "BRCAmt", names(colors.customtype2), colors.customtype2, unitlabels)
plot_differential_enrichment("Run10", "remap", "BRCAwt", "BRCAmt", names(colors.customtype2), colors.customtype2, unitlabels)
plot_differential_enrichment("Run10", "chromvar", "Naive", "NACT", names(colors.customtype2), colors.customtype2, unitlabels)
plot_differential_enrichment("Run10", "remap", "Naive", "NACT", c("Epithelial_cells"), colors.customtype2["Epithelial_cells"], unitlabels)
plot_differential_enrichment("Run10", "chromvar", "Sensitive", "Resistant", names(colors.customtype2), colors.customtype2, unitlabels)
plot_differential_enrichment("Run10", "remap", "Sensitive", "Resistant", c("Epithelial_cells"), colors.customtype2["Epithelial_cells"], unitlabels)
plot_differential_enrichment("Run10", "chromvar", "FT", "Naive", names(colors.customtype2), colors.customtype2, unitlabels)
plot_differential_enrichment("Run10", "remap", "FT", "Naive", c("Epithelial_cells"), colors.customtype2["Epithelial_cells"], unitlabels)

plot_differential_enrichment("Run10_filter2", "remap", "Naive", "NACT", c("Epithelial_cells"), colors.customtype2["Epithelial_cells"], unitlabels)
plot_differential_enrichment("Run10_filter2", "remap", "Sensitive", "Resistant", c("Epithelial_cells"), colors.customtype2["Epithelial_cells"], unitlabels)
plot_differential_enrichment("Run10_filter3", "remap", "Naive", "NACT", c("Epithelial_cells"), colors.customtype2["Epithelial_cells"], unitlabels)
plot_differential_enrichment("Run10_filter3", "remap", "Sensitive", "Resistant", c("Epithelial_cells"), colors.customtype2["Epithelial_cells"], unitlabels)

Tumor$remap$data[is.na(Tumor$remap$data)] <- 0
Tumor$remap$data[Tumor$remap$data == Inf] <- 0
differential_enrichment(Tumor, "remap", "customtype2", "subgroup2", "Naive", "NACT", "Run10test")
plot_differential_enrichment("Run10test", "remap", "Naive", "NACT", c("Epithelial_cells"), colors.customtype2["Epithelial_cells"], unitlabels)

################################# Define Persister Signature #########################################

Epithelial <- subset(filtered, subset = customtype2 == "Epithelial_cells")
TumorEpithelial <- subset(Epithelial, subset = group == "Tumor")

NaivNACT_R <- read.table("Run10_differential_remap_NACT_Naive.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct_1", "pct_2", "p_val_adj")) %>% filter(p_val_adj < 0.05) %>% drop_na()
signature <- NaivNACT_R %>% filter(celltype == "Epithelial_cells", p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene)
signature <- intersect(signature, rownames(filtered$RNA))
labs <- c("SOX17", "MECOM", "PAX8", "WT1", "FOXK2", "SOX2", "POU5F1", "NANOG", "ESR1", "PGR")

signature_expr <- TumorEpithelial$RNA$counts[signature,] %>% t() %>% data.frame() %>% mutate(group = TumorEpithelial$subgroup2) %>%
	pivot_longer(!group, names_to = "gene", values_to = "expr") %>%
	mutate(expressed = case_when(expr > 0 ~ 1, TRUE ~ 0))
NeoExprSignature <- signature_expr %>% group_by(gene, group) %>% 
	summarize(percexpr = sum(expressed)*100/n()) %>%
	pivot_wider(names_from = "group", values_from = "percexpr") %>%
	filter(NACT > 1) %>% pull(gene)

ord_table <- TumorEpithelial$remap$data[signature,] %>% t() %>% data.frame() %>% mutate(group = TumorEpithelial$subgroup2) %>%
	pivot_longer(!group, names_to = "gene", values_to = "score") %>%
	pivot_wider(names_from = "group", values_from = "score", values_fn = {mean}) %>%
	arrange(NACT-Naive) %>% mutate(is_signature = (NACT > 0 & Naive < 0 & gene %in% NeoExprSignature))
ord <- ord_table %>% pull(gene)
signature <- ord_table %>% filter(is_signature) %>% pull(gene) %>% sort()

SensResi_R <- read.table("Run10_differential_remap_Resistant_Sensitive.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj")) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::select(gene, celltype, avg_diff) %>% drop_na() %>% dplyr::filter(celltype =="Epithelial_cells")  %>% dplyr::select(gene, avg_diff)
NaivNACT_R <- read.table("Run10_differential_remap_NACT_Naive.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj")) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::select(gene, celltype, avg_diff) %>% drop_na() %>% dplyr::filter(celltype =="Epithelial_cells") %>% dplyr::select(gene,  avg_diff)

sharedEnriched_NvsN_SvsR <- dplyr::inner_join(NaivNACT_R, SensResi_R, by = join_by(gene == gene)) %>% drop_na() %>%
	dplyr::filter((avg_diff.x > 0 & avg_diff.y > 0) | (avg_diff.x < 0 & avg_diff.y < 0 )) 

Double_sig <- sharedEnriched_NvsN_SvsR %>% filter(gene %in% signature) %>% dplyr::pull(gene) %>% sort()

Epith_Tumor_figR <- read.table("Run10_Epith_Tumor_figR_motif_relaxed_sub.tsv", header = TRUE) %>% 
	set_colnames(c("DORC", "Motif", "Score")) %>%
	mutate(Score = as.numeric(Score))
signature_remap <- TumorEpithelial$remap$data[Double_sig,] %>% t() %>% data.frame() %>% mutate(group = TumorEpithelial$subgroup2) %>%
	pivot_longer(!group, names_to = "gene", values_to = "score") %>%
	left_join(ord_table %>% dplyr::select(gene, is_signature)) %>%
	mutate(gene = factor(gene, levels = ord)) %>%
	filter(is_signature) %>%
	mutate(FigRinfo = factor(case_when(gene %in% Epith_Tumor_figR$Motif ~ "RegTF", TRUE ~ "DBP"), levels = c("RegTF", "DBP")))
labeltable <- signature_remap %>% filter(group == "NACT", gene %in% OvSp) %>% group_by(gene) %>% summarize(score = mean(score)) %>%
	mutate(FigRinfo = factor(case_when(gene %in% Epith_Tumor_figR$Motif ~ "RegTF", TRUE ~ "DBP"), levels = c("RegTF", "DBP")))

pdf(file=paste0('Run10_remap_distribution.pdf'), height = 5, width = 5)
p2 <- ggplot(data = signature_remap, aes(x = gene, y = score)) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	# stat_summary(aes(color = group), fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange", size = 0.5, linewidth = 0.1, position=position_dodge(width = 0.5))+
	stat_summary(aes(color = group), fun.data="mean_sdl",  fun.args = list(mult=1), geom = "point", size = 1.5, position=position_dodge(width = 0.5))+
	# geom_text_repel(data = labeltable, mapping = aes(label = gene), size=3, hjust = 0, force = 1.5, nudge_y = 4, direction = "x", segment.size = 0.2) +
	geom_text_repel(data = labeltable, mapping = aes(label = gene), size=4, hjust = 0, force = 1.5, nudge_y = 1, direction = "x", segment.size = 0.2) +
	scale_color_manual(values = c("#0C8A8A", "#407280")) + 
	coord_cartesian(ylim = c(-3.8, 3.8)) +
	facet_grid(cols = vars(FigRinfo), scales = "free", space = "free") +
	theme_cowplot() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
p2
dev.off()

signature_remap %>% dplyr::select(gene, FigRinfo) %>% distinct() %>% pull(FigRinfo) %>% table()

Epith_Tumor_figR <- read.table("Run10_Epith_Tumor_figR_motif_relaxed_sub.tsv", header = TRUE) %>% 
	set_colnames(c("DORC", "Motif", "Score")) %>%
	mutate(Score = as.numeric(Score))
signature_remap <- TumorEpithelial$remap$data[signature,] %>% t() %>% data.frame() %>% mutate(group = TumorEpithelial$subgroup2) %>%
	pivot_longer(!group, names_to = "gene", values_to = "score") %>%
	left_join(ord_table %>% dplyr::select(gene, is_signature)) %>%
	mutate(gene = factor(gene, levels = ord)) %>%
	filter(is_signature) %>%
	mutate(FigRinfo = factor(case_when(gene %in% Epith_Tumor_figR$Motif ~ "RegTF", TRUE ~ "DBP"), levels = c("RegTF", "DBP")))
# labeltable <- signature_remap %>% filter(group == "NACT", gene %in% labs) %>% group_by(gene) %>% summarize(score = mean(score)) %>%
	# mutate(FigRinfo = factor(case_when(gene %in% Epith_Tumor_figR$Motif ~ "RegTF", TRUE ~ "DBP"), levels = c("RegTF", "DBP")))

pdf(file=paste0('Run10_remap_distribution_partial.pdf'), height = 6, width = 8)
p2 <- ggplot(data = signature_remap, aes(x = gene, y = score)) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	# stat_summary(aes(color = group), fun.data="mean_sdl",  fun.args = list(mult=1), geom = "pointrange", size = 0.1, position=position_dodge(width = 0.5))+
	stat_summary(aes(color = group), fun.data="mean_sdl",  fun.args = list(mult=1), geom = "point", size = 1.5, position=position_dodge(width = 0.5))+
	# geom_text_repel(data = labeltable, mapping = aes(label = gene), size=3, hjust = 0, force = 1.5, nudge_y = 4, direction = "x", segment.size = 0.2) +
	scale_color_manual(values = c("#0C8A8A", "#407280")) + 
	coord_cartesian(ylim = c(-3.8, 3.8)) +
	facet_grid(cols = vars(FigRinfo), scales = "free", space = "free") +
	theme_cowplot() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
p2
dev.off()
signature_remap %>% dplyr::select(gene, FigRinfo) %>% distinct() %>% pull(FigRinfo) %>% table()
signature_remap %>% dplyr::select(gene, FigRinfo) %>% distinct() %>% arrange(FigRinfo) %>% 
	write.table("Run10_Fig3d.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

################################### Figures and Venn Diagrams for Schematic Diagrams ################################

SensResi_R <- read.table("Run10_differential_remap_Resistant_Sensitive.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj")) %>% dplyr::filter(p_val_adj < 0.05) %>% drop_na() %>% dplyr::filter(celltype =="Epithelial_cells") %>% mutate(neg_log10_adj_pval = -log10(p_val_adj)) %>% mutate(neg_log10_adj_pval = case_when(neg_log10_adj_pval == Inf ~ 310, TRUE ~ neg_log10_adj_pval))
NaivNACT_R <- read.table("Run10_differential_remap_NACT_Naive.tsv", sep = "\t", skip = 1) %>% set_colnames(c("celltype", "gene", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj")) %>% dplyr::filter(p_val_adj < 0.05) %>% drop_na() %>% dplyr::filter(celltype =="Epithelial_cells") %>% mutate(neg_log10_adj_pval = -log10(p_val_adj)) %>% mutate(neg_log10_adj_pval = case_when(neg_log10_adj_pval == Inf ~ 310, TRUE ~ neg_log10_adj_pval))
pdf(file=paste0('Run10_epigenetic_signature_define.pdf'), height = 3, width = 8)
p1 <- ggplot(data = NaivNACT_R, aes(x = avg_diff, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(size=0.8, stroke=0.1, shape=16, color = colors.customtype2["Epithelial_cells"]), dpi = 400) +
	coord_cartesian(xlim = c(-3, 3), ylim = c(0, 310)) +
	theme_cowplot()
p2 <- ggplot(data = SensResi_R, aes(x = avg_diff, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(size=0.8, stroke=0.1, shape=16, color = colors.customtype2["Epithelial_cells"]), dpi = 400) +
	coord_cartesian(xlim = c(-3, 3), ylim = c(0, 310)) +
	theme_cowplot()
p1 | p2
dev.off()

pdf(file=paste0('Run10_epigenetic_signature_define1.pdf'), height = 4, width = 6)
# v <- venneuler(c(A=669, B=974, "A&B"=594))
s1 <- list(NACT = (NaivNACT_R %>% filter(avg_diff > 0) %>% pull(gene)), Resistant = (SensResi_R %>% filter(avg_diff > 0) %>% pull(gene)))
plot(euler(s1, shape = "ellipse"), quantities = TRUE)
dev.off()

intersect(NaivNACT_R %>% filter(avg_diff > 0) %>% pull(gene), SensResi_R %>% filter(avg_diff > 0) %>% pull(gene)) %>% length()
setdiff(NaivNACT_R %>% filter(avg_diff > 0) %>% pull(gene), SensResi_R %>% filter(avg_diff > 0) %>% pull(gene)) %>% length()
setdiff(SensResi_R %>% filter(avg_diff > 0) %>% pull(gene), NaivNACT_R %>% filter(avg_diff > 0) %>% pull(gene)) %>% length()
OvSp <- c("WT1", "EMX2", "SOX17", "MEIS1", "BHLHE41", "PAX8", "ESR1", "ZNF503", "MECOM", "TGIF2", "NR2F6", "PBX1", "ZNF217", "PLSCR1")
intersect(OvSp, intersect(NaivNACT_R %>% filter(avg_diff > 0) %>% pull(gene), SensResi_R %>% filter(avg_diff > 0) %>% pull(gene)))

########################### Plot Persister Scores Violin #########################################

Tumor_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & group == "Tumor"))
remove_factors <- union(names(which(rowSums(is.na(Tumor_Epi$remap$data))>0)), names(which(rowSums(Tumor_Epi$remap$data == Inf)>0)))
remove_factors <- union(remove_factors, setdiff(rownames(filtered$remap), rownames(filtered$RNA)))
remove_factors <- setdiff(remove_factors, "SOX2-OV")
DefaultAssay(Tumor_Epi) <- "remap"
Tumor_Epi <- DietSeurat(Tumor_Epi, assays = "remap")
Tumor_Epi$remap$data <- as(Tumor_Epi$remap$data, 'sparseMatrix')
Tumor_Epi <- subset(Tumor_Epi, features = sort(setdiff(rownames(Tumor_Epi), remove_factors)))
Tumor_Epi$sample <- factor(Tumor_Epi$sample, levels = c("PH1199", "PH1230", "PH1239", "PH1243", "PH1302", "PH1303", "PH1214", "PH1238", "PH1252"))

Tumor_Epi <- AddModuleScore(
  object = Tumor_Epi,
  features = list(Double_sig),
  assay = "remap",
  name = 'sig_remap',
  ctrl = 10,
  seed = 123
)

# table1 <- Tumor_Epi@meta.data %>% dplyr::select(sample, subgroup3, subgroup5, sig_remap1, epicluster) %>%
	# filter(sig_remap1 != Inf & sig_remap1 != -Inf)
table1 <- Tumor_Epi@meta.data %>% dplyr::select(sample, subgroup3, sig_remap1, epicluster) %>%
	filter(sig_remap1 != Inf & sig_remap1 != -Inf)

pdf(file=paste0(runid, "_VlnPlot_persister_ModuleScore.pdf"), width=14, height=6)
p3 <- ggplot(table1, aes(x=subgroup3, y=sig_remap1)) + 
    geom_violin(colour = "black") +
	stat_compare_means(aes(label = ..p.signif..), comparisons = list(c("Sensitive", "Resistant"), c("Resistant", "NACT"))) +
	ggrastr::rasterise(geom_jitter(size=0.001, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
	stat_summary(fun=mean, geom="point", size=0.7, color="red") +
    coord_cartesian(ylim = c(-5, 11)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 <- ggplot(table1, aes(x=sample, y=sig_remap1)) + 
    geom_violin() +
	ggrastr::rasterise(geom_jitter(size=0.1, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
	stat_summary(fun=mean, geom="point", size=0.7, color="#F61A23") +
    facet_grid(cols = vars(subgroup3), scales = "free", space = "free") +
	coord_cartesian(ylim = c(-11, 17)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p5 <- ggplot(table1, aes(x=epicluster, y=sig_remap1)) + 
    geom_violin(colour = "black") +
	ggrastr::rasterise(geom_jitter(size=0.001, stroke=0.1, shape=16), dpi=400) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
	stat_summary(fun=mean, geom="point", size=0.7, color="#F61A23") +
    coord_cartesian(ylim = c(-5, 11)) +
	theme_cowplot() + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 + p4 + p5 + plot_layout(widths = c(2, 4, 3))
dev.off()

############################ Correlation analysis ##################################
########################## Correlations with FULL cleanlist (ReMap) ################
DBFs <- sort(rownames(filtered$remap))
cleanlist <- sort(intersect(DBFs, rownames(filtered$RNA)))
filtered$remap$data[is.na(filtered$remap$data)] <- 0
filtered$remap$data[filtered$remap$data == Inf] <- 0

library(correlation)
DBFs <- read.csv("defined_DBFs_65.csv")
TFs <- read.csv("defined_TFs_113.csv")
signature_details <- bind_rows(TFs, DBFs)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
col_fun2 = gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
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
load("Run10_remap_smooth_groups_ALL.RData")

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
	ha = HeatmapAnnotation(isSig = (colnames(orderedcorrs) %in% Double_sig), col = list(isSig = c(`TRUE` = "black", `FALSE` = "white")), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
	# ra = rowAnnotation(isSig = (colnames(orderedcorrs) %in% Double_sig), col = list(isSig = c(`TRUE` = "black", `FALSE` = "white")), simple_anno_size = unit(5, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
	pdf(file=paste0("Run10_remap_Pearson_", label, ".pdf"), width=4, height=4)
		plot(Heatmap(matrix=orderedcorrs, show_row_dend=FALSE, show_column_dend=FALSE, cluster_rows=FALSE, cluster_columns=FALSE, col=col_fun, show_row_names=FALSE, show_column_names=FALSE, show_heatmap_legend=FALSE, left_annotation=rowAnnotation(link = anno_mark(at = match(labs, colnames(orderedcorrs)), labels = labs, which = "row", side = "left", labels_gp = gpar(fontsize = 60))), top_annotation=ha, use_raster = TRUE, raster_by_magick = TRUE))
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

pearsons_Control_Epi <- forcor_Control_Epi %>% do(data.frame(t(cor(., .))))
pearsons_Treated_Epi <- forcor_Treated_Epi %>% do(data.frame(t(cor(., .))))
pearsons_PH27C_Epi <- forcor_PH27C_Epi %>% do(data.frame(t(cor(., .))))
pearsons_PH626C_Epi <- forcor_PH626C_Epi %>% do(data.frame(t(cor(., .))))
pearsons_PH27T_Epi <- forcor_PH27T_Epi %>% do(data.frame(t(cor(., .))))
pearsons_PH626T_Epi <- forcor_PH626T_Epi %>% do(data.frame(t(cor(., .))))

plot_corr_heatmap(pearsons_Control_Epi, "ALL_Control_Epi_cc")
plot_corr_heatmap(pearsons_Treated_Epi, "ALL_Treated_Epi_cc")
plot_corr_heatmap(pearsons_PH27C_Epi, "ALL_PH27C_Epi_cc")
plot_corr_heatmap(pearsons_PH626C_Epi, "ALL_PH626C_Epi_cc")
plot_corr_heatmap(pearsons_PH27T_Epi, "ALL_PH27T_Epi_cc")
plot_corr_heatmap(pearsons_PH626T_Epi, "ALL_PH626T_Epi_cc")

corrs <- pearsons_NACT_Epi
dists <- as.dist(1 - corrs)
tree <- hclust(dists, method="complete")
dend <- as.dendrogram(tree)
labels(dend)[!(labels(dend) %in% Double_sig)] <- ""
# orderedcorrs <- corrs[order.dendrogram(dend), order.dendrogram(dend)]

pdf(file="Run10_remap_Pearson_ALL_NACT_Epi_dend.pdf", width=200, height=6)
plot(dend)
dev.off()

corrs <- pearsons_PH27T_Epi
dists <- as.dist(1 - corrs)
tree <- hclust(dists, method="complete")
dend <- as.dendrogram(tree)
labels(dend)[!(labels(dend) %in% Double_sig)] <- ""
# orderedcorrs <- corrs[order.dendrogram(dend), order.dendrogram(dend)]

pdf(file="Run10_remap_Pearson_ALL_PH27T_Epi_cc_dend.pdf", width=200, height=6)
plot(dend)
dev.off()

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

FT_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & subgroup2 == "FT"))
Naive_Epi <- subset(filtered, subset = (customtype2 == "Epithelial_cells" & subgroup2 == "Naive"))

forcor_FT_Epi <- smooth_scores(FT_Epi)
forcor_Naive_Epi <- smooth_scores(Naive_Epi)
save(forcor_FT_Epi, forcor_Naive_Epi, file="Run10_remap_smooth_fig2.RData")
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

################################ Compare with known genes / Upset Plots ############################

degsr_NACT_Naive <- read.table(file = "Run10_degsr_NACT_Naive.tsv", sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells")
degsr_Resistant_Sensitive <- read.table(file = "Run10_degsr_Resistant_Sensitive.tsv", sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells")

Verhaak_BadProg <- c("OR3A1", "HIF3A", "MARK4", "CCDC49", "PLXNA1", "MED25", "DCP2", "TNFSF11", "KPNA3", "B4GALT5", "CTNNBL1", "PAK4", "MUC4", "GAGE1", "CALB1", "MLL4", "DLGAP4", "SUPT5H", "FLJ10241", "DBF4B", "ACTR5", "CLTCL1", "KRT4", "RP11-125A7.3", "ANXA4", "GRB7", "SASH1", "PERLD1", "MTRF1", "DHX34", "PHF20", "PPM2C", "APC", "ZFP36", "EFNA1", "TTC31", "ERBB2", "CDH16", "RAB4B", "ZNF611", "GTF2F2", "ZNF331", "CASC3", "IRF2BP1", "RGS19", "ELMO2", "XRCC4", "C19ORF7", "WIPF2", "LTBP1", "RNF44", "RPL23", "TGIF2", "MAPT", "DMWD", "UTP6", "NEK3", "GMEB2", "SPAG5", "PTPN1", "ATP1A3", "BLCAP", "DHX35", "STRN4", "HIGD1B", "USH1C", "TBP", "TRIM13", "ITPKC", "KRT13", "ADA", "CEMP1", "SLC25A30", "HHLA2", "SAMD4B", "YTHDC2", "ZHX3", "AKT2", "GEMIN7", "SH3TC2", "STK4", "ANKRD46", "PSMC4", "HSPBAP1", "NUFIP1", "TCF15", "C20ORF4", "GRWD1", "FLJ20323", "C20ORF3", "SMG5", "C17ORF63", "KCNG1", "LRRC31", "DMPK", "CCNF", "ACVR2A", "FBXL18", "CAMK2G", "EML2", "MGRN1", "ZNF574", "MRPS31", "NLK", "C20ORF121", "PABPC1", "PLAC4", "CLK2")

Verhaak_GoodProg <- c("ALDH5A1", "GJB1", "AADAC", "DMC1", "GEMIN8", "ZFR", "HSPA1L", "PIGP", "DAP", "IL2", "PCK2", "GMPR", "HIST1H2BE", "CXCL9", "CECR1", "GSTZ1", "ETV4", "THEM2", "HLA-DOB", "RBP1", "HIRIP3", "FJX1", "FAH", "CD27", "HLA-DMB", "RTN3", "SHMT2", "CUTA", "ARL6IP5", "CREB3", "FLOT1", "WARS", "TXNDC13", "PHKG2", "CLEC11A", "CXCL13", "HIST1H2BO", "EDNRB", "HIST1H2AM", "BCL7C", "HIST1H2AB", "MTHFS", "SLC1A4", "HIST1H3E", "MAST4", "TRIM27", "GAD1", "TCP11", "SLC7A11", "FAM8A1", "HLCS", "UBD", "MBNL3", "ICOS", "BTN2A2", "ZNF76", "BSCL2", "TMEM30B", "WWOX", "HYAL2", "ID4", "EFTUD1", "BMPR2", "STK16", "PHLDA3", "GALNT6", "DNAJC4", "MATN2", "PNPLA4", "CUL4B", "HIST1H4H", "LEF1", "HLA-DOA", "SPRY2", "BTN3A1", "C14ORF159", "HIST1H2BI", "ZNF184", "NAB1", "C6ORF64", "FLJ14213", "NDUFC2", "CSRP1", "COMMD4", "NOTCH4")

Verhaak_differentiated <- c("FUT2", "STEAP3", "LCN2", "HRBL", "MYO5C", "ALS2CL", "ANXA1", "PTGER2", "MVP", "SSH3", "C1ORF116", "DTX4", "MLPH", "TNFRSF14", "MGLL", "CLU", "TMEM24", "SLC37A1", "RIN1", "SCGB2A2")
Verhaak_immunoreactive <- c("FBXO5", "LAG3", "AIM2", "CXCL10", "IFI30", "CASP1", "ECGF1", "LAIR1", "IL21R", "SLC31A2", "PDZK1IP1", "SAMSN1", "FCER1G", "RARRES3", "ADAMDEC1", "SLA", "LAPTM5", "UBE2L6", "KMO", "CIITA", "GIMAP5", "APOBEC3G")
Verhaak_mesenchymal <- c("PCOLCE", "MXRA8", "MAPRE2", "CXCR7", "MMP14", "PDPN", "PLEKHO1", "TPST1", "MARCKS", "CALD1", "TMEPAI", "NNMT", "GAS7", "WIPF1", "FSCN1", "NBL1", "LAMB1", "CACNA1C", "DLC1", "RUNX1", "KIAA0247", "PALLD", "RHOBTB3", "NUAK1", "MYO1B", "COL5A1", "COL8A1", "SRPX2", "LGALS1", "CTGF", "DDR2", "EVI2A", "MCAM", "FN1", "KDELC1", "STCH", "F2R")
Verhaak_proliferative <- c("COPS3", "EFS", "MARCKSL1", "RAD51AP1", "PKIA", "SALL2", "NETO2", "SPATS2", "MAPRE1", "UCHL1", "MFAP2", "TCF7L1", "LRP4", "BEX1", "KIF1A", "COL4A6", "TRO", "BLMH", "STMN1", "SMARCD1", "C5ORF13")

Millstein_OS <-c("TAP1", "ZFHX4", "CXCL9", "FBN1", "PTGER3")

Millstein_OS <- c("CDH1", "STK16", "MAL", "GJB1", "TESK1", "PTH2R", "DNAJC9", "SRI", "WWP1", "AKT1S1", "MYOD1", "NF1", "CPNE1", "BNIP3L", "NUCB2", "AADAC", "MITF", "CEACAM5", "GMNN", "ATP5A1", "C19orf12", "BRCA2", "GMPR", "SMARCA4", "PCDH9", "MAK", "PCK2", "GFRA1", "BAALC", "RNASEL", "B4GALT5", "MYC", "LPAR3", "DUSP4", "CDKN3", "E2F6", "FAM58A", "PARP4", "KRT6", "FOXJ1", "HSP90AA1", "USP8", "HIF1A", "SMO", "CTLA4", "CCNE1", "ASRGL1", "FGFR1", "IL22", "PPP2R4", "ZFHX4", "SPTLC2", "PD.1", "FOXP3", "RB1", "CXCL10", "PAX2", "OPA1", "EZR", "OASL", "RARRES1", "ANXA4", "HBB", "SERPINE1", "MINPP1", "APPL2", "MRPS27", "MDM2", "CDK6", "CRABP2", "SHPRH", "WDR91", "NTRK2", "GUSB", "OR1G1", "TAP1", "ESR2", "IGHM", "CX3CR1", "MRE11A", "VCAN", "GTF2H5", "MEST", "IGF2", "KDM5D", "TCF7L1", "VSIG4", "NF2", "FOXRED2", "MAP2K4", "ADH1B", "TBX2", "NUAK2", "ESD", "FGF1", "ZNHIT2", "PGRA", "KLHL7", "SOX17", "TSHR", "CXCL9")

Matondo_Chemoresponse <- c("ABCB9", "ACTA2", "ACTC1", "ACTG2", "ADH1C", "ALDH1A3", "AOC3", "APOD", "ASAP3", "ASPN", "BICD1", "CABYR", "CES1", "CILP", "COL10A1", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A3", "COL8A1", "CRISPLD2", "CRYAB", "CSAD", "CTGF", "CYP26A1", "CYR61", "DSC3", "DUSP1", "ECM2", "EGR2", "EPOR", "F3", "FBLN5", "FBN1", "FGF7", "FILIP1L", "FOS", "GPR124", "HOPX", "HOXA10", "IFFO1", "IGF1", "IGFBP7", "INHBA", "ITIH3", "JAM3", "KCNJ8", "KCNN4", "KDELR3", "KLK8", "LAMA1", "LAMA4", "LAMB1", "LMOD1", "LTBP2", "LUM", "MATN3", "MCM10", "MMP11", "MN1", "NBL1", "NETO2", "NR4A1", "NT5E", "NUAK1", "OMD", "PALLD", "PDGFC", "PDGFRA", "PDGFRB", "PDLIM7", "PDPN", "PKP2", "PLAU", "PLS3", "PMEPA1", "RAI14", "RASSF2", "ROR2", "RTP4", "RUNX1T1", "SCG2", "SEMA3C", "SFRP4", "SLC7A6", "SPINK1", "SST", "TAGLN", "THBS1", "THBS2", "TIMP3", "TLK1", "WNT6", "ZFHX4", "ZGPAT")

Reddy_OvMTF <- c("WT1", "EMX2", "SOX17", "MEIS1", "BHLHE41", "PAX8", "ESR1", "ZNF503", "MECOM", "TGIF2", "NR2F6", "PBX1", "ZNF217", "PLSCR1", "ADNP", "AEBP1", "ARID1A", "ATF4", "ATF5", "ATF6B", "BAZ2A", "BCLAF1", "BHLHE40", "CIC", "EGR1", "ELF3", "FOS", "FOSL2", "FOXP4", "GATAD2A", "GLYR1", "GPBP1L1", "GTF2I", "GTF3A", "HIF1A", "HMGA1", "HMGN3", "HSF1", "JUN", "JUNB", "JUND", "KDM2A", "KDM5B", "KLF5", "KLF6", "MAZ", "MYC", "NFE2L1", "NFE2L2", "NFIX", "NME2", "PA2G4", "PBX2", "RBCK1", "RELA", "REPIN1", "RERE", "SKI", "SON", "SOX4", "SP1", "SRCAP", "SREBF2", "STAT1", "STAT2", "STAT3", "STAT6", "TFDP1", "TP53", "TSC22D1", "UBP1", "UBTF", "USF2", "WHSC1", "XBP1", "YBX1", "ZNF146", "ZNF207", "ZNF664")
OvSp <- c("WT1", "EMX2", "SOX17", "MEIS1", "BHLHE41", "PAX8", "ESR1", "ZNF503", "MECOM", "TGIF2", "NR2F6", "PBX1", "ZNF217", "PLSCR1")
MultiOv <- c("ADNP", "AEBP1", "ARID1A", "ATF4", "ATF5", "ATF6B", "BAZ2A", "BCLAF1", "BHLHE40", "CIC", "EGR1", "ELF3", "FOS", "FOSL2", "FOXP4", "GATAD2A", "GLYR1", "GPBP1L1", "GTF2I", "GTF3A", "HIF1A", "HMGA1", "HMGN3", "HSF1", "JUN", "JUNB", "JUND", "KDM2A", "KDM5B", "KLF5", "KLF6", "MAZ", "MYC", "NFE2L1", "NFE2L2", "NFIX", "NME2", "PA2G4", "PBX2", "RBCK1", "RELA", "REPIN1", "RERE", "SKI", "SON", "SOX4", "SP1", "SRCAP", "SREBF2", "STAT1", "STAT2", "STAT3", "STAT6", "TFDP1", "TP53", "TSC22D1", "UBP1", "UBTF", "USF2", "WHSC1", "XBP1", "YBX1", "ZNF146", "ZNF207", "ZNF664")

TCGA_CommonlyMutated <- c("TP53", "BRCA1", "CSMD3", "NF1", "CDK12", "FAT3", "GABRA6", "BRCA2", "RB1")

Shannon_ChemosensML <- c("CYTH3", "GALNT3", "S100A14", "ERI1")

Fang_Gu_deM_PFS <- c("AAED1", "ABI3BP", "ABLIM2", "ABRACL", "ACKR3", "ADAMTS16", "ADAMTS5", "AHSG", "AKAP10", "ALKBH8", "AMY2B", "ANKRD36", "AP3B2", "APAF1", "APITD1", "APOD", "PRR5-ARHGAP8", "ARHGAP8", "ARL4C", "ARMC9", "ASNS", "ATAT1", "ATF7", "ATP6V0E2", "BAK1", "BTBD8", "C15orf65", "C1orf112", "C1orf162", "C4orf3", "C7orf73", "C8A", "CA13", "CARS2", "CCDC112", "CCSAP", "CDS1", "CDYL", "CEP112", "CEP41", "CFB", "CHIC1", "CLASP2", "CLDN15", "CLEC1A", "CNIH4", "COMMD8", "CRYBG3", "CWC27", "CYP3A43", "DCAF4L1", "DGCR14", "DNAH8", "DPT", "DSTNP2", "DYNC1I1", "EIF2B1", "EMP1", "ETV1", "EYA3", "FAIM", "FAM133B", "FAM200A", "GLRX2", "HAT1", "HEATR6", "HFE", "HHAT", "HIST1H1D", "ICK", "IGDCC4", "IGFBP1", "IGSF10", "IQCK", "IQSEC2", "ISCA1", "ITIH3", "KBTBD2", "KCNA1", "KCNJ14", "KIF3A", "KLRD1", "CCBL2", "LIX1", "LSM10", "MIR647", "MRE11A", "MRPL19", "MRTO4", "MSH3", "MTFMT", "MYSM1", "NAA15", "NAA25", "NARF-IT1", "NDC1", "NDUFS3", "NEURL1B", "NIT2", "NRK", "NUDCD1", "NUDT16", "NUPL1", "OFD1", "OLFML3", "PAX1", "PBX2", "PGGT1B", "PINK1-AS", "PODXL", "PPP2R2D", "PRKG2", "PRR5L", "PSMG4", "RABGGTA", "RBM44", "RGMA", "RINT1", "RNF138", "RRP15", "RSC1A1", "SCARNA6", "SF3B4", "SFRP2", "SIX4", "SLC22A5", "SLIT2-IT1", "SNHG10", "SNHG17", "SNORA1", "SNORA12", "SNORA32", "SNORA49", "SNORD37", "SNORD67", "SNORD73A", "STK17A", "SUPT20H", "SYNGAP1", "TADA2B", "TARSL2", "TAS2R3", "TENM3", "TFEB", "TGFB3", "THAP6", "TMCC3", "TMEM218", "ZMYM6NB", "TMEM86B", "TMEM99", "TOMM40L", "TRAM1L1", "TRAPPC2", "TREX1", "TRPC1", "TSNAX", "TSNAX-DISC1", "TSSC2", "TTC30A", "UAP1", "VNN3", "VPS33B", "VTA1", "WIPF2", "WRB", "WRNIP1", "ZBTB46", "ZFP37", "ZIC3", "ZNF121", "ZNF17", "ZNF204P", "ZNF234", "ZNF350", "ZNF470", "ZNF518B", "ZNF655", "ZNF684", "ZNF684", "ZNF708", "ZNF718", "ZNF776", "ZNF816", "ZNF844", "ZNF85", "ZNRD1", "ZNRF3", "ZSCAN31", "RP1-283E3.8", "RP11-793H13.10", "RP11-282O18.3", "RP4-758J18.13", "RP11-122K13.12", "RP11-203J24.9", "AC007390.5", "AC021593.1", "CTC-429L19.3", "RP11-872D17.8", "RP13-39P12.2", "AC002398.9", "RP11-690G19.3", "RP13-516M14.1")

Fang_Gu_UP_PFS <- c("C5AR1", "OSTM1", "SNORA12", "SNORD67", "BTG2", "AC021593.1", "SNN", "NDEL1", "SMOX", "CSDC2", "HRH2", "COL9A2", "ATP6V1E1", "AC125232.1", "RGL2", "PELI1", "SNORD35B", "SNORA32", "MCTP1", "RIPK3", "NR2F6", "GRB2", "NAA60", "FKBP1B", "SAMSN1", "PTGER2", "FAM179A", "NAA60", "MICAL1", "FXYD2", "MT-TM", "SRPK2", "TELO2", "CCBL2", "FBXL16", "RN7SKP48", "RN7SL521P", "SNORD35A", "SYTL3", "PNOC", "MED17", "NUPL1", "ISLR", "CCL5", "MT-TI", "PBX2", "AKAP10", "GLTSCR1L", "TRPC1", "NFE2", "RASSF10", "SNORA4", "VNN3", "SNORA19", "PACSIN2", "RP11-347C12.2", "C1orf162", "SNORD73A", "MRPS18C", "TMCC3", "FOXJ2", "TRAFD1", "AP006621.6", "SNORA49", "RN7SKP124", "RNF141", "FAM136A", "PLA2G4A", "ZDHHC5", "CORO7", "ZNF14", "SLC45A4", "DCAF4L1", "ZNF668", "RPS3A", "LAMTOR1", "HELB", "SNX15", "PPDPF", "RP11-399J13.3", "KLRD1", "NUDT16", "SNORD105", "CBWD1", "TSC22D2", "GPR75-ASB3", "RHBDF1", "PNPLA8", "ZNF845", "DIRC2", "ILKAP", "EIF1AX", "CORO1B", "FAM115C", "RNPS1", "UPK3B", "ZNF547", "PRR5-ARHGAP8", "SNORD87", "IRF5", "HERPUD2", "HIST1H1D", "ZNF708", "UNK", "RAB3A", "FAM104A", "PREP", "MRPS31", "METRNL", "RP13-516M14.1", "LGR6", "ZMYM6NB", "SPAG8", "RP11-894J14.5", "ARL4C", "OGFOD2", "HEATR6", "RP11-452L6.1", "RP11-815J21.1", "WDR65", "ACVRL1", "HIST1H2AJ", "SHPK", "ZNF655", "SCN2B", "FTSJ1", "SNORD18C", "SLC25A12", "BCL9", "CLDN15", "LYPD1", "RN7SL483P", "FHL3", "TADA2B", "PAXIP1", "MFSD10", "VPS53", "EPC2", "CCZ1", "TSNAX-DISC1", "SORCS3", "RRS1", "BOC", "ZC3HC1", "CDK2AP1", "DTX2P1", "MSLN", "UBIAD1", "APCDD1", "DNAJB6", "SNORD11", "KLHL11", "C5orf51", "PPIL6", "UBE2D4", "SCARNA6", "WWC3", "RP11-690G19.3", "ATP6V0E2", "SERTM1", "ATP8B2", "CDKN2D", "ZNF17", "SULF2", "VRK3", "ZER1", "RBM44", "SF3B4", "LINC00639", "SLC22A5", "C16orf59", "SYNGAP1", "NUP62", "DGCR14", "METRN", "ZBTB46", "NXF1", "TM2D2", "ABI3BP", "CRIP2", "DTD2", "MKL1", "SNORD121A", "HSPA1L", "DOT1L", "KIAA1731", "STXBP3", "APAF1", "C11orf30", "RP1-142L7.9", "RP4-639F20.1", "POM121C", "HIST1H1C", "GPR141", "ZNF226", "NFATC2", "SNORD5", "DNAH7", "HIST1H3J", "GLIS2", "AC002398.9", "ZNF582", "BAK1", "PGP", "RP4-758J18.13", "TGFBR3", "IGSF10", "RP11-268G12.1", "PLEKHF2", "CLTCL1", "ZNF446", "VAMP2", "WIPF2", "CLEC1A", "RGCC", "TRAF6", "RP11-712P20.2", "FAM222B", "TMEM136", "ABCA11P", "TIMM17B", "C6orf136", "RP11-793H13.10", "VCPKMT", "RP11-761B3.1", "ARHGAP8", "RTN4IP1", "TRAM1L1", "RP11-175P13.3", "ING4", "PPM1L", "STON1", "SNORA1", "AXIN2", "RIN1", "TMEM86B", "GSPT2", "KRT18P31", "RP3-368A4.5", "SNRNP48", "CTD-2228K2.5", "HID1", "SNHG10", "RP11-430C7.4", "LSM10", "RP11-77P16.4", "RN7SL735P", "RP11-815J4.6", "OSBPL7", "POLN", "RP11-357C3.3", "CA5B", "ARHGEF6", "RP11-3J1.1", "FRMD8", "MRPS10", "NDUFS3", "E2F8", "ITM2A", "S1PR1", "PVRIG2P", "RP11-293G6__A.3", "TMEM8B", "RN7SL35P", "MGAT3", "GCNT1P1", "ARMC2", "ALOX15", "MCF2L", "CBX8", "GOLGA8R", "RN7SL246P", "SNORD123", "TSNAXIP1", "SFRP2", "LAT", "RP11-564D11.3", "CDC42BPG", "SLC35E2", "TMEM130", "CDKL4", "SYNGR3", "NSUN5P1", "RN7SL577P", "RP11-453F18__B.1", "CCDC102A", "LRCH4", "CHCHD4", "RP11-1109F11.5", "NME8", "MIR647", "RPL7AP34", "COX6B1P4")

Fang_Gu_deM <- c("AASDHPPT", "ADAMTS13", "ART3", "BCL2L2", "BREA2", "C3orf35", "C5orf45", "C7orf55", "C9orf131", "CCND3P1", "CD180", "CLDN20", "CNTF", "DFFBP1", "DTX2P1-UPK3BP1-PMS2P11", "DUSP28", "FAM13A-AS1", "FBXW4P1", "FIZ1", "FTH1P16", "GVQW1", "HAMP", "HIST1H1T", "HMOX1", "KLRAP1", "LCMT1-AS2", "LDHAL6FP", "LETM2", "LINC01596", "LOC100129215", "LOC100130849", "LOC101927550", "LOC101928524", "LOC101930028", "LOC102723335", "LOC641367", "MC1R", "MEF2BNBP1", "mir-1205", "mir-1237", "mir-3130", "mir-4484", "mir-548", "mir-922", "MIR3164", "MIR4326", "MIR4655", "MIR4737", "MOV10", "OR2A1-AS1", "PHLDA1", "PINK1-AS", "POLR3G", "PRKXP1", "PTGES3P1", "PTPRA", "RBFADN", "RN7SL263P", "RN7SL307P", "RNA5SP260", "RNA5SP434", "RNA5SP488", "RNA5SP74", "RNF152P1", "RNU6-254P", "RNU6-520P", "RNU7-107P", "RPL35AP31", "RRM1-AS1", "SAPCD2P3", "SEMA6A-AS1", "SIRT7", "SLC25A1P5", "SNORD88A", "SNORD88B", "SPCS2P1", "SPDYE1", "STK31", "SYCE1L", "SYCP3", "TADA2B", "TAS2R13", "TECRP1", "TMEM209", "TMEM86B", "TREX1", "TTLL3", "TUBA3C/TUBA3D", "TUBB1", "TUBGCP2", "UBA5", "ZBTB20-AS1", "ZDHHC20-IT1", "ZMYND10")

stress_associated <- c("CEBPB", "CEBPD", "FOS", "IL6", "JUN", "JUNB", "MCL1", "MYC", "SOCS3", "ATF3", "DUSP1", "EGR1", "FOSB", "CEBPA", "DDIT3", "EGR2", "CDKN1A", "GADD45B", "TNF","HES1", "HBEGF", "BCL6", "NR4A1", "DUSP6", "GADD45G", "ID2", "NFKBIA", "PLK3", "SNAI2", "CREB5", "HLA-G", "HIST1H2BC", "HIST1H2BG", "CALML3", "SNAI1")

# Zhang et al 2022
EMT_associated <- c("DST","EGFR","EPHB2","ITGA6","PIK3CD","PLEC","SFN","SMAD3","CLTC", "EIF2AK2","ITGA2","KPNA1","STAT1","BCAR1","CCND1","COL1A2","ITCH", "MMP1","MMP12","PML","TNC","CDK6","FADD","FAS","MX1","E2F3","BIRC2","CAPN2","VASP","MMP9","DFFA","KPNB1","ADAMTS5","MMP10","AP2B1","ASAP1","SDCBP","CYP1A1","CYP1B1","YP3A5","MAFG","SLC7A11","IGF2")

# Xie et al 2022
Quiescence_associated <- c("TTYH1", "AQP4", "GFAP", "RGMA", "EDNRB", "CLU", "SLC1A3", "TSC22D4", "ATP1B2", "KCNN3", "FAM107A", "HEPN1", "APOE", "NMB", "CD99", "B2M", "TIMP3", "C1ORF61", "ZFP36L1", "NCAN", "ZFP36L2", "CRYAB", "IGFBP7", "PTN", "TMEM47", "DKK3", "GAS2L1", "PLCD3", "HEPACAM", "FGFR3", "GABARAPL2", "ALDOC", "SOX2", "PEA15", "PLPP3", "CXCL14", "SLC6A11", "GLIS3", "RARRES3", "HES5", "ATP1A2", "HES1", "CNN3", "QKI", "ID4", "GPR137B", "CRB2", "GRAMD3", "ID3", "DOK5", "PLTP", "NDP", "PCDH9", "RFX4", "PHKG1", "GLUD1", "EEPD1", "HLA.E", "HES4", "SLC4A4", "GLUL", "NDRG2", "LIFR", "MIR99AHG", "ADCYAP1R1", "ZBTB20", "NPC2", "PBXIP1", "NFIA", "ST5", "ID2", "AQP1", "RNF213", "CALD1", "GATM", "FAM69C", "PLA2G16", "TOB2", "ADAMTS6", "PAX6", "KIF1B", "PNISR", "FGFR2", "WIPF2", "DAZAP2", "SLC1A2", "C6ORF1", "DOCK7", "HSPB1", "SLC6A9", "CDC42EP4", "NEAT1", "P2RY1", "BMPR1B", "FMN2", "SARAF", "DTNA", "RSRP1", "LRP10", "PTTG1IP", "PTPRF", "POLR2J3", "NAT8L", "EMP3", "EFHC1", "ARHGEF26", "F3", "PNRC1", "FSTL1", "BDH2", "KAT2A", "PDLIM3", "PIK3C2A", "DDX17", "ALDH6A1", "APC", "WASL", "ACYP2")

differentiated <- c("B3GNT7", "B3GNT8", "GCNT3", "MUC15", "MUC4", "MUC5B", "ST3GAL4", "LCN2", "S100A8", "S100A9", "ALOX5", "AOC1", "B2M", "C3", "CEACAM1", "CEACAM6", "MGST1", "OSTF1", "PPBP", "S100P", "SLPI", "TMEM173", "DSC2", "EVPL", "IVL", "KRT13", "KRT4", "KRT6A", "KRT6C", "HRASLS2", "LPCAT4", "PLBD1", "RARRES3", "CDKN2B", "STEAP4", "TGFB1", "CFB", "CCL28", "CX3CL1", "CXCL10")

Proliferative_DNA_repair <- c("AURKA", "AURKB", "BIRC5", "BUB1", "BUB1B", "CCNA2", "CCNB1", "CCNB2", "CDC20", "CDC25C", "CDC45", "CDC6", "CDCA8", "CDK1", "CDT1", "CENPA", "CENPE", "CENPF", "CENPH", "CENPK", "CENPM", "CENPN", "CENPU", "CKAP5", "DHFR", "E2F1", "H2AFV", "H2AFZ", "HAUS1", "HIST1H2BH", "HIST1H3G", "HIST1H4C", "HIST2H2AC", "KIF18A", "KIF20A", "KIF23", "KIF2C", "LIG1", "LMNB1", "MAD2L1", "MCM10", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "MCM8", "NDC80", "NEK2", "NUF2", "ORC6", "PCNA", "PLK1", "PLK4", "POLD1", "POLD2", "PRIM1", "PTTG1", "RFC2", "RFC3", "RFC4", "RRM2", "SGO1", "SGO2", "SKA1", "SKA2", "SKP2", "SPC24", "SPC25", "TK1", "TMPO", "TPX2", "UBE2C", "VRK1", "ZWILCH", "ZWINT", "BRCA1", "CHEK1", "CLSPN", "RHNO1", "RMI2", "FANCG", "HMGB2", "KIF4A", "RRM1", "TYMS", "RACGAP1", "CBX5", "BRCA2", "FANCD2", "RAD51", "RAD51AP1", "FANCB", "FANCI", "UBE2T", "KIF11", "KIF15", "KIF20B", "KIF22", "KIFC1", "HIST1H1A", "HIST1H1B", "GGH", "DTYMK", "AKR1B1")

cytokine_apoptosis <- c("CCL20", "CXCL1", "IL1R2", "IL1RN", "SEC11C", "SPCS3", "BIK", "BIRC3", "CDKN2A", "DAB2", "GJB2")

AntigenPresentation <- c("CTSL", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DRA", "HLA-DRB5", "HSP90AB1", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "PDIA3", "CPE", "CTSF", "C1R", "CXADR", "SFTPD", "TUBB2A", "TUBB6", "FOLR1", "GAS6", "NPC2", "RUNX1", "TGFBR2", "GNAI1", "RBPJ", "IL11RA", "HIST1H2AE", "DNAJC3", "FBXO2", "HERPUD1", "HSP90B1", "SKP1", "UBE2D1", "UBQLN1", "NECTIN2", "SDC2", "CCND3", "DNAJB9", "ERP27", "SIAH1", "PPP2CA", "BAG3", "HSPA12A", "BMP2", "LATS2", "PPP2R2B", "PRKCI", "WNT6", "IFITM1", "IFITM2", "CRMP1", "DPYSL3", "PLXNA2", "MYL9", "TNNC1", "TNNI1", "TNNT2", "TLE1", "TLE4", "APP", "PELI2", "VAMP2", "DNAJA4", "PLAT", "ACADVL", "TSPYL2", "NR1D2", "PURA", "GATA6", "TWIST1", "RING1", "BAMBI", "ATG101", "RAPGEF2", "NR4A3", "SOX4", "SLC22A18")

Interferon_signaling <- c("GBP4", "IFI27", "IFI35", "IFIT1", "ISG15", "OAS1", "OAS2", "STAT2", "TRIM22", "TAP1", "ITGAV")

RNA_processing <- c("BMS1", "DCAF13", "MPHOSPH10", "PNO1", "ACIN1", "ROCK1", "TJP1", "ESCO1", "WAPL", "EIF1AX", "PNN", "RNPS1", "DYNC1LI1", "BDP1", "CRCP", "PAFAH1B1", "PRPF40A", "SFSWAP", "KMT2A", "SETD2")

Proteasomal_degradation <- c("AIMP1", "EIF3I", "MRPL13", "MRPL15", "MRPL16", "MRPL2", "MRPL28", "MRPL3", "MRPS7", "RPN2", "SRP9", "TSFM", "PSMA4", "PSMB5", "PSMB6", "PSMB7", "PSMC5", "PSMD8", "ACTL6A", "BANF1", "PPIA", "XRCC6", "VDAC3", "FH", "MDH2", "SDHB", "TPI1", "ECHS1", "HSD17B10", "POLR2G", "NDUFB6", "HNRNPA2B1", "HNRNPA3", "SNRPD3", "PRDX3", "PTMA", "CCT7", "VBP1", "COPE")

TCA_cycle <- c("GLO1", "IDH3B", "NDUFA10", "NDUFB5", "PDHA1", "SUCLG1", "TRAP1", "UQCRC1", "UQCRC2", "MDH1", "ALDH7A1", "DECR1", "ECH1", "ECI2", "PCCB", "SHMT1", "HACD3", "PTGR1", "MGST2", "PRDX6")

# cat(setdiff(Verhaak_BadProg108, rownames(filtered)))

NACT <- degsr_NACT_Naive %>% filter(avg_log2FC > 0) %>% pull(gene)
Naive <- degsr_NACT_Naive %>% filter(avg_log2FC < 0) %>% pull(gene)
Resistant <- degsr_Resistant_Sensitive %>% filter(avg_log2FC > 0) %>% pull(gene)
Sensitive <- degsr_Resistant_Sensitive %>% filter(avg_log2FC < 0) %>% pull(gene)

######

listInput4a <- list(NACT = NACT, Naive = Naive, Verhaak_BadProg = Verhaak_BadProg, Verhaak_GoodProg = Verhaak_GoodProg, Verhaak_differentiated = Verhaak_differentiated, Verhaak_immunoreactive = Verhaak_immunoreactive, Verhaak_mesenchymal = Verhaak_mesenchymal, Verhaak_proliferative = Verhaak_proliferative, Millstein_OS = Millstein_OS, Matondo_Chemoresponse = Matondo_Chemoresponse, OvSp = OvSp, MultiOv = MultiOv, stress_associated = stress_associated, Proliferative_DNA_repair = Proliferative_DNA_repair)

temp1 <- union(Naive, NACT)
temp2 <- Reduce(union, listInput4a[3:14])
length(intersect(temp1, temp2))

m = make_comb_mat(listInput4a)
m = m[comb_size(m) >= 2 & comb_degree(m) != 1]
m = m[str_starts(comb_name(m), "10|01")]

pdf(file=paste0('Run10_upset4a.pdf'), height = 5, width = 7)
# UpSet(m, set_order = order(set_size(m), decreasing=TRUE), comb_order = order(comb_size(m), decreasing=TRUE), top_annotation = upset_top_annotation(m, add_numbers = T), right_annotation = upset_right_annotation(m, add_numbers = T))
UpSet(m, set_order = order(set_size(m), decreasing=TRUE), comb_order = c(20, 16, 17, 15, 12, 19, 13, 18, 14, 3, 9, 4, 8, 5, 11, 1, 2, 6, 7, 10), top_annotation = upset_top_annotation(m, add_numbers = T), right_annotation = upset_right_annotation(m, add_numbers = T))
dev.off()

# temp <- sapply(order(comb_size(m), decreasing=TRUE), function(x){extract_comb(m, comb_name(m)[x])})
# cat(unlist(temp[1]), sep = "\n")
# extract_comb(m, comb_name(m)[12]) # grey70
# comb_name(m)[order(comb_size(m), decreasing=TRUE)]

degsr <- read.table(file = paste0("Run10", "_degsr_", "NACT", "_", "Naive", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells") %>%
	mutate(ovps = case_when(gene %in% union(extract_comb(m, comb_name(m)[16]), 
									  union(extract_comb(m, comb_name(m)[8]), extract_comb(m, comb_name(m)[3]))) ~ "Millstein_OS", 
							gene %in% union(extract_comb(m, comb_name(m)[17]), 
									  union(extract_comb(m, comb_name(m)[9]), extract_comb(m, comb_name(m)[1]))) ~ "Matondo_Chemoresponse",
							gene %in% extract_comb(m, comb_name(m)[18]) ~ "OvSp", 
							gene %in% union(extract_comb(m, comb_name(m)[19]), 
									  union(extract_comb(m, comb_name(m)[10]), extract_comb(m, comb_name(m)[2]))) ~ "MultiOv", 
							gene %in% union(extract_comb(m, comb_name(m)[12]), extract_comb(m, comb_name(m)[4])) ~ "Verhaak_BadProg",
							gene %in% union(extract_comb(m, comb_name(m)[13]), extract_comb(m, comb_name(m)[3])) ~ "Verhaak_GoodProg", 
							gene %in% extract_comb(m, comb_name(m)[15]) ~ "Verhaak_proliferative",
							gene %in% extract_comb(m, comb_name(m)[5]) ~ "Verhaak_differentiated",
							gene %in% extract_comb(m, comb_name(m)[6]) ~ "Verhaak_immunoreactive",
							gene %in% union(extract_comb(m, comb_name(m)[14]), 
									  union(extract_comb(m, comb_name(m)[7]), extract_comb(m, comb_name(m)[1]))) ~ "Verhaak_mesenchymal",
							gene %in% extract_comb(m, comb_name(m)[20]) ~ "Proliferative_DNA_repair", 
							gene %in% union(extract_comb(m, comb_name(m)[11]), extract_comb(m, comb_name(m)[2])) ~ "stress_associated",
							TRUE ~ "none"))
cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#860C1E", "#984EA3", "#FF7F00", "#F781BF", "#FFBF00", "#061A74", "#4eddff", "#AC8EFF", "#97D540")
names(cols) <- c("Millstein_OS", "Matondo_Chemoresponse", "OvSp", "MultiOv", "Verhaak_BadProg", "Verhaak_GoodProg", "Verhaak_proliferative", "Verhaak_differentiated", "Verhaak_immunoreactive", "Verhaak_mesenchymal", "Proliferative_DNA_repair", "stress_associated")

temp1 <- degsr %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))
    
pdf(file=paste0('Run10_upset4a_degsr.pdf'), height = 6, width = 8)
ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps == "none"), color = "grey70", size = 1, stroke=0.1, shape = 16), dpi = 400) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps != "none"), aes(color = ovps), size = 2, stroke=0.1, shape = 16), dpi = 400) +
	# geom_text_repel(data = temp1 %>% filter(ovps != "none"), aes(label = gene, segment.color=ovps), size = 2, max.overlaps = Inf, segment.size = 0.5, force = 3, fontface = "italic") +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	coord_cartesian(xlim = c(-10, 10)) +
	scale_color_manual(values = cols, aesthetics = c("color", "segment.color")) +
	labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(legend.title=element_blank())
dev.off()

temp1 <- read.table(file = paste0("Run10", "_degsr_", "NACT", "_", "Naive", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells") %>%
	mutate(ovps = case_when(gene %in% union(extract_comb(m, comb_name(m)[16]), 
									  union(extract_comb(m, comb_name(m)[8]), extract_comb(m, comb_name(m)[3]))) ~ "Millstein_OS", 
							TRUE ~ "none")) %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))
pdf(file=paste0('Run10_upset4a_degsr_Millstein_OS.pdf'), height = 6, width = 7)
ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps == "none"), color = "grey70", size = 1, stroke=0.1, shape = 16), dpi = 400) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps != "none"), aes(color = ovps), size = 2, stroke=0.1, shape = 16), dpi = 400) +
	geom_text_repel(data = temp1 %>% filter(ovps != "none"), aes(label = gene, segment.color=ovps), size = 4, max.overlaps = Inf, segment.size = 0.5, force = 3, fontface = "italic") +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	coord_cartesian(xlim = c(-10, 10)) +
	scale_color_manual(values = cols, aesthetics = c("color", "segment.color")) +
	labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(legend.position="none")
dev.off()
temp1 <- read.table(file = paste0("Run10", "_degsr_", "NACT", "_", "Naive", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells") %>%
	mutate(ovps = case_when(gene %in% union(extract_comb(m, comb_name(m)[17]), 
									  union(extract_comb(m, comb_name(m)[9]), extract_comb(m, comb_name(m)[1]))) ~ "Matondo_Chemoresponse",
							TRUE ~ "none")) %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))
pdf(file=paste0('Run10_upset4a_degsr_Matondo_Chemoresponse.pdf'), height = 6, width = 7)
ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps == "none"), color = "grey70", size = 1, stroke=0.1, shape = 16), dpi = 400) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps != "none"), aes(color = ovps), size = 2, stroke=0.1, shape = 16), dpi = 400) +
	geom_text_repel(data = temp1 %>% filter(ovps != "none"), aes(label = gene, segment.color=ovps), size = 4, max.overlaps = Inf, segment.size = 0.5, force = 3, fontface = "italic") +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	coord_cartesian(xlim = c(-10, 10)) +
	scale_color_manual(values = cols, aesthetics = c("color", "segment.color")) +
	labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(legend.position="none")
dev.off()
temp1 <- read.table(file = paste0("Run10", "_degsr_", "NACT", "_", "Naive", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells") %>%
	mutate(ovps = case_when(gene %in% union(extract_comb(m, comb_name(m)[12]), extract_comb(m, comb_name(m)[4])) ~ "Verhaak_BadProg",
							gene %in% union(extract_comb(m, comb_name(m)[13]), extract_comb(m, comb_name(m)[3])) ~ "Verhaak_GoodProg", 
							TRUE ~ "none")) %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))
pdf(file=paste0('Run10_upset4a_degsr_Verhaak_GoodBadProg.pdf'), height = 6, width = 7)
ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps == "none"), color = "grey70", size = 1, stroke=0.1, shape = 16), dpi = 400) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps != "none"), aes(color = ovps), size = 2, stroke=0.1, shape = 16), dpi = 400) +
	geom_text_repel(data = temp1 %>% filter(ovps != "none"), aes(label = gene, segment.color=ovps), size = 4, max.overlaps = Inf, segment.size = 0.5, force = 3, fontface = "italic") +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	coord_cartesian(xlim = c(-10, 10)) +
	scale_color_manual(values = cols, aesthetics = c("color", "segment.color")) +
	labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(legend.position="none")
dev.off()
temp1 <- read.table(file = paste0("Run10", "_degsr_", "NACT", "_", "Naive", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells") %>%
	mutate(ovps = case_when(gene %in% extract_comb(m, comb_name(m)[15]) ~ "Verhaak_proliferative",
							gene %in% extract_comb(m, comb_name(m)[5]) ~ "Verhaak_differentiated",
							gene %in% extract_comb(m, comb_name(m)[6]) ~ "Verhaak_immunoreactive",
							gene %in% union(extract_comb(m, comb_name(m)[14]), 
									  union(extract_comb(m, comb_name(m)[7]), extract_comb(m, comb_name(m)[1]))) ~ "Verhaak_mesenchymal",
							TRUE ~ "none")) %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))
pdf(file=paste0('Run10_upset4a_degsr_Verhaak_subtypes.pdf'), height = 6, width = 7)
ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps == "none"), color = "grey70", size = 1, stroke=0.1, shape = 16), dpi = 400) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps != "none"), aes(color = ovps), size = 2, stroke=0.1, shape = 16), dpi = 400) +
	geom_text_repel(data = temp1 %>% filter(ovps != "none"), aes(label = gene, segment.color=ovps), size = 4, max.overlaps = Inf, segment.size = 0.5, force = 3, fontface = "italic") +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	coord_cartesian(xlim = c(-10, 10)) +
	scale_color_manual(values = cols, aesthetics = c("color", "segment.color")) +
	labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(legend.position="none")
dev.off()
temp1 <- read.table(file = paste0("Run10", "_degsr_", "NACT", "_", "Naive", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells") %>%
	mutate(ovps = case_when(gene %in% extract_comb(m, comb_name(m)[18]) ~ "OvSp", 
							gene %in% union(extract_comb(m, comb_name(m)[19]), 
									  union(extract_comb(m, comb_name(m)[10]), extract_comb(m, comb_name(m)[2]))) ~ "MultiOv",
							TRUE ~ "none")) %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))
pdf(file=paste0('Run10_upset4a_degsr_Reddy.pdf'), height = 6, width = 7)
ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps == "none"), color = "grey70", size = 1, stroke=0.1, shape = 16), dpi = 400) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps != "none"), aes(color = ovps), size = 2, stroke=0.1, shape = 16), dpi = 400) +
	geom_text_repel(data = temp1 %>% filter(ovps != "none"), aes(label = gene, segment.color=ovps), size = 4, max.overlaps = Inf, segment.size = 0.5, force = 3, fontface = "italic") +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	coord_cartesian(xlim = c(-10, 10)) +
	scale_color_manual(values = cols, aesthetics = c("color", "segment.color")) +
	labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(legend.position="none")
dev.off()
temp1 <- read.table(file = paste0("Run10", "_degsr_", "NACT", "_", "Naive", ".tsv"), sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>%
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells") %>%
	mutate(ovps = case_when(gene %in% extract_comb(m, comb_name(m)[20]) ~ "Proliferative_DNA_repair", 
							gene %in% union(extract_comb(m, comb_name(m)[11]), extract_comb(m, comb_name(m)[2])) ~ "stress_associated",
							TRUE ~ "none")) %>% 
	mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
	mutate(neg_log10_adj_pval = -log10(p_val_adj))
pdf(file=paste0('Run10_upset4a_degsr_Zhang.pdf'), height = 6, width = 7)
ggplot(temp1, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps == "none"), color = "grey70", size = 1, stroke=0.1, shape = 16), dpi = 400) +
	ggrastr::rasterise(geom_point(data = temp1 %>% filter(ovps != "none"), aes(color = ovps), size = 2, stroke=0.1, shape = 16), dpi = 400) +
	geom_text_repel(data = temp1 %>% filter(ovps == "stress_associated"), aes(label = gene, segment.color=ovps), size = 4, max.overlaps = Inf, segment.size = 0.5, force = 3, fontface = "italic") +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	coord_cartesian(xlim = c(-10, 10)) +
	scale_color_manual(values = cols, aesthetics = c("color", "segment.color")) +
	labs(x = "log2(Fold change of average expression)", y = "-log10(Adjusted P value)") +
	guides(color = guide_legend(override.aes = list(size=5))) +
	theme_cowplot() + theme(legend.position="none")
dev.off()

################################# Explore targets of signature (FigR analysis) #########################
Epith_Tumor_figR <- read.table("Run10_Epith_Tumor_figR_motif_relaxed_sub.tsv", header = TRUE) %>% 
	set_colnames(c("DORC", "Motif", "Score")) %>%
	mutate(Score = as.numeric(Score))
degsr_FN <- read.table(file = "Run10_degsr_Naive_FT.tsv", sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells")
degsr_NN <- read.table(file = "Run10_degsr_NACT_Naive.tsv", sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells")
degsr_SR <- read.table(file = "Run10_degsr_Resistant_Sensitive.tsv", sep = "\t", skip = 1) %>%
    set_colnames(c("celltype", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) %>% 
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells")
dars_NN <- read.table(file = "Run10_differential_remap_NACT_Naive.tsv", sep = "\t", header=TRUE) %>% 
	filter(p_val_adj < 0.05 & celltype == "Epithelial_cells")
dars_SR <- read.table(file = "Run10_differential_remap_Resistant_Sensitive.tsv", sep = "\t", header=TRUE) %>% 
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

pdf(file="Run10_degsr_NACT_Naive_targets_partial.pdf", width=8, height=8)
p1 <- ggplot(degsr_sig, aes(x = avg_log2FC, y = neg_log10_adj_pval)) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(is.na(PropAct)), size = 1.5, stroke=0.1, shape = 16, color = "grey80"), dpi = 400) +
	ggrastr::rasterise(geom_point(data = degsr_sig %>% filter(!is.na(PropAct)), mapping = aes(color = PropAct), size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	# scale_color_distiller(palette = "RdYlGn", direction = 1) +
	scale_color_gradient2(low = "#d73027", mid = "#FEAE2F", high = "#1a9850", midpoint = 0.5) +
	theme_cowplot()
(p1 + annotation_custom(ggplotGrob(DOWN_pie_degsr), xmin = -10, xmax = -7, ymin = 220, ymax = 250) + annotation_custom(ggplotGrob(UP_pie_degsr), xmin = 7, xmax = 10, ymin = 220, ymax = 250))
dev.off()

write.table(temp1, "Run10_degsr_NACT_Naive_targets_partial.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

temp1 %>% filter(Expr == "UP" & Reg != "NoInfo") %>% pull(gene) %>% cat(sep = ",")
temp1 %>% filter(Expr == "DOWN" & Reg != "NoInfo") %>% pull(gene) %>% cat(sep = ",")

############################# GSEA / gprofiler2 ###############################

signature <- Double_sig
Epith_Tumor_figR <- read.table("Run10_Epith_Tumor_figR_motif_relaxed_sub.tsv", header = TRUE) %>% 
	set_colnames(c("DORC", "Motif", "Score")) %>%
	mutate(Score = as.numeric(Score))
degsr_NN <- read.table(file = "Run10_degsr_NACT_Naive.tsv", sep = "\t", skip = 1) %>%
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

rightgreengenes <- temp1 %>% filter(Expr == "UP" & Reg == "MoreAct") %>% pull(gene)
leftredgenes <- temp1 %>% filter(Expr == "DOWN" & Reg == "MoreRep") %>% pull(gene)

gostres1 <- gost(query = leftredgenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))
gostres2 <- gost(query = rightgreengenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))

# highlight_terms <- c("response to stimulus", "developmental process", "cell communication", "Transcriptional misregulation in cancer", "Pathways in cancer", "Interleukin-4 and Interleukin-13 signaling", "cell division", "DNA repair", "Cell Cycle")
highlight_terms <- gostres1$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres1$result <- gostres1$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))
highlight_terms <- gostres2$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres2$result <- gostres2$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))

pdf(file="Run10_degsr_NACT_Naive_targets_partial_GOST.pdf", width=8, height=8)
p1 <- gostplot_MOD(gostres1, capped = FALSE, interactive = FALSE)
p2 <- gostplot_MOD(gostres2, capped = FALSE, interactive = FALSE)
p1 | p2
dev.off()

#### S v R

leftredgenes <- read.table("SvRleftredgenes.txt") %>% pull(V1)
rightgreengenes <- read.table("SvRrightgreengenes.txt") %>% pull(V1)

gostres1 <- gost(query = leftredgenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))
gostres2 <- gost(query = rightgreengenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))

# highlight_terms <- c("response to stimulus", "developmental process", "cell communication", "Transcriptional misregulation in cancer", "Pathways in cancer", "Interleukin-4 and Interleukin-13 signaling", "cell division", "DNA repair", "Cell Cycle")
highlight_terms <- gostres1$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres1$result <- gostres1$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))
highlight_terms <- gostres2$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres2$result <- gostres2$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))

pdf(file="Run10_degsr_Resistant_Sensitive_targets_partial_GOST.pdf", width=16, height=8)
p1 <- gostplot_MOD(gostres1, capped = FALSE, interactive = FALSE)
p2 <- gostplot_MOD(gostres2, capped = FALSE, interactive = FALSE)
p1 | p2
dev.off()

#### C v T

leftredgenes <- read.table("CvTleftredgenes.txt") %>% pull(V1)
rightgreengenes <- read.table("CvTrightgreengenes.txt") %>% pull(V1)

gostres1 <- gost(query = leftredgenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))
gostres2 <- gost(query = rightgreengenes, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))

# highlight_terms <- c("response to stimulus", "developmental process", "cell communication", "Transcriptional misregulation in cancer", "Pathways in cancer", "Interleukin-4 and Interleukin-13 signaling", "cell division", "DNA repair", "Cell Cycle")
highlight_terms <- gostres1$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres1$result <- gostres1$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))
highlight_terms <- gostres2$result %>% group_by(source) %>% slice_head(n = 10) %>% pull(term_name)
gostres2$result <- gostres2$result %>% mutate(term_name = case_when(term_name %in% highlight_terms ~ term_name, TRUE ~ ""))

pdf(file="Run10_degsr_Treated_Control_targets_partial_GOST.pdf", width=16, height=8)
p1 <- gostplot_MOD(gostres1, capped = FALSE, interactive = FALSE)
p2 <- gostplot_MOD(gostres2, capped = FALSE, interactive = FALSE)
p1 | p2
dev.off()


######################### Cutsite Dist ################################

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
gene_bodies <- read.table("/mforge/research/labs/experpath/maia/m237371/OvCancerMULT/features.tsv", sep = "\t") %>% 
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

gene_windows <- read.table("/mforge/research/labs/experpath/maia/m237371/OvCancerMULT/features.tsv", sep = "\t") %>% 
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
