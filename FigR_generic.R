
args = commandArgs(trailingOnly=TRUE)
wd <- args[1]
obj_rds_filename <- args[2]
prefix <- args[3]
assay <- args[4]
genome <- args[5]
figrStep <- as.numeric(args[6])

setwd(wd)
print(paste0(Sys.time(), ":: Reading input object ... "))
obj <- readRDS(obj_rds_filename)

library(FigR)
library(doParallel)
library(fastmatch)

ATAC_se <- SummarizedExperiment(assays=SimpleList(counts=obj@assays[[assay]]@counts), rowRanges=granges(obj[[assay]]), colData=obj@meta.data)
RNA_sm <- obj[["RNA"]]$counts

if (figrStep == 1){
    print(paste0(Sys.time(), ":: Running runGenePeakcorr() ... "))
    cisCor <- runGenePeakcorr(ATAC.se = ATAC_se, RNAmat = RNA_sm, genome = genome, nCores = 14, p.cut=NULL, keepPosCorOnly=FALSE, windowPadSize = 150000)
    saveRDS(cisCor, paste0(prefix, "_cisCor_relaxed.rds"))
    print(paste0(Sys.time(), ":: Finished runGenePeakcorr() !!!"))
}

if (figrStep == 2){
    print(paste0(Sys.time(), ":: Building FigR GRN pre-requisites ... "))
    cisCor <- readRDS(paste0(prefix, "_cisCor_relaxed.rds"))
    cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)
    
    dorcGenes <- cisCor.filt %>% dorcJPlot(cutoff=2, returnGeneList = TRUE)
    dorcMat <- getDORCScores(ATAC_se,dorcTab=cisCor.filt,geneList=dorcGenes,nCores=14)
    cellkNN <- FNN::get.knn(obj[["pca"]][[]],k=50)$nn.index
    rownames(cellkNN) <- colnames(dorcMat)
    
    print(paste0(Sys.time(), ":: Smoothing dorcMat ... "))
    dorcMat.smooth <- smoothScoresNN(NNmat=cellkNN,mat=dorcMat,nCores=14)
    print(paste0(Sys.time(), ":: Smoothing RNAmat ... "))
    RNAmat.smooth <- smoothScoresNN(NNmat=cellkNN,mat=RNA_sm,nCores=14)
    print(paste0(Sys.time(), ":: Saving figR_GRN.RData ... "))
    save(cellkNN, dorcMat.smooth, RNAmat.smooth, file = paste0(prefix, "_GRN_relaxed.RData"))
    print(paste0(Sys.time(), ":: Saved figR_GRN.RData ... "))
}

if (figrStep == 3){
    load(paste0(prefix, "_GRN_relaxed.RData"))
    cisCor <- readRDS(paste0(prefix, "_cisCor_relaxed.rds"))
    cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)

    print(paste0(Sys.time(), ":: Building FigR GRN (Motif) ... "))
    figR.d <- runFigRGRN(ATAC.se = ATAC_se, dorcTab = cisCor.filt, genome = genome, dorcMat = dorcMat.smooth, dorcK = 150, rnaMat = RNAmat.smooth, nCores = 14)
    write.table(figR.d, file = paste0(prefix, "_figR_motif_relaxed.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    print(paste0(Sys.time(), ":: Finished FigR GRN (Motif) !!!"))
}

# if (figrStep == 4){
    # source("/mforge/research/labs/experpath/maia/m237371/Generics/runFigRGRN_TFBS.R")
    # load(paste0(prefix, "_figR_GRN.RData"))
    # cisCor <- readRDS(paste0(prefix, "_cisCor.rds"))
    # cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)
    
    # print(paste0(Sys.time(), ":: Reading remap nr (big)data ... "))
    # if (genome == "hg38"){
        # remapPath <- "/mforge/research/labs/experpath/maia/m237371/public_data/remap2022_nr_macs2_hg38_v1_0.bed"
    # }else{
        # remapPath <- "/mforge/research/labs/experpath/maia/m237371/public_data/remap2022_nr_macs2_mm10_v1_0.bed"
    # }
    # TFBS <- read.table(remapPath) %>% 
        # select(V1:V4) %>% 
        # set_colnames(c("chr", "start", "end", "name")) %>%
        # separate(col = "name", sep = ":", into = c("TF", NA))

    # print(paste0(Sys.time(), ":: Building FigR GRN (TFBS) ... "))
    # figR.d <- runFigRGRN_TFBS(ATAC.se = ATAC_se, dorcTab = cisCor.filt, genome = genome, dorcMat = dorcMat.smooth, dorcK = 150, rnaMat = RNAmat.smooth, nCores = 14, TFBS = TFBS)
    # write.table(figR.d, file = paste0(prefix, "_figR_TFBS.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# }
