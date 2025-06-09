
args = commandArgs(trailingOnly=TRUE)
wd <- args[1]
filtered_rds_filename <- args[2]

setwd(wd)

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
# pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
pwm_set <- readRDS("/mforge/research/labs/experpath/maia/m237371/public_data/cisBP_human_pfms_2021.rds")
library(BiocParallel)
register(SerialParam())

filtered <- readRDS(filtered_rds_filename)

###### Motif ######

# DefaultAssay(filtered) <- "mATAC"
# filtered <- AddMotifs(filtered, genome = genome, pfm = pwm_set)
# filtered <- RunChromVAR(filtered, genome = genome, assay = "mATAC")

# saveRDS(filtered, filtered_rds_filename)

###### Remap ######

remap <- read.table("/mforge/research/labs/experpath/maia/m237371/public_data/myremap.bed") %>% 
    select(V1:V4) %>% 
    set_colnames(c("chr", "start", "end", "TF"))

DefaultAssay(filtered) <- "mATAC"
MULT_peaks <- data.frame(Peak = rownames(filtered)) %>% separate(col = "Peak", sep = "-", into = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()

TFs <- unique(remap$TF)
TF_features <- lapply(TFs, function(cur_TF){
    cur_peaks <- makeGRangesFromDataFrame(remap %>% filter(TF == cur_TF))
    cur_peaks <- cur_peaks[as.vector(seqnames(cur_peaks) %in% standardChromosomes(cur_peaks))]
    ovp_features <- suppressWarnings(GRangesToString(MULT_peaks[case_when(is.na(findOverlaps(MULT_peaks, cur_peaks, select = "first", ignore.strand=TRUE)) ~ FALSE, TRUE ~ TRUE)]))
    return(ovp_features)
})
names(TF_features) <- paste0("remap_", TFs)
remap <- NULL
filtered <- AddChromatinModule(filtered, features = TF_features, genome = genome, assay = "mATAC")

remapscores <- filtered@meta.data %>% select(starts_with("remap-")) %>% t()
rownames(remapscores) <- str_sub(rownames(remapscores), 7L, -1L)
filtered[["remap"]] <- CreateAssayObject(data = remapscores)
filtered@meta.data <- filtered@meta.data %>% select(!starts_with("remap-"))

saveRDS(filtered, filtered_rds_filename)

