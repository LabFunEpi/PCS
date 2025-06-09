args = commandArgs(trailingOnly=TRUE)
wd <- args[1]
obj_rds_filename <- args[2]
ident1 <- args[3] # "Tumor"
ident2 <- args[4] # "FT"
celltype <- args[5] # "customtype"
group <- args[6] # "group"
prefix <- args[7] # "Run1"

setwd(wd)

obj <- readRDS(obj_rds_filename)

statuses <- c(ident1, ident2)

obj$celltype.group <- paste(obj@meta.data[,celltype], obj@meta.data[,group], sep = ".")
Idents(obj) <- "celltype.group"
DefaultAssay(obj) <- "RNA"

types <- unique(obj@meta.data[,celltype])
samples <- obj@meta.data %>% filter(!!as.name(group) %in% statuses) %>% select(sample) %>% unlist() %>% unname() %>% unique()

degs <- lapply(
    types, 
    FUN = function(ct){
        curridents <- unique(Idents(obj))
        if (!(paste0(ct, ".", statuses[[1]]) %in% curridents) | !(paste0(ct, ".", statuses[[2]]) %in% curridents)) return(data.frame())
        if (length(WhichCells(obj, idents = paste0(ct, ".", statuses[[1]]))) < 3 | length(WhichCells(obj, idents = paste0(ct, ".", statuses[[2]]))) < 3) return(data.frame())
        temp1 <- FindMarkers(obj, ident.1 = paste0(ct, ".", statuses[[1]]), ident.2 = paste0(ct, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
        return(temp1)
    }
)
names(degs) <- types
degs <- bind_rows(degs, .id = "celltype")
write.table(degs, file = paste0(prefix, "_degs_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

obj <- DietSeurat(obj, assays = "RNA")
oneouts <- lapply(samples, function(sam){subset(obj, subset = sample != sam)})

degs_robust <- lapply(
    oneouts, 
    FUN = function(oneout){
        curridents <- unique(Idents(oneout))
        ctsp_degs_list <- lapply(
            types,
            FUN = function(ct){
                print(ct)
                if (!(paste0(ct, ".", statuses[[1]]) %in% curridents) | !(paste0(ct, ".", statuses[[2]]) %in% curridents)) return(c())
                if (length(WhichCells(oneout, idents = paste0(ct, ".", statuses[[1]]))) < 3 | length(WhichCells(oneout, idents = paste0(ct, ".", statuses[[2]]))) < 3) return(c())
                temp1 <- FindMarkers(oneout, ident.1 = paste0(ct, ".", statuses[[1]]), ident.2 = paste0(ct, ".", statuses[[2]]), min.pct = 0) %>% rownames_to_column("gene")
                return((temp1 %>% filter(p_val_adj < 0.05))$gene)
            }
        )
        names(ctsp_degs_list) <- types
        return(ctsp_degs_list)
    }
)
names(degs_robust) <- samples
degs_robust <- lapply(degs_robust %>% purrr::transpose(), function(x){Reduce(intersect, x)})

degsr <- bind_rows(lapply(degs_robust, function(gene){data.frame(gene)}), .id = "celltype") %>% left_join(degs)
write.table(degsr, file = paste0(prefix, "_degsr_", statuses[[1]], "_", statuses[[2]], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#######################################################################