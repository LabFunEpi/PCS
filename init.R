library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(cicero)
library(Matrix)
# library(Matrix.utils)
library(GenomeInfoDb)
# library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
library(hdf5r)
library(GenomicRanges)
library(future)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
# library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)

library(ggnewscale)
library(ggpubr)
library(patchwork)
library(scales)
library(cowplot)
library(matrixStats)
library(magrittr)

library(plyr)
library(edgeR)
library(DESeq2)
library(Rsubread)
library(ggrepel)

library(tidyverse)

library(ggsci)
library(gplots)
library(VennDiagram)
library(SeuratDisk)
library(Azimuth)
library(ggrastr)
library(RColorBrewer)
# library(proxyC)

library(fastmatch)
library("data.table")
library(ggbeeswarm)
library(eulerr)
library(coin)

library(glmGamPoi)

library(celldex)
library(SingleR)
library(harmony)
library(ComplexHeatmap)
library(BuenColors)
library(msigdbr)
library(singleseqgset)
library(ggVennDiagram)
library(ggallin)
library(ggbreak)
library(FigR)
library(doParallel)
library(PerformanceAnalytics)
# library(ggpmisc)
library(dendextend)
library(cluster)
library(factoextra)
library(ggforce)
library(clusterCrit)
library(circlize)
library(scico)
library(network)
library(sna)
library(ggnetwork)
library(numbat)
library(gprofiler2)

average.expression <- function(object) {
  data <- GetAssayData(object, slot = "data", assay = "RNA")
  data <- expm1(data)
  idents <- as.vector(Idents(object))
  ident.names <- unique(idents)
  m <- list()
  for (i in 1:length(ident.names)) {
    if(length(which(idents == ident.names[i])) > 1){
      m[[i]] <- Matrix::rowMeans(data[, which(idents == ident.names[i])])
    }else{
      m[[i]] <- data[, which(idents == ident.names[i])]
    }
  }
  result <- do.call(cbind, m)
  colnames(result) <- ident.names
  result <- log1p(result)
  return(result)
}

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
umaps_generic <- function(obj, colors_list, rnatag, atactag, wnntag, ptsize = 0.5){
    p1 <- ggplot(data.frame(obj[[rnatag]][[]], label = factor(Idents(obj))), aes_string(x=paste0(rnatag, "_1"), y=paste0(rnatag, "_2"), color="label")) + 
        ggrastr::rasterise(geom_point(size=ptsize, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = colors_list) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() + theme(legend.title=element_blank())
    p2 <- ggplot(data.frame(obj[[atactag]][[]], label = factor(Idents(obj))), aes_string(x=paste0(atactag, "_1"), y=paste0(atactag, "_2"), color="label")) + 
        ggrastr::rasterise(geom_point(size=ptsize, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = colors_list) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() + theme(legend.position = "none")
    p3 <- ggplot(data.frame(obj[[wnntag]][[]], label = factor(Idents(obj))), aes_string(x=paste0(wnntag, "_1"), y=paste0(wnntag, "_2"), color="label")) + 
        ggrastr::rasterise(geom_point(size=ptsize, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = colors_list) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() + theme(legend.position = "none")
    return(p1 + p2 + p3 + guide_area() + plot_layout(ncol = 2, guides = "collect"))
}

recalculate <- function(seurat_object, rnaassay = "RNA", atacassay = "ATAC"){
    DefaultAssay(seurat_object) <- rnaassay
    seurat_object <- seurat_object %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA(npcs = 50)
    seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:50, reduction.name = "rnaumap", reduction.key = "rnaumap_")
    DefaultAssay(seurat_object) <- atacassay
    seurat_object <- seurat_object %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 20) %>% RunSVD()
    seurat_object <- RunUMAP(seurat_object, reduction = 'lsi', dims = 2:50, reduction.name = "atacumap", reduction.key = "atacumap_")
    seurat_object <- seurat_object %>% FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    seurat_object <- RunUMAP(seurat_object, nn.name = "weighted.nn", reduction.name = "wnnumap", reduction.key = "wnnumap_")
    seurat_object <- seurat_object %>% FindClusters(graph.name = "wsnn", resolution = 0.05, verbose = FALSE)
    return(seurat_object)
}

umaps_wrapped <- function(seurat_object, leftcolors = "hpca_sub", middlecolors = "sample", rightcolors = "seurat_clusters", rnadimreduc = "rnaumap", atacdimreduc = "atacumap", wnndimreduc = "wnnumap", leftcolorpal = NULL, rightcolorpal = NULL, ptsize = 0.5){
    Idents(seurat_object) <- leftcolors
    if (is.null(leftcolorpal)){
        pal <- colors.hpca_sub
    }else{
        pal <- leftcolorpal
    }
    plots1 <- umaps_generic(seurat_object, pal, rnadimreduc, atacdimreduc, wnndimreduc, ptsize = ptsize)
    Idents(seurat_object) <- middlecolors
    plots2 <- umaps_generic(seurat_object, getPalette(length(unique(Idents(seurat_object)))), rnadimreduc, atacdimreduc, wnndimreduc, ptsize = ptsize)
    Idents(seurat_object) <- rightcolors
    if (is.null(rightcolorpal)){
        pal <- getPalette(length(unique(Idents(seurat_object))))
    }else{
        pal <- rightcolorpal
    }
    plots3 <- umaps_generic(seurat_object, pal, rnadimreduc, atacdimreduc, wnndimreduc, ptsize = ptsize)
    return(plots1 | plots2 | plots3)
}

umaps_1dimreduc <- function(obj, colors_list, dimreductag){
    p1 <- ggplot(data.frame(obj[[dimreductag]][[]], label = factor(Idents(obj))), aes_string(x=paste0(dimreductag, "_1"), y=paste0(dimreductag, "_2"), color="label")) + 
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = colors_list) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot()
    return(p1)
}

umaps_wrapped_1dimreduc <- function(seurat_object, dimreductag){
    Idents(seurat_object) <- "hpca_sub"
    plots1 <- umaps_1dimreduc(seurat_object, colors.hpca_sub, dimreductag)
    Idents(seurat_object) <- "sample"
    plots2 <- umaps_1dimreduc(seurat_object, getPalette(length(unique(Idents(seurat_object)))), dimreductag)
    Idents(seurat_object) <- "seurat_clusters"
    plots3 <- umaps_1dimreduc(seurat_object, getPalette(length(unique(Idents(seurat_object)))), dimreductag)
    return(plots1 | plots2 | plots3)
}

##### https://github.com/dgrtwo/drlib/blob/master/R/reorder_within.R
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}
scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}
#####

ExportGroupBW  <- function(
    object,
    assay = NULL,
    group.by = NULL,
    idents = NULL,
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    cutoff = NULL,
    genome = NULL,
    outdir = NULL,
    verbose=TRUE
) {
  # Check if temporary directory exist
  if (!dir.exists(outdir)){
    dir.create(outdir)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) { 
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/") 
    return(NULL) 
  }
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  DefaultAssay(object = object) <- assay
  group.by <- SetIfNull(x = group.by, y = 'ident')
  Idents(object = object) <- group.by
  idents <- SetIfNull(x = idents, y = levels(x = object))
  GroupsNames <- names(x = table(object[[group.by]])[table(object[[group.by]]) > minCells])
  GroupsNames <- GroupsNames[GroupsNames %in% idents]
  # Check if output files already exist
  lapply(X = GroupsNames, FUN = function(x) {
    fn <- paste0(outdir, .Platform$file.sep, x, ".bed")
    if (file.exists(fn)) {
      message(sprintf("The group \"%s\" is already present in the destination folder and will be overwritten !",x))
      file.remove(fn)
    }
  })      
  # Splitting fragments file for each idents in group.by
  SplitFragments(
    object = object,
    assay = assay,
    group.by = group.by,
    idents = idents,
    outdir = outdir,
    file.suffix = "",
    append = TRUE,
    buffer_length = 256L,
    verbose = verbose
  )
  # Column to normalized by
  if(!is.null(x = normMethod)) {
    if (tolower(x = normMethod) %in% c('rc', 'ncells', 'none')){
      normBy <- normMethod
    } else{
      normBy <- object[[normMethod, drop = FALSE]]
    }
  }
  availableChr <- names(x = seqlengths(genome))
  chromLengths <- seqlengths(genome)
  chromSizes <- GRanges(
    seqnames = availableChr,
    ranges = IRanges(
      start = rep(1, length(x = availableChr)),
      end = as.numeric(x = chromLengths)
      )
    )
  if (verbose) {
    message("Creating tiles")
  }
  # Create tiles for each chromosome, from GenomicRanges
  tiles <- unlist(
    x = slidingWindows(x = chromSizes, width = tileSize, step = tileSize)
  )
  if (verbose) {
    message("Creating bigwig files at ", outdir)
  }
  # Run the creation of bigwig for each cellgroups
  if (nbrOfWorkers() > 1) { 
    mylapply <- future_lapply 
  } else { 
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply) 
  }
  
  covFiles <- mylapply(
    GroupsNames,
    FUN = CreateBWGroup,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    tileSize,
    normMethod,
    cutoff,
    outdir
  )
  return(covFiles)
}

CreateBWGroup <- function(
    groupNamei,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    tileSize,
    normMethod,
    cutoff,
    outdir
) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) { 
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/") 
    return(NULL) 
  }
  normMethod <- tolower(x = normMethod)
  # Read the fragments file associated to the group
  fragi <- rtracklayer::import(
    paste0(outdir, .Platform$file.sep, groupNamei, ".bed"), format = "bed"
  )
  cellGroupi <- unique(x = fragi$name)
  # Open the writing bigwig file
  covFile <- file.path(
    outdir,
    paste0(groupNamei, "-TileSize-",tileSize,"-normMethod-",normMethod,".bw")
  )
  
  covList <- lapply(X = seq_along(availableChr), FUN = function(k) {
    fragik <- fragi[seqnames(fragi) == availableChr[k],]
    tilesk <- tiles[BiocGenerics::which(S4Vectors::match(seqnames(tiles), availableChr[k], nomatch = 0) > 0)]
    if (length(x = fragik) == 0) {
      tilesk$reads <- 0
      # If fragments
    } else {
      # N Tiles
      nTiles <- chromLengths[availableChr[k]] / tileSize
      # Add one tile if there is extra bases
      if (nTiles%%1 != 0) {
        nTiles <- trunc(x = nTiles) + 1
      }
      # Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(fragik)$name, cellGroupi)
      
      # For each tiles of this chromosome, create start tile and end tile row,
      # set the associated counts matching with the fragments
      mat <- sparseMatrix(
        i = c(trunc(x = start(x = fragik) / tileSize),
              trunc(x = end(x = fragik) / tileSize)) + 1,
        j = as.vector(x = c(matchID, matchID)),
        x = rep(1, 2*length(x = fragik)),
        dims = c(nTiles, length(x = cellGroupi))
      )
      
      # Max count for a cells in a tile is set to cutoff
      if (!is.null(x = cutoff)){
        mat@x[mat@x > cutoff] <- cutoff
      }
      # Sums the cells
      mat <- rowSums(x = mat)
      tilesk$reads <- mat
      # Normalization
      if (!is.null(x = normMethod)) {
        if (normMethod == "rc") {
          tilesk$reads <- tilesk$reads * 10^4 / length(fragi$name)
        } else if (normMethod == "ncells") {
          tilesk$reads <- tilesk$reads / length(cellGroupi)
        } else if (normMethod == "none") {
        } else {
          if (!is.null(x = normBy)){
            tilesk$reads <- tilesk$reads * 10^4 / sum(normBy[cellGroupi, 1])
          }
        }
      }
    }
    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]
    tilesk
  })
  
  names(covList) <- availableChr
  covList <- as(object = covList, Class = "RleList")
  rtracklayer::export.bw(object = covList, con = covFile)
  return(covFile)
}

###################################################################################
########################### Third Party funciton MODs #############################
###################################################################################

gostplot_MOD <- function(gostres, capped = TRUE, interactive = TRUE, pal = c("GO:MF" = "#dc3912",
                                                                         "GO:BP" = "#ff9900",
                                                                         "GO:CC" = "#109618",
                                                                         "KEGG" = "#dd4477",
                                                                         "REAC" = "#3366cc",
                                                                         "WP" = "#0099c6",
                                                                         "TF" = "#5574a6",
                                                                         "MIRNA" = "#22aa99",
                                                                         "HPA" = "#6633cc",
                                                                         "CORUM" = "#66aa00",
                                                                         "HP" = "#990099")
){
  # gostres is the GOSt response list (contains results and metadata)
  # This function will plot only the sources that were asked from the query

  if( is.null(pal) ){
    pal <- c("GO:MF" = "#dc3912",
             "GO:BP" = "#ff9900",
             "GO:CC" = "#109618",
             "KEGG" = "#dd4477",
             "REAC" = "#3366cc",
             "WP" = "#0099c6",
             "TF" = "#5574a6",
             "MIRNA" = "#22aa99",
             "HPA" = "#6633cc",
             "CORUM" = "#66aa00",
             "HP" = "#990099")
  }

  if (!("result" %in% names(gostres))) stop("Name 'result' not found from the input")
  if (!("meta" %in% names(gostres))) stop("Name 'meta' not found from the input")

  source_order <- logpval <- term_id <- opacity <- NULL
  term_size <- term_name <- p_value <- term_size_scaled <- NULL

  df <- gostres$result
  # Order of data sources comes metadata
  meta <- gostres$meta

  # make sure all the essential column names are there
  essential_names <- c("source_order", "term_size", "term_name", "term_id", "source", "significant")

  if (!(all(essential_names %in% colnames(df)))) stop(paste("The following columns are missing from the result:", paste0(setdiff(essential_names, colnames(df)), collapse = ", ")))

  if (!any(grepl("p_value", colnames(df)))) stop("Column 'p_value(s)' is missing from the result")

  # nr_of_terms for every source
  widthscale <- unlist(lapply(meta$query_metadata$sources, function(x) meta$result_metadata[[x]][["number_of_terms"]]))
  names(widthscale) <- meta$query_metadata$sources # all the sources that were queried

  # Define the start positions for sources in the plot

  # start position for every term
  space <- 1000 # space between different sources
  starts <- c()
  start <- 1
  starts[1] <- start

  if(!length(widthscale) < 2) {
    for(idx in 2:length(widthscale)){
      starts[idx] <- starts[idx - 1] + space + widthscale[idx - 1]
    }
  }

  names(starts) <- names(widthscale)

  # Make sure that all the sources have colors

  if (is.null(names(pal))){
    names(pal) = meta$query_metadata$sources[1:length(pal)]
  }

  sourcediff = setdiff(meta$query_metadata$sources, names(pal))
  colors = grDevices::colors(distinct = TRUE)[grep('gr(a|e)y|white|snow|khaki|lightyellow', grDevices::colors(distinct = TRUE), invert = TRUE)]

  if (length(sourcediff) > 0){
    use_cols = sample(colors, length(sourcediff), replace = FALSE)
    pal[sourcediff] <- use_cols
  }

  # If multiquery
  if("p_values" %in% colnames(df)){
    p_values <- query <- significant <- NULL
    # spread the data frame to correct form
    df$query <- list(names(meta$query_metadata$queries))
    df <- tidyr::unnest(data = df, cols = c(p_values, query, significant))
    df <- dplyr::rename(df, p_value = p_values)
  }

  # Set sizescale of circles
  logScale <- function(input, input_start = 1, input_end = 50000, output_start = 2, output_end = 10){
    m = (output_end - output_start)/(log(input_end) - log(input_start))
    b = -m*log(input_start) + output_start
    output = m * log(input) + b
    return(output)
  }

  # Scale x positions
  xScale <- function(input, input_start = 1, input_end = sum(widthscale) + (length(widthscale) - 1)*space, output_start = 2, output_end = 200){
    m = (output_end - output_start)/(input_end - input_start)
    b = -m*input_start + output_start
    output = m * input + b
    return(output)
  }

  # Add values to df needed for plotting
  # add -log10 pval
  df$logpval <- -log10(df$p_value)
  df$opacity <- ifelse(df$significant, 0.8, ifelse(df$p_value == 1, 0, 0.3))
  df$term_size_scaled = logScale(df$term_size)
  # add x axis position
  df <- df %>% dplyr::group_by(source) %>% dplyr::mutate(order = xScale(source_order, input_start = 1, input_end = widthscale[source], output_start = starts[source], output_end = starts[source] + widthscale[source]))
  df$order <- xScale(df$order)

  if (capped) {
    df$logpval[df$logpval > 16] <- 17
    ymin <- -1
    ymax <- 18.5
    ticklabels <- c("0", "2", "4", "6", "8", "10", "12", "14", ">16")
    tickvals <- c(0, 2, 4, 6, 8, 10, 12, 14, 16)
  } else {
    ymin <- -1
    ymax <- ceiling(max(df$logpval)) + 5
    ticklabels <- ggplot2::waiver()
    tickvals <- ggplot2::waiver()
  }

  if (interactive){
    # Start plotting
    sd <- crosstalk::SharedData$new(df, key = ~term_id)
  } else {
    sd <- df
  }

  p <- ggplot2::ggplot(data = sd, ggplot2::aes(x = order, y = logpval, text = paste(term_id, paste0('(', term_size,')'), '<br>', term_name, '<br>', formatC(p_value, format = "e", digits = 3)))) +
    ggplot2::geom_point(ggplot2::aes(color = source, size = term_size_scaled, alpha = opacity), show.legend = FALSE) +
    geom_text_repel(aes(label = term_name), min.segment.length = unit(0, 'lines'), show.legend = FALSE, max.overlaps = 10, force = 5) +
	ggplot2::facet_wrap(~ query, ncol = 1, scales = "free_x", shrink = FALSE) +
    ggplot2::ylab("-log10(p-adj)") +
    theme_cowplot() +
    ggplot2::theme(legend.position='none',
                   panel.border=ggplot2::element_blank(),
                   strip.text=ggplot2::element_text(size=12, colour = "darkgrey"),
                   strip.background=ggplot2::element_blank(),
                   axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_text(size = 8, angle=45, hjust = 1),
                   axis.ticks.x=ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_text(angle = 0, hjust = 0, size = 10),
                   plot.margin = ggplot2::margin(t = 0, r = 5, b = 20, l = 20, unit = "pt"),
                   axis.title.y = ggplot2::element_text(size = 10, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))
    ) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_alpha(range = c(0, 0.8), limits = c(0, 0.8)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax),
                                labels = ticklabels,
                                breaks = tickvals) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 210),
                                breaks = (xScale(starts) + xScale(starts + widthscale))/2,
                                labels = names(widthscale))

  for (s in names(widthscale)){
    xstart = xScale(starts[s])
    xend = xScale(starts[s] + widthscale[s])
    p <- p + ggplot2::annotate("segment", x = xstart, xend = xend, y = -1, yend = -1,
                               size = 3, colour = pal[s])
  }

  if (capped){
    p <- p + ggplot2::annotate(geom = "text", x = 180,
                               y = 16.2, label = "values above this threshold are capped", size = 2, color = "grey") +
      ggplot2::geom_hline(yintercept = 16, linetype = "dashed", size = 0.2, color = "grey")
  }

  if (interactive){
    p <- p %>% plotly::ggplotly(tooltip = "text")
    p <- p %>% plotly::highlight(on = "plotly_click", off = "plotly_doubleclick", dynamic = FALSE, persistent = FALSE)
  }

  return(p)
}


SetIfNull <- function(x, y) {
    if (is.null(x = x)) {
        return(y)
    } else {
        return(x)
    }
}

SingleDimPlot <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  shape.by = NULL,
  alpha = 1,
  alpha.by = NULL,
  order = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  raster = NULL,
  raster.dpi = NULL
) {
  if ((nrow(x = data) > 1e5) & is.null(x = raster)){
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data) > 1e5)
  pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)

  if (!is.null(x = cells.highlight) && pt.size != AutoPointSize(data = data, raster = raster) && sizes.highlight != pt.size && isTRUE(x = raster)) {
    warning("When `raster = TRUE` highlighted and non-highlighted cells must be the same size. Plot will use the value provided to 'sizes.highlight'.")
  }

  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2)
      stop("'raster.dpi' must be a two-length numeric vector")
  }
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    if (inherits(x = cells.highlight, what = "data.frame")) {
      stop("cells.highlight cannot be a dataframe. ",
           "Please supply a vector or list")
    }
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight %||% pt.size,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% '#C3C3C3',
      pt.size = pt.size,
      raster = raster
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(!is.na(x = data[, col.by]), data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning(
      "Cannot find alpha variable ",
      alpha.by,
      " in data, setting to NULL",
      call. = FALSE,
      immediate. = TRUE
    )
    alpha.by <- NULL
  }

  plot <- ggplot(data = data)
  plot <- if (isTRUE(x = raster)) {
    plot + geom_scattermore(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      pointsize = pt.size,
      alpha = alpha,
      pixels = raster.dpi
    )
  } else {
    plot + ggrastr::rasterise(geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by,
        alpha = alpha.by
      ),
      size = pt.size,
      alpha = alpha
    ), dpi = 400)
  }
  plot <- plot +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(color = NULL, title = col.by) +
    CenterTitle()
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

FeaturePlot <- function(
  object,
  features,
  dims = c(1, 2),
  cells = NULL,
  cols = if (blend) {
    c('lightgrey', '#ff0000', '#00ff00')
  } else {
    c('lightgrey', 'blue')
  },
  pt.size = NULL,
  alpha = 1,
  order = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  keep.scale = "feature",
  shape.by = NULL,
  slot = 'data',
  blend = FALSE,
  blend.threshold = 0.5,
  label = FALSE,
  label.size = 4,
  label.color = "black",
  repel = FALSE,
  ncol = NULL,
  coord.fixed = FALSE,
  by.col = TRUE,
  sort.cell = deprecated(),
  interactive = FALSE,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
) {
  if (isTRUE(x = interactive)) {
    return(IFeaturePlot(
      object = object,
      feature = features[1],
      dims = dims,
      reduction = reduction,
      slot = slot
    ))
  }
  # Set a theme to remove right-hand Y axis lines
  # Also sets right-hand Y axis text label formatting
  no.right <- theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 7)
    )
  )
  # Get the DimReduc to use
  reduction <- reduction %||% DefaultDimReduc(object = object)
  # Figure out blending stuff
  if (isTRUE(x = blend) && length(x = features) != 2) {
    abort(message = "Blending feature plots only works with two features")
  }
  # Set color scheme for blended FeaturePlots
  if (isTRUE(x = blend)) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(
      EXPR = as.character(x = length(x = cols)),
      '0' = {
        warn(message = "No colors provided, using default colors")
        default.colors
      },
      '1' = {
        warn(message = paste(
          "Only one color provided, assuming",
          sQuote(x = cols),
          "is double-negative and augmenting with default colors"
        ))
        c(cols, default.colors[2:3])
      },
      '2' = {
        warn(message = paste(
          "Only two colors provided, assuming specified are for features and agumenting with",
          sQuote(default.colors[1]),
          "for double-negatives",
        ))
        c(default.colors[1], cols)
      },
      '3' = cols,
      {
        warn(message = "More than three colors provided, using only first three")
        cols[1:3]
      }
    )
  }
  if (isTRUE(x = blend) && length(x = cols) != 3) {
    abort("Blending feature plots only works with three colors; first one for negative cells")
  }
  # Name the reductions
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% Cells(x = object[[reduction]])
  # Get plotting data
  data <- FetchData(
    object = object,
    vars = c(dims, 'ident', features),
    cells = cells,
    slot = slot
  )
  # Check presence of features/dimensions
  if (ncol(x = data) < 4) {
    abort(message = paste(
      "None of the requested features were found:",
      paste(features, collapse = ', '),
      "in slot ",
      slot
    ))
  } else if (!all(dims %in% colnames(x = data))) {
    abort(message = "The dimensions requested were not found")
  }
  features <- setdiff(x = names(x = data), y = c(dims, 'ident'))
  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    abort(
      message = "There must be the same number of minimum and maximum cuttoffs as there are features"
    )
  }
  names(x = min.cutoff) <- names(x = max.cutoff) <- features
  brewer.gran <- ifelse(
    test = length(x = cols) == 1,
    yes = brewer.pal.info[cols, ]$maxcolors,
    no = length(x = cols)
  )
  # Apply cutoffs
  for (i in seq_along(along.with = features)) {
    f <- features[i]
    data.feature <- data[[f]]
    min.use <- SetQuantile(cutoff = min.cutoff[f], data = data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[f], data = data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    if (brewer.gran != 2) {
      data.feature <- if (all(data.feature == 0)) {
        rep_len(x = 0, length.out = length(x = data.feature))
      } else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = 2
        )))
      }
    }
    data[[f]] <- data.feature
  }
  # Figure out splits (FeatureHeatmap)
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(
      EXPR = split.by,
      ident = Idents(object = object)[cells, drop = TRUE],
      object[[split.by, drop = TRUE]][cells, drop = TRUE]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  # Set shaping variable
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  # Make list of plots
  plots <- vector(
    mode = "list",
    length = ifelse(
      test = blend,
      yes = 4,
      no = length(x = features) * length(x = levels(x = data$split))
    )
  )
  # Apply common limits
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  # Set blended colors
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(
      two.colors = cols[2:3],
      col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(
      color.matrix[, 1],
      color.matrix[1, ],
      as.vector(x = color.matrix)
    )
  }
  # Make the plots
  for (i in 1:length(x = levels(x = data$split))) {
    # Figure out which split we're working with
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, , drop = FALSE]
    # Blend expression values
    if (isTRUE(x = blend)) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[, features]) == 0]
      if (length(x = no.expression) != 0) {
        abort(message = paste(
          "The following features have no value:",
          paste(no.expression, collapse = ', ')
        ))
      }
      data.plot <- cbind(data.plot[, c(dims, 'ident')], BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    # Make per-feature plots
    for (j in 1:length(x = features)) {
      feature <- features[j]
      # Get blended colors
      if (isTRUE(x = blend)) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, 'ident', feature, shape.by)]
      # Make the plot
      plot <- SingleDimPlot(
        data = data.single,
        dims = dims,
        col.by = feature,
        order = order,
        pt.size = pt.size,
        alpha = alpha,
        cols = cols.use,
        shape.by = shape.by,
        label = FALSE,
        raster = raster,
        raster.dpi = raster.dpi
      ) +
        scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) +
        theme_cowplot() +
        CenterTitle()
        # theme(plot.title = element_text(hjust = 0.5))
      # Add labels
      if (isTRUE(x = label)) {
        plot <- LabelClusters(
          plot = plot,
          id = 'ident',
          repel = repel,
          size = label.size,
          color = label.color
        )
      }
      # Make FeatureHeatmaps look nice(ish)
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, colour = 'black'))
        # Add title
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        # Add second axis
        if (j == length(x = features) && !blend) {
          suppressMessages(
            expr = plot <- plot +
              scale_y_continuous(
                sec.axis = dup_axis(name = ident),
                limits = ylims
              ) +
              no.right
          )
        }
        # Remove left Y axis
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        # Remove bottom X axis
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      # Add colors scale for normal FeaturePlots
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warn(message = paste0(
              "All cells have the same value (",
              unique.feature.exp,
              ") of ",
              dQuote(x = feature)
            ))
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else{
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(
            expr = plot + scale_color_gradientn(
              colors = cols.grad,
              guide = "colorbar"
            )
          )
        }
      }
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" && !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)))
      }
      # Add coord_fixed
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      # I'm not sure why, but sometimes the damn thing fails without this
      # Thanks ggplot2
      plot <- plot
      # Place the plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  # Add blended color key
  if (isTRUE(x = blend)) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(
          blend.legend +
            scale_y_continuous(
              sec.axis = dup_axis(name = ifelse(
                test = length(x = levels(x = data$split)) > 1,
                yes = levels(x = data$split)[ii],
                no = ''
              )),
              expand = c(0, 0)
            ) +
            labs(
              x = features[1],
              y = features[2],
              title = if (ii == 1) {
                paste('Color threshold:', blend.threshold)
              } else {
                NULL
              }
            ) +
            no.right
        ),
        after = 4 * ii - 1
      ))
    }
  }
  # Remove NULL plots
  plots <- Filter(f = Negate(f = is.null), x = plots)
  # Combine the plots
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(
    test = is.null(x = split.by) || isTRUE(x = blend),
    yes = ncol,
    no = length(x = features)
  )
  legend <- if (isTRUE(x = blend)) {
    'none'
  } else {
    split.by %iff% 'none'
  }
  # Transpose the FeatureHeatmap matrix (not applicable for blended FeaturePlots)
  if (isTRUE(x = combine)) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(
        X = plots,
        FUN = function(x) {
          return(suppressMessages(
            expr = x +
              theme_cowplot() +
              ggtitle("") +
              scale_y_continuous(sec.axis = dup_axis(name = ""), limits = ylims) +
              no.right
          ))
        }
      )
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(
          expr = plots[[i]] +
            scale_y_continuous(
              sec.axis = dup_axis(name = features[[idx]]),
              limits = ylims
            ) +
            no.right
        )
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) == 1)) {
        plots[[i]] <- plots[[i]] +
          ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] +
            ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(
        what = rbind,
        args = split(
          x = 1:length(x = plots),
          f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features))
        )
      ))]
      # Set ncol to number of splits (nrow) and nrow to number of features (ncol)
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == 'none') {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == 'none') {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" && !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      plots <- suppressMessages(plots & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)))
    }
  }
  return(plots)
}

###################################################################################
########################### Gene Sets #############################################
###################################################################################

Verhaak_BadProg <- c("OR3A1", "HIF3A", "MARK4", "CCDC49", "PLXNA1", "MED25", "DCP2", "TNFSF11", "KPNA3", "B4GALT5", "CTNNBL1", "PAK4", "MUC4", "GAGE1", "CALB1", "MLL4", "DLGAP4", "SUPT5H", "FLJ10241", "DBF4B", "ACTR5", "CLTCL1", "KRT4", "RP11-125A7.3", "ANXA4", "GRB7", "SASH1", "PERLD1", "MTRF1", "DHX34", "PHF20", "PPM2C", "APC", "ZFP36", "EFNA1", "TTC31", "ERBB2", "CDH16", "RAB4B", "ZNF611", "GTF2F2", "ZNF331", "CASC3", "IRF2BP1", "RGS19", "ELMO2", "XRCC4", "C19ORF7", "WIPF2", "LTBP1", "RNF44", "RPL23", "TGIF2", "MAPT", "DMWD", "UTP6", "NEK3", "GMEB2", "SPAG5", "PTPN1", "ATP1A3", "BLCAP", "DHX35", "STRN4", "HIGD1B", "USH1C", "TBP", "TRIM13", "ITPKC", "KRT13", "ADA", "CEMP1", "SLC25A30", "HHLA2", "SAMD4B", "YTHDC2", "ZHX3", "AKT2", "GEMIN7", "SH3TC2", "STK4", "ANKRD46", "PSMC4", "HSPBAP1", "NUFIP1", "TCF15", "C20ORF4", "GRWD1", "FLJ20323", "C20ORF3", "SMG5", "C17ORF63", "KCNG1", "LRRC31", "DMPK", "CCNF", "ACVR2A", "FBXL18", "CAMK2G", "EML2", "MGRN1", "ZNF574", "MRPS31", "NLK", "C20ORF121", "PABPC1", "PLAC4", "CLK2")

Verhaak_GoodProg <- c("ALDH5A1", "GJB1", "AADAC", "DMC1", "GEMIN8", "ZFR", "HSPA1L", "PIGP", "DAP", "IL2", "PCK2", "GMPR", "HIST1H2BE", "CXCL9", "CECR1", "GSTZ1", "ETV4", "THEM2", "HLA-DOB", "RBP1", "HIRIP3", "FJX1", "FAH", "CD27", "HLA-DMB", "RTN3", "SHMT2", "CUTA", "ARL6IP5", "CREB3", "FLOT1", "WARS", "TXNDC13", "PHKG2", "CLEC11A", "CXCL13", "HIST1H2BO", "EDNRB", "HIST1H2AM", "BCL7C", "HIST1H2AB", "MTHFS", "SLC1A4", "HIST1H3E", "MAST4", "TRIM27", "GAD1", "TCP11", "SLC7A11", "FAM8A1", "HLCS", "UBD", "MBNL3", "ICOS", "BTN2A2", "ZNF76", "BSCL2", "TMEM30B", "WWOX", "HYAL2", "ID4", "EFTUD1", "BMPR2", "STK16", "PHLDA3", "GALNT6", "DNAJC4", "MATN2", "PNPLA4", "CUL4B", "HIST1H4H", "LEF1", "HLA-DOA", "SPRY2", "BTN3A1", "C14ORF159", "HIST1H2BI", "ZNF184", "NAB1", "C6ORF64", "FLJ14213", "NDUFC2", "CSRP1", "COMMD4", "NOTCH4")

Verhaak_differentiated <- c("FUT2", "STEAP3", "LCN2", "HRBL", "MYO5C", "ALS2CL", "ANXA1", "PTGER2", "MVP", "SSH3", "C1ORF116", "DTX4", "MLPH", "TNFRSF14", "MGLL", "CLU", "TMEM24", "SLC37A1", "RIN1", "SCGB2A2")
Verhaak_immunoreactive <- c("FBXO5", "LAG3", "AIM2", "CXCL10", "IFI30", "CASP1", "ECGF1", "LAIR1", "IL21R", "SLC31A2", "PDZK1IP1", "SAMSN1", "FCER1G", "RARRES3", "ADAMDEC1", "SLA", "LAPTM5", "UBE2L6", "KMO", "CIITA", "GIMAP5", "APOBEC3G")
Verhaak_mesenchymal <- c("PCOLCE", "MXRA8", "MAPRE2", "CXCR7", "MMP14", "PDPN", "PLEKHO1", "TPST1", "MARCKS", "CALD1", "TMEPAI", "NNMT", "GAS7", "WIPF1", "FSCN1", "NBL1", "LAMB1", "CACNA1C", "DLC1", "RUNX1", "KIAA0247", "PALLD", "RHOBTB3", "NUAK1", "MYO1B", "COL5A1", "COL8A1", "SRPX2", "LGALS1", "CTGF", "DDR2", "EVI2A", "MCAM", "FN1", "KDELC1", "STCH", "F2R")
Verhaak_proliferative <- c("COPS3", "EFS", "MARCKSL1", "RAD51AP1", "PKIA", "SALL2", "NETO2", "SPATS2", "MAPRE1", "UCHL1", "MFAP2", "TCF7L1", "LRP4", "BEX1", "KIF1A", "COL4A6", "TRO", "BLMH", "STMN1", "SMARCD1", "C5ORF13")

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

PI3K <- c("ACACA", "ACTR2", "ACTR3", "ADCY2", "GRK2", "AKT1", "AKT1S1", "AP2M1", "ARF1", "ARHGDIA", "ARPC3", "ATF1", "CAB39", "CAB39L", "CALR", "CAMK4", "CDK1", "CDK2", "CDK4", "CDKN1A", "CDKN1B", "CFL1", "CLTC", "CSNK2B", "CXCR4", "DAPP1", "DDIT3", "DUSP3", "E2F1", "ECSIT", "EGFR", "EIF4E", "FASLG", "FGF17", "FGF22", "FGF6", "GNA14", "GNGT1", "GRB2", "GSK3B", "HRAS", "HSP90B1", "IL2RG", "IL4", "IRAK4", "ITPR2", "LCK", "MAP2K3", "MAP2K6", "MAP3K7", "MAPK1", "MAPK10", "MAPK8", "MAPK9", "MAPKAP1", "MKNK1", "MKNK2", "MYD88", "NCK1", "NFKBIB", "NGF", "NOD1", "PAK4", "PDK1", "PFN1", "PIK3R3", "PIKFYVE", "PIN1", "PITX2", "PLA2G12A", "PLCB1", "PLCG1", "PPP1CA", "PPP2R1B", "PRKAA2", "PRKAG1", "PRKAR2A", "PRKCB", "PTEN", "PTPN11", "RAC1", "RAF1", "RALB", "RIPK1", "RIT1", "RPS6KA1", "RPS6KA3", "RPTOR", "SFN", "SLA", "SLC2A1", "SMAD2", "SQSTM1", "STAT2", "TBK1", "THEM4", "TIAM1", "TNFRSF1A", "TRAF2", "TRIB3", "TSC2", "UBE2D3", "UBE2N", "VAV3", "YWHAB")

Notch <- c("APH1A", "ARRB1", "CCND1", "CUL1", "DLL1", "DTX1", "DTX2", "DTX4", "FBXW11", "FZD1", "FZD5", "FZD7", "HES1", "HEYL", "JAG1", "KAT2A", "LFNG", "MAML2", "NOTCH1", "NOTCH2", "NOTCH3", "PPARD", "PRKCA", "PSEN2", "PSENEN", "RBX1", "SAP30", "SKP1", "ST3GAL6", "TCF7L2", "WNT2", "WNT5A")

## Lahtinen et al. (Supplementary Tables)
overrep_cols = c("ID_Reactome_pathway", "Reactome_pathway", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "gene_list")
Adaptive = read.table("Adaptive.txt", sep = "\t") %>% set_colnames(overrep_cols)
Adaptive = lapply(Adaptive$gene_list, function(gl){str_split_1(gl, ",")}) %>% unlist() %>% unique() %>% sort()

Evolving = read.table("Evolving.txt", sep = "\t") %>% set_colnames(overrep_cols)
Evolving = lapply(Evolving$gene_list, function(gl){str_split_1(gl, ",")}) %>% unlist() %>% unique() %>% sort()

Maintaining = read.table("Maintaining.txt", sep = "\t") %>% set_colnames(overrep_cols)
Maintaining = lapply(Maintaining$gene_list, function(gl){str_split_1(gl, ",")}) %>% unlist() %>% unique() %>% sort()

postNACT = read.table("postNACT.txt", sep = "\t") %>% set_colnames(overrep_cols)
postNACT = lapply(postNACT$gene_list, function(gl){str_split_1(gl, ",")}) %>% unlist() %>% unique() %>% sort()

Relapse = read.table("Relapse.txt", sep = "\t") %>% set_colnames(overrep_cols)
Relapse = lapply(Relapse$gene_list, function(gl){str_split_1(gl, ",")}) %>% unlist() %>% unique() %>% sort()
