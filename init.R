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


# library(garnett)
# library(scibet)
# library(valr)

# library(clustree)

# library(Polychrome)
# library(ggokabeito)

stat_smooth_func <- function(mapping = NULL, data = NULL,
                        geom = "smooth", position = "identity",
                        ...,
                        method = "auto",
                        formula = y ~ x,
                        se = TRUE,
                        n = 80,
                        span = 0.75,
                        fullrange = FALSE,
                        level = 0.95,
                        method.args = list(),
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE,
                        xpos = NULL,
                        ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}

StatSmoothFunc <- ggproto("StatSmooth", Stat,
                      setup_params = function(data, params) {
                        # Figure out what type of smoothing to do: loess for small datasets,
                        # gam with a cubic regression basis for large data
                        # This is based on the size of the _largest_ group.
                        if (identical(params$method, "auto")) {
                          max_group <- max(table(data$group))
                          
                          if (max_group < 1000) {
                            params$method <- "loess"
                          } else {
                            params$method <- "gam"
                            params$formula <- y ~ s(x, bs = "cs")
                          }
                        }
                        if (identical(params$method, "gam")) {
                          params$method <- mgcv::gam
                        }
                        
                        params
                      },
                      
                      compute_group = function(data, scales, method = "auto", formula = y~x,
                                               se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                               xseq = NULL, level = 0.95, method.args = list(),
                                               na.rm = FALSE, xpos=NULL, ypos=NULL) {
                        if (length(unique(data$x)) < 2) {
                          # Not enough data to perform fit
                          return(data.frame())
                        }
                        
                        if (is.null(data$weight)) data$weight <- 1
                        
                        if (is.null(xseq)) {
                          if (is.integer(data$x)) {
                            if (fullrange) {
                              xseq <- scales$x$dimension()
                            } else {
                              xseq <- sort(unique(data$x))
                            }
                          } else {
                            if (fullrange) {
                              range <- scales$x$dimension()
                            } else {
                              range <- range(data$x, na.rm = TRUE)
                            }
                            xseq <- seq(range[1], range[2], length.out = n)
                          }
                        }
                        # Special case span because it's the most commonly used model argument
                        if (identical(method, "loess")) {
                          method.args$span <- span
                        }
                        
                        if (is.character(method)) method <- match.fun(method)
                        
                        base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                        model <- do.call(method, c(base.args, method.args))
                        
                        m = model
                        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                         list(a = format(coef(m)[[1]], digits = 3), 
                                              b = format(coef(m)[[2]], digits = 3), 
                                              r2 = format(summary(m)$r.squared, digits = 3)))
                        func_string = as.character(as.expression(eq))
                        
                        if(is.null(xpos)) xpos = min(data$x)*0.9
                        if(is.null(ypos)) ypos = max(data$y)*0.9
                        data.frame(x=xpos, y=ypos, label=func_string)
                        
                      },
                      
                      required_aes = c("x", "y")
)

# DittoSeq-v1.4 Colors (based on Okabe-Ito colors)
OkabeIto40 <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#666666",
    "#AD7700", "#1C91D4", "#007756", "#D5C711",
    "#005685", "#A04700", "#B14380", "#4D4D4D",
    "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
    "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C",
    "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E",
    "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3",
    "#8A5F00", "#1674A9", "#005F45", "#AA9F0D",
    "#00446B", "#803800", "#8D3666", "#3D3D3D")
# pdf(file='OkabeIto40.pdf', width=5, height=5)
# show_col(OkabeIto40)
# dev.off()

# https://stackoverflow.com/questions/39694490/highlighting-individual-axis-labels-in-bold-using-ggplot2
colorado <- function(src, boulder) {
  if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
  src_levels <- levels(src)                                 # retrieve the levels in their order
  brave <- boulder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
  if (all(brave)) {                                         # if so
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels)) # then find out where they are
    b_vec <- rep("plain", length(src_levels))               # make'm all plain first
    b_vec[b_pos] <- "bold"                                  # make our targets bold
    b_vec                                                   # return the new vector
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}

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
