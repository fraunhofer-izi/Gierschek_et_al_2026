
sc_gating = function(obj, obj.l, model) {
  model_name <- if (isS4(model)) model@name else model$name
  output_name <- paste0("SCGATE_", model_name)

  df_list = parallel::mclapply(obj.l, function(x) {
    x = scGate(
      x, model = model, assay = "RNA", slot = "data",
      output.col.name = output_name, ncores = 10
    )

    matched_cols <- grep(paste0("^", output_name), colnames(x@meta.data), value = TRUE)

    if (length(matched_cols) == 0) {
      stop(paste("No columns matching", output_name, "found in metadata."))
    }

    df <- x@meta.data[, matched_cols, drop = FALSE]
    df$barcode <- rownames(df)
    return(df)
  }, mc.cores = 10)

  df = do.call("rbind", df_list)
  rownames(df) = df$barcode
  df$barcode = NULL
  obj = AddMetaData(obj, df)
  
  return(obj)
}


process_cellranger_metrics <- function(seurat_path, batch_name) {
  cellranger.dirs <- list.dirs(path = seurat_path, full.names = TRUE, recursive = FALSE)
  cellranger.samples <- basename(cellranger.dirs)

  fltrd.dirs <- paste0(cellranger.dirs, "/metrics_summary.csv")
  names(fltrd.dirs) <- cellranger.samples

  metrics_list <- lapply(fltrd.dirs, function(f) {
    if (file.exists(f)) {
      df <- read.csv(f, stringsAsFactors = FALSE)
      df$Sample <- basename(dirname(dirname(f)))
      return(df)
    } else {
      warning(paste("File not found:", f))
      return(NULL)
    }
  })

  metrics_summary <- do.call(rbind, metrics_list)

  required_metrics <- c(
    "Cells",
    "Mean reads per cell",
    "Median UMI counts per cell",
    "Median genes per cell"
  )

  metrics_filtered <- metrics_summary %>%
    filter(Metric.Name %in% required_metrics) %>%
    filter(Category != "Library") %>%
    mutate(
      Sample_ID = rownames(.),
      Metric.Value = as.numeric(gsub(",", "", Metric.Value)),
      Sample_Clean = sub("\\.\\d+$", "", Sample_ID),
      Batch = batch_name  # <-- Add batch column here
    ) %>%
    group_by(Sample_Clean, Metric.Name, Batch) %>%
    summarise(Mean_Value = mean(Metric.Value, na.rm = TRUE), .groups = "drop")

  return(metrics_filtered)
}



#' Split a Seurat Object into a List of Subsets Based on Metadata
#'
#' This function splits a `Seurat` object into a list of smaller `Seurat` objects based on a grouping variable 
#' (e.g., sample ID, condition). Splitting is done in parallel using `BiocParallel`.
#'
#' @param object A `Seurat` object to be split.
#' @param split.by Character. Name of the metadata column to split by (default: `"orig.ident"`).
#' @param threads Integer. Number of parallel threads to use (default: 5). If `NULL`, uses one thread per group.
#'
#' @return A named list of `Seurat` objects, each corresponding to a level of the grouping variable.
#'
#' @examples
#' split_list <- Split_Object(seurat_object, split.by = "sample", threads = 4)
#'
#' @importFrom Seurat FetchData subset
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
Split_Object = function(object, split.by = "orig.ident", threads = 5) {

  groupings <- FetchData(object = object, vars = split.by)[, 1]
  groupings <- unique(x = as.character(x = groupings))
  names(groupings) = groupings

  if (is.null(threads)) {
    bpparam = BiocParallel::MulticoreParam(workers = length(groupings))
  } else {
    bpparam = BiocParallel::MulticoreParam(workers = threads)
  }

  obj.list = BiocParallel::bplapply(groupings, function(grp) {
    cells <- which(x = object[[split.by, drop = TRUE]] == grp)
    cells <- colnames(x = object)[cells]
    se = subset(x = object, cells = cells)
    se@meta.data = droplevels(se@meta.data)
    se
  }, BPPARAM = bpparam)

  return(obj.list)
}

#' Generate a QC Violin Plot of Features per Cell
#'
#' This function creates a violin plot showing the distribution of a specified feature 
#' (typically number of detected genes or UMIs) per cell across different groups 
#' (e.g., samples or conditions). Useful for quality control in single-cell RNA-seq analysis.
#'
#' @param obj A `Seurat` object, a list of `Seurat` objects, or a `data.frame` containing cell-level metadata.
#' @param .features Character. Name of the feature to plot (e.g., `"nFeature_RNA"`).
#' @param .group.by Character. Metadata column to group cells by on the x-axis (e.g., `"orig.ident"`).
#' @param plot_title Character. Title of the plot.
#' @param x_axis_label Character or NULL. Label for the x-axis. If NULL, uses `.group.by`.
#' @param y_axis_label Character. Label for the y-axis (default: `"Features"`).
#' @param low_cutoff Numeric or NULL. Optional lower horizontal cutoff line (e.g., gene/cell threshold).
#' @param high_cutoff Numeric or NULL. Optional upper horizontal cutoff line.
#'
#' @return A `ggplot` object displaying a violin plot of feature distribution per group.
#'
#' @examples
#' qc_vln_plot_cell(obj = seurat_object, .features = "nFeature_RNA", .group.by = "orig.ident")
#'
#' @import ggplot2
#' @importFrom dplyr select
#' @export
qc_vln_plot_cell = function(
    obj = se.meta,
    .features = "nFeature_RNA",
    .group.by = "orig.ident",
    plot_title = "Genes Per Cell",
    x_axis_label = NULL,
    y_axis_label = "Features",
    low_cutoff = NULL,
    high_cutoff = NULL
){
  
  if(!is.list(obj)) {
    df = obj@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]])
    df
  }
  if(class(obj) == "list") {
    l = lapply(obj, function(x){
      df = x@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]])
    })
    df = do.call("rbind", l)
    df
  }
  if(class(obj) == "data.frame") {
    df = obj %>% dplyr::select(.data[[.group.by]], .data[[.features]])
    df
  }
  
  ggplot(data = df, mapping = aes(x = .data[[.group.by]], y = .data[[.features]])) +
    geom_violin(
      size = .1,
      width = 1,
      scale = "area",
      na.rm = TRUE
    ) +
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, colour = "#0077BB", size = .2) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = rel(1.2))
    ) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "#BB5566", size = .3) +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title)
}

#' Count cells per sample in a list of Seurat objects
#'
#' This function counts the number of cells per sample (based on 'orig.ident') in a list of Seurat objects.
#' Optionally, it can add the counts as a new column to an existing data frame.
#'
#' @param obj A list of Seurat objects.
#' @param count.base Optional data frame to which the counts will be added as a new column.
#' @param col.name Optional character string specifying the name of the new column in count.base.
#' @return A data frame with cell counts per sample, or count.base with an added column of counts.
#' @examples
#' counts <- count_cells_per_sample(obj = seurat_list)
#' counts <- count_cells_per_sample(obj = seurat_list, count.base = meta_df, col.name = "cell_count")
count_cells_per_sample = function(obj = NULL, count.base = NULL, col.name = NULL){

  l = lapply(obj, function(x){
    x@meta.data %>% dplyr::select(orig.ident)
  })
  df = do.call("rbind", l)
  df = df %>%  dplyr::count(orig.ident)
  if(is.null(count.base)) {
    return(df)
  } else {
    count.base[[col.name]] = df$n[match(count.base$orig.ident, df$orig.ident)]
    return(count.base)
  }
}



#' Perform clustering on a Seurat object and return cluster assignments
#'
#' This function normalizes the data, identifies variable features, scales the data,
#' runs PCA, finds neighbors, and clusters the cells using Seurat's workflow.
#' It returns the cluster assignments for each cell.
#'
#' @param sobj A Seurat object containing single-cell data.
#' @return A vector of cluster assignments (seurat_clusters) for each cell.
#' @examples
#' clusters <- get_soup_groups(seurat_object)
#' @export
get_soup_groups <- function(sobj){
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(
    object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst'
  )
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  return(sobj@meta.data[['seurat_clusters']])
}


#' Perform XSoupe.
#'
#' @param dir Character. Path to the directory containing input data.
#' @param id Character or numeric. Identifier for the dataset or sample.
#' @return A soup object constructed from the provided directory and id.
#' @export
make_soup <- function(dir, id){
  message("Correcting for ambient RNA contamination using SoupX for sample: ", id)

  raw_counts_path <- paste0(dir, "sample_raw_feature_bc_matrix")
  filtered_counts_path <- paste0(dir, "sample_filtered_feature_bc_matrix")

  message("Load RAW counts from:", raw_counts_path)
  raw_counts = Read10X(raw_counts_path)
  message("Load filtered counts from:", filtered_counts_path)
  fltrd_counts = Read10X(filtered_counts_path)
  
  message("Check if genes are common between raw and filtered counts")
  common_genes <- intersect(rownames(raw_counts), rownames(fltrd_counts))
  raw_counts <- raw_counts[common_genes, , drop = FALSE]
  fltrd_counts <- fltrd_counts[common_genes, , drop = FALSE]

  message("Create Seurat object from filtered counts")
  seu_obj_filtered = CreateSeuratObject(counts = fltrd_counts, project = id)
  seu_obj_filtered = RenameCells(seu_obj_filtered, new.names = gsub("multi_", "", colnames(seu_obj_filtered)))
  seu_obj_filtered@meta.data$orig.ident = gsub("multi_", "", id)
  message("Calculate soup groups")
  seu_obj_filtered$soup_group <- get_soup_groups(seu_obj_filtered)

  message("Process with soupX")
  sc = SoupChannel(raw_counts, fltrd_counts)
  sc = setClusters(sc, seu_obj_filtered$soup_group)  
  sc <- tryCatch({
    autoEstCont(sc, doPlot = FALSE)
  }, error = function(e) {
    message("autoEstCont failed: ", e$message)
    message("Falling back to manual contamination = 5%")
    setContaminationFraction(sc, 0.05)
  })
  out = adjustCounts(sc, roundToInt = TRUE)

  message("Keep original counts and return Seurat object")
  seu_obj_filtered[["original.counts"]] <- CreateAssayObject(counts = fltrd_counts)
  seu_obj_filtered[["RNA"]]$counts <- out
  seu_obj_filtered

  return(seu_obj_filtered)
} 


#' Perform Doublet Detection with scDblFinder
#'
#' This function applies the scDblFinder algorithm to a Seurat object to identify potential doublets in single-cell RNA-seq data.
#' It extracts count data from the Seurat object, runs scDblFinder, and adds the resulting doublet scores and classifications
#' as metadata to the Seurat object.
#'
#' @param seu.obj A Seurat object containing single-cell RNA-seq data.
#'
#' @return A Seurat object with additional metadata columns:
#'   \describe{
#'     \item{scDblFinder_score}{The doublet score assigned by scDblFinder.}
#'     \item{scDblFinder_class}{The doublet classification ("doublet" or "singlet") assigned by scDblFinder.}
#'   }
#'
#' @details
#' The scDblFinder algorithm performs simulation-based doublet detection, informed by clustering, library sizes,
#' and neighborhoods in dimensionality-reduced space (PCA or UMAP).
#'
#' @seealso \code{\link[scDblFinder]{scDblFinder}}, \code{\link[Seurat]{AddMetaData}}
#'
#' @examples
#' \dontrun{
#'   seu <- perform_scDblFinder(seu)
#' }
perform_scDblFinder <- function(seu.obj){
  sce = scDblFinder(GetAssayData(seu.obj, layer = "counts"), dbr.sd=1)
  # simulation-based doublet detection, informed by:
  # Clustering, Library sizes, Neighborhoods in dimensionality-reduced space (PCA or UMAP)
  df = data.frame(sce@colData) %>% dplyr::select(scDblFinder.score, scDblFinder.class)
  colnames(df) = c("scDblFinder_score", "scDblFinder_class")
  seu.obj = AddMetaData(seu.obj, df)
  return(seu.obj)
}

#' Detect doublets in a Seurat object using SCDS methods
#'
#' This function takes a Seurat object, converts it to a SingleCellExperiment,
#' and applies the SCDS doublet detection methods: cxds, bcds, and cxds_bcds_hybrid.
#' The resulting doublet scores and calls are added to the Seurat object's metadata.
#'
#' @param seu.obj A Seurat object containing single-cell RNA-seq data.
#'
#' @return A Seurat object with additional metadata columns:
#'   \describe{
#'     \item{cxds_score}{Doublet score from the cxds method}
#'     \item{bcds_score}{Doublet score from the bcds method}
#'     \item{hybrid_score}{Doublet score from the cxds_bcds_hybrid method}
#'     \item{cxds_call}{Doublet call (logical) from the cxds method}
#'     \item{bcds_call}{Doublet call (logical) from the bcds method}
#'     \item{hybrid_call}{Doublet call (logical) from the hybrid method}
#'   }
#'
#' @importFrom scds cxds bcds cxds_bcds_hybrid
#' @importFrom SingleCellExperiment as.SingleCellExperiment
#' @importFrom dplyr select
#' @importFrom Seurat AddMetaData
#' @export
scds_doublets <- function(seu.obj) {
    message("Detecting doublets using SCDS methods...")
    sce <- as.SingleCellExperiment(seu.obj)
    message("Using cxds_bcds_hybrid")
    sce <- scds::cxds_bcds_hybrid(sce, estNdbl = TRUE)

    df <- as.data.frame(colData(sce)) %>%
      dplyr::select(hybrid_score,
                    hybrid_call)

    seu.obj <- AddMetaData(seu.obj, df)
    message("Done")
    return(seu.obj)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
integration = function(
    obj,
    obj.l = NULL,
    no.ftrs = 500,
    .assay = "RNA",
    max.cl = 1,
    min.cells.per.sample = 25,
    threads = 5,
    .nbr.dims = 15,
    do.cluster = T,
    do.dimreduc = T,
    cc.regr = F,
    hvg.union = T,
    custom.features = NULL,
    run.integration = T,
    harmony.group.vars = NULL,
    obj.split.by = "orig.ident",
    .perp = 50,
    dmreduc.dims = NULL,
    n.neighbors = 30,
    min.dist = 0.3) {

  library(SignatuR)
  library(parallel)
  library(BiocParallel)
  library(harmony)

  if (any(!"readgmt" %in% installed.packages())) {
    Sys.unsetenv("GITHUB_PAT")
    devtools::install_github("jhrcook/readgmt")
  }
  library(readgmt)


  start.time <- Sys.time()

  if(is.null(dmreduc.dims)) {
    dmreduc.dims = .nbr.dims
  }

  DefaultAssay(obj) = "RNA"
  obj@meta.data = droplevels(obj@meta.data)
  obj = DietSeurat(obj, counts = TRUE, data = TRUE)
  obj = NormalizeData(obj)

  # Gene categories to exclude from variable genes
  bl <- c(
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Mito)[[1]]
    #SignatuR::GetSignature(SignatuR$Hs$Compartments$Immunoglobulins)[[1]]
    #SignatuR::GetSignature(SignatuR$Hs$Compartments$TCR)[[1]]
    # SignatuR::GetSignature(SignatuR$Hs$Blocklists)[[1]]
  )
  # bl = c(bl, c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "XIST")) # gender genes
  # bl <- unique(bl)

  if (hvg.union == T) {

    if(is.null(obj.l)) {
      print("Split object")
      obj.l = Split_Object(obj, split.by = obj.split.by, threads = threads)
    }

    select.bool = unlist(lapply(obj.l, function(x){ncol(x) >= min.cells.per.sample}))
    print(table(select.bool))
    obj.l = obj.l[select.bool]
    length(obj.l)

    print("HVG")
    obj.l = parallel::mclapply(obj.l, function(x) {
      x = x[!rownames(x) %in% bl, ]
      x = FindVariableFeatures(x, selection.method = "vst",  assay = .assay, verbose = FALSE)
      x
    }, mc.cores = threads)

    features = SelectIntegrationFeatures(object.list = obj.l, nfeatures = no.ftrs)

    VariableFeatures(obj) = features
    rm(obj.l); gc()

  } else if (!is.null(custom.features)) {
    VariableFeatures(obj) = custom.features
  } else {
    obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = no.ftrs, assay = .assay)
  }

  if (cc.regr == T) {
    obj = ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), assay = .assay)
    # obj = ScaleData(obj, vars.to.regress = c("Perc_of_mito_genes"), assay = .assay)
  } else {
    obj = ScaleData(obj, assay = .assay)
  }

  obj = RunPCA(obj, assay = .assay)
  p1 <- plot(ElbowPlot(obj, ndims = 50))

  # Assuming you've run PCA already on 'obj'
  pca_stdev <- obj[["pca"]]@stdev

  # Calculate variance explained (%)
  var_explained <- (pca_stdev^2 / sum(pca_stdev^2)) * 100

  # Cumulative variance explained (optional)
  cum_var_explained <- cumsum(var_explained)

  # Create data frame for plotting
  pc_df <- data.frame(
    PC = 1:length(var_explained),
    VarianceExplained = var_explained,
    CumulativeVariance = cum_var_explained
  )

  # Plot with ggplot2
  library(ggplot2)
  p1 <- ggplot(pc_df, aes(x = PC, y = VarianceExplained)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_line(aes(y = CumulativeVariance), color = "red", size = 1) +
    geom_point(aes(y = CumulativeVariance), color = "red") +
    labs(
      title = "Variance Explained by Principal Components",
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    theme_minimal()
  ggsave("PCA_variance_explained.png", plot = p1, width = 8, height = 6)

  print(paste0("### PCs used for harmony clustering and dim reduc: ", .nbr.dims, " ###"))

  if(run.integration == T){
    obj = RunHarmony(
      obj, group.by.vars = harmony.group.vars, # theta = c(2,3),
      reduction.use ='pca', dims.use = 1:.nbr.dims, max_iter = 15, ncores = threads
    )
    comp.wrk = 'harmony'
  } else {
    comp.wrk = 'pca'
  }

  if (do.cluster == T) {
    print("Find clusters")
    obj = FindNeighbors(obj, reduction = comp.wrk, dims = 1:.nbr.dims)

    reso = seq(0,max.cl,.1)
    names(reso) = reso
    suppressWarnings({
      suppressMessages({
        findclusters.res = parallel::mclapply(reso, function(x) {
          FindClusters(obj, resolution = x)@meta.data[, "seurat_clusters", drop = F]
        }, mc.cores = length(reso))
      })
    })
    res.names = names(findclusters.res)
    findclusters.res = do.call("cbind", findclusters.res)
    colnames(findclusters.res) = paste0("RNA_snn_res.", res.names)
    stopifnot(identical(rownames(obj@meta.data), rownames(findclusters.res)))
    obj = AddMetaData(obj, findclusters.res)
  }
  if (do.dimreduc == T) {
    # print("tSNE")
    # obj = RunTSNE(
    #   obj, reduction = comp.wrk, dims = 1:dmreduc.dims, seed.use = 1234,
    #   nthreads = threads, tsne.method = "FIt-SNE"
    # )
    print("UMAP")
    obj = RunUMAP(
      obj, reduction = comp.wrk, dims = 1:dmreduc.dims, seed.use = 1234,
      min.dist = min.dist, n.neighbors = n.neighbors
    )
  }

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)

  return(obj)
}

label_cells_rm = function(obj, nFeature_low_cutoff, nFeature_high_cutoff, nCount_low_cutoff, nCount_high_cutoff, mt_cutoff, complx_cutoff) {
  obj@meta.data = obj@meta.data %>% mutate(
    KEEP_CELL = case_when(
      (nFeature_RNA < nFeature_low_cutoff) | (nFeature_RNA > nFeature_high_cutoff) |
        (nCount_RNA < nCount_low_cutoff) | (nCount_RNA > nCount_high_cutoff) |
        (Perc_of_mito_genes > mt_cutoff)  | (log10GenesPerUMI < complx_cutoff) ~ FALSE,
      TRUE ~ TRUE
    )
  )
  obj
}
