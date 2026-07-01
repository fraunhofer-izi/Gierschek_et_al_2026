#' Run Differential Expression Analysis Between Two Groups
#'
#' This function performs differential gene expression analysis
#' from a Seurat object, comparing two groups (e.g., CAR vs NT) using the MAST test. It
#' adjusts for technical covariates, saves results of signficiant genes, and returns the top 10 up- and downregulated genes.
#'
#' @param seu.ob A Seurat object containing cells to subset and analyze.
#' @param group_map A named list defining group assignments. Each name (e.g., "CAR", "NT") 
#'   corresponds to a vector of `orig.ident` sample names.
#' @param group_labels A character vector of length 2 defining which groups to compare (default: c("CAR", "NT")).
#' @param output_filename The name of the CSV file to save the full DEG result.
#' @param min_pct Minimum fraction of cells expressing the gene in either group to be tested (default: 0.2).
#' @param logfc_threshold Log2 fold change threshold to call a gene up-/downregulated (default: 1.0).
#' @param padj_threshold Adjusted p-value threshold to call significance (default: 0.05).
#' @param result_dir Directory path where DEG results will be saved (default: "./results/DEG/").
#'
#' @return A named list with:
#'   \item{deg_results}{A data.frame of all genes tested with log2FC, adjusted p-values, and significance labels.}
#'   \item{top_genes}{A data.frame of the top 10 upregulated and top 10 downregulated genes.}
#'
#' @examples
#' group_map <- list(
#'   CAR = c("SmS_CD123-OCIAML2", "SmS_CD33-OCIAML2", "Prod_CD123-OCIAML2"),
#'   NT = "NTNK-OCIAML2"
#' )
#' results <- run_deg_analysis(
#'   seu.ob = seu.ob,
#'   group_map = group_map,
#'   output_filename = "deg_mast_CAR_vs_NT.csv"
#' )
#'
#' @export

run_deg_analysis <- function(
  seu.obj,
  group_map,
  group_labels = c("CAR", "NT"),
  output_filename,
  min_pct = 0.2,
  logfc_threshold = log2(1.25),
  padj_threshold = 0.05,
  result_dir = "../results/DEG/",
  downsample_to_smallest = FALSE,
  downsample_n = NULL  # new argument
) {
  # Subset NK cells based on input group_map keys
  subset_ids <- unlist(group_map)
  seu.ob_subset <- subset(seu.obj, subset = orig.ident %in% subset_ids)
  # seu.ob_subset$orig.ident <- droplevels(seu.ob_subset$orig.ident)
  
  # Map orig.ident to comparison group
  seu.ob_subset$comparison_group <- as.character(seu.ob_subset$orig.ident)
  for (grp in names(group_map)) {
    seu.ob_subset$comparison_group[seu.ob_subset$comparison_group %in% group_map[[grp]]] <- grp
  }

  Idents(seu.ob_subset) <- "comparison_group"
  # 🧪 Optional downsampling
  if (downsample_to_smallest || !is.null(downsample_n)) {
    group_sizes <- table(seu.ob_subset$comparison_group)
    
    if (!is.null(downsample_n)) {
      # Downsample to user-specified number
      target_n <- downsample_n
      message("Downsampling each group to fixed size: ", target_n, " cells.")
    } else {
      # Downsample to smallest group
      target_n <- min(group_sizes)
      message("Downsampling each group to smallest group size: ", target_n, " cells.")
    }
    
    # Warn if any group has fewer cells than target_n
    too_small <- names(group_sizes[group_sizes < target_n])
    if (length(too_small) > 0) {
      warning("The following groups have fewer cells than the target (", target_n, "): ",
              paste(too_small, collapse = ", "), ". They will be kept as-is (no upsampling).")
    }

    balanced_cells <- unlist(lapply(group_labels, function(lbl) {
      cells <- WhichCells(seu.ob_subset, idents = lbl)
      if (length(cells) > target_n) {
        sample(cells, size = target_n)
      } else {
        cells
      }
    }))
    seu.ob_subset <- subset(seu.ob_subset, cells = balanced_cells)
  }
  # Run differential expression analysis
  deg_mast <- FindMarkers(
    seu.ob_subset,
    ident.1 = group_labels[1],
    ident.2 = group_labels[2],
    test.use = "MAST",
    latent.vars = c("nCount_RNA"),
    min.pct = min_pct
  )

  # Add gene column
  deg_mast$gene <- rownames(deg_mast)

  # Annotate significance
  deg_mast$significance <- "Not Significant"
  deg_mast$significance[deg_mast$avg_log2FC > logfc_threshold & deg_mast$p_val_adj < padj_threshold] <- "Upregulated"
  deg_mast$significance[deg_mast$avg_log2FC < -logfc_threshold & deg_mast$p_val_adj < padj_threshold] <- "Downregulated"

  # Save CSV
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }
  write.csv(deg_mast, file = file.path(result_dir, output_filename), row.names = TRUE)

  # Get top genes
  top_up <- deg_mast[deg_mast$significance == "Upregulated", ]
  top_up <- top_up[order(-top_up$avg_log2FC), ][1:min(10, nrow(top_up)), ]

  top_down <- deg_mast[deg_mast$significance == "Downregulated", ]
  top_down <- top_down[order(top_down$avg_log2FC), ][1:min(10, nrow(top_down)), ]

  top_genes <- rbind(top_up, top_down)

  return(list(
    deg_results = deg_mast,
    top_genes = top_genes
  ))
}


#' Generate a Volcano Plot from DEG Results with Up/Downregulated Counts
#'
#' Creates a volcano plot from differential expression results. Top genes can be labeled,
#' and the title will include counts of up- and downregulated genes automatically.
#'
#' @param deg_results A data frame with DEG results, including `avg_log2FC`, `p_val_adj`, `gene`, and `significance`.
#' @param top_genes A data frame of genes to label in the plot (e.g., top up- and downregulated).
#' @param title Base title of the plot (counts are appended automatically).
#' @param logfc_threshold Numeric. Threshold for vertical lines (default: 1.0).
#' @param pval_threshold Numeric. Adjusted p-value threshold for horizontal line (default: 0.05).
#' @param xlim Optional numeric vector for x-axis limits (default: c(-4, 4)).
#' @param ylim Optional numeric vector for y-axis limits (default: NULL).
#'
#' @return A ggplot2 object representing the volcano plot.
#' @export
plot_volcano <- function(deg_results, top_genes,
                         title = "Volcano Plot",
                         logfc_threshold = 1.0,
                         pval_threshold = 0.05,
                         xlim = c(-4, 4),
                         ylim = NULL,
                         col_up = "#8d8dff",
                         col_down = "#56e473",
                         col_ns = "grey",
                         base.size = 14) {

  # Count up/downregulated genes
  upregulated_genes <- sum(deg_results$significance == "Upregulated", na.rm = TRUE)
  downregulated_genes <- sum(deg_results$significance == "Downregulated", na.rm = TRUE)

  # Append counts to title
  full_title <- paste0(title, " (up=", upregulated_genes, ", down=", downregulated_genes, ")")

  # Create plot
  p <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", color = "black") +
    scale_color_manual(values = c("Upregulated" = col_up,
                                  "Downregulated" = col_down,
                                  "Not Significant" = col_ns)) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      max.overlaps = 10,
      size = base.size / 3,
      show.legend = FALSE
    ) +
    labs(
      title = full_title,
      x = "Log2 Fold Change",
      y = "-log10(Adjusted P-value)",
      color = "Significance"
    ) +
    theme_bw(base_size = base.size) +   # match GSEA style
    theme(
      # legend.position = c(0.05, 0.95),
      legend.position = "top",
      legend.direction = "horizontal",
      # legend.justification = c("left", "top"),
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
      legend.title = element_text(size = base.size * 0.7, face = "bold"),
      legend.text = element_text(size = base.size * 0.65),
      plot.title = element_text(hjust = 0.5, size = base.size * 1.1),
      axis.title.x = element_text(size = base.size),     # larger x-axis title
      axis.title.y = element_text(size = base.size),     # larger y-axis title
      axis.text.x  = element_text(size = base.size * 0.8),
      axis.text.y  = element_text(size = base.size * 0.8)
    )

  # Apply axis limits
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }

  return(p)
}



mytheme = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(colour = "black"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = .3),
      axis.line.x.top = element_blank(),
      axis.line.y.right = element_blank(),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(
        hjust = 0, face = "plain", colour = "black", size = rel(1)
      ),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}