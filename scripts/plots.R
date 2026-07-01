generate_genes_umis_plot <- function(se.meta,
                                     batch_name = NULL,
                                     output_path = "results/qc") {

  # Extract metadata
  qc_box_data <- se.meta@meta.data %>%
    dplyr::select(orig.ident, nCount_RNA, nFeature_RNA)

  # UMIs per cell
  p_box_nCount <- ggplot(
    qc_box_data,
    aes(x = orig.ident, y = nCount_RNA)
  ) +
    geom_boxplot(outlier.size = 0.5, fill = "steelblue") +
    labs(
      title = "UMIs per Cell (After QC)",
      x = "Sample",
      y = "nCount_RNA (UMIs)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  # Genes per cell
  p_box_nFeature <- ggplot(
    qc_box_data,
    aes(x = orig.ident, y = nFeature_RNA)
  ) +
    geom_boxplot(outlier.size = 0.5, fill = "tomato") +
    labs(
      title = "Genes per Cell (After QC)",
      x = "Sample",
      y = "nFeature_RNA (Genes)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  # Combine plots
  plot <- cowplot::plot_grid(
    p_box_nCount,
    p_box_nFeature,
    labels = "AUTO",
    ncol = 2
  )

  # Create filename
  filename <- if (!is.null(batch_name)) {
    paste0("genes_and_umis_", batch_name, ".png")
  } else {
    "genes_and_umis.png"
  }

  output_file <- file.path(output_path, filename)

  message("Saving plot to: ", output_file)

  ggsave2(
    filename = here(output_file),
    plot = plot,
    width = 300,
    height = 150,
    dpi = 300,
    bg = "white",
    units = "mm"
  )

  return(plot)
}

generate_cellranger_metrics_plot <- function(
  combined_metrics,
  output_path = "results/qc",
  file_name = "metrics_comparison_across_samples.png"
) {

  p_metrics <- ggplot(
  combined_metrics,
  aes(
  x = Sample_Clean,
  y = Mean_Value,
  fill = Batch
  )
  ) +
  geom_bar(
  stat = "identity",
  position = position_dodge()
  ) +
  facet_wrap(
  ~ Metric.Name,
  scales = "free_y"
  ) +
  labs(
  title = "Per-Sample Metrics",
  x = "Sample",
  y = "Value",
  fill = "Batch"
  ) +
  theme_minimal() +
  theme(
  axis.text.x = element_text(
  angle = 45,
  hjust = 1
  )
  )

  output_file <- file.path(
  output_path,
  file_name
  )

  message("Saving plot to: ", output_file)

  ggsave2(
  filename = here(output_file),
  plot = p_metrics,
  width = 250,
  height = 180,
  dpi = 300,
  units = "mm",
  bg = "white"
  )

return(p_metrics)
}


generate_qc_plots <- function(se.meta, batch_name, manifest, output_path="results/qc") {
  cutoff <- manifest$qc_cutoffs[[batch_name]]

  p1 <- qc_vln_plot_cell(se.meta, low_cutoff = cutoff$nFeature_low, high_cutoff = cutoff$nFeature_high)
  p2 <- qc_vln_plot_cell(se.meta, .features = "nCount_RNA", plot_title = "UMIs per Cell", y_axis_label = "UMIs",
                         low_cutoff = cutoff$nCount_low, high_cutoff = cutoff$nCount_high)
  p3 <- qc_vln_plot_cell(se.meta, .features = "Perc_of_mito_genes", plot_title = "Mito Gene % per Cell",
                         y_axis_label = "% Mito Gene Counts", high_cutoff = cutoff$mt_cutoff)
  p4 <- qc_vln_plot_cell(se.meta, .features = "log10GenesPerUMI", plot_title = "Cell Complexity",
                         y_axis_label = "log10(Genes) / log10(UMIs)", high_cutoff = cutoff$complx_cutoff)

  plot <- cowplot::plot_grid(p1, p2, p3, p4, scale = 0.9, nrow = 1,
                             labels = "AUTO", label_fontface = "bold", label_size = 14)

  output_path <- file.path(output_path, paste0("stats_tech_per_cell_", batch_name, ".png"))
  message("Saving plot to: ", output_path)
  ggsave2(
    filename = here(output_path),
    plot,
    width = 180,
    height = 80,
    dpi = 300,
    bg = "white",
    units = "mm",
    scale = 2
  )

  return(plot)
}



#' Plot GSEA NES barplots grouped by biological themes
#'
#' @export
plot_gsea_bars_with_theme <- function(
    df,
    top_n = 10,
    fdr = 0.05,
    wrap = 35,
    col_up = "#8D8DFF",
    col_down = "#56e473",
    base.size = 14,
    grouping = c("facet", "bracket"),
    label_up = "Up (CD123-CAR)",
    label_down = "Up (CD33-CAR)"
) {

  grouping <- match.arg(grouping)

  # ------------------------------------------------------------
  # Filter + prepare
  # ------------------------------------------------------------
  df1 <- df %>%
    dplyr::filter(!is.na(NES)) %>%
    dplyr::filter(is.na(padj) | padj <= fdr) %>%
    dplyr::mutate(

      dir = dplyr::if_else(
        NES >= 0,
        label_up,
        label_down
      ),

      term_clean = term %>%
        gsub("^GOBP_", "", .) %>%
        gsub("_", " ", .),

      term_wrapped = stringr::str_wrap(
        term_clean,
        width = wrap
      ),

      # IMPORTANT:
      # preserve manual \n linebreaks already present
      Theme_wrapped = if ("Theme" %in% colnames(.))
        Theme else NA_character_
    )

  # ------------------------------------------------------------
  # Select top pathways
  # ------------------------------------------------------------
  up <- df1 %>%
    dplyr::filter(NES > 0) %>%
    dplyr::slice_max(
      NES,
      n = top_n
    )

  down <- df1 %>%
    dplyr::filter(NES < 0) %>%
    dplyr::slice_min(
      NES,
      n = top_n
    )

  sel <- dplyr::bind_rows(
    down,
    up
  ) %>%
    dplyr::arrange(NES)

  if (nrow(sel) == 0) {
    stop("No terms pass filters.")
  }

  # ------------------------------------------------------------
  # Dynamic axis scaling
  # ------------------------------------------------------------
  max_nes <- max(sel$NES, na.rm = TRUE)
  min_nes <- min(sel$NES, na.rm = TRUE)

  lim <- max(abs(sel$NES), na.rm = TRUE) * 1.05

  if (min_nes >= 0) {

    x_limits <- c(
      0,
      max_nes * 1.05
    )

  } else if (max_nes <= 0) {

    x_limits <- c(
      min_nes * 1.05,
      0
    )

  } else {

    x_limits <- c(
      -lim,
      lim
    )
  }

  # ------------------------------------------------------------
  # Main plot
  # ------------------------------------------------------------
  p <- ggplot(
    sel,
    aes(
      x = NES,
      y = forcats::fct_reorder(
        term_wrapped,
        NES
      ),
      fill = dir
    )
  ) +

    geom_col() +

    geom_vline(
      xintercept = 0,
      color = "grey40"
    ) +

    scale_fill_manual(
      values = setNames(c(col_down, col_up), c(label_down, label_up)),
      breaks = c(label_down, label_up)
    )+

    scale_x_continuous(
      limits = x_limits
    ) +

    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      fill = NULL
    ) +

    theme_bw(
      base_size = base.size
    ) +

    theme(

      legend.position = "top",

      legend.margin = margin(
        t = -5,
        b = -5
      ),

      legend.box.margin = margin(
        b = -5
      ),

      legend.justification = "center",

      panel.grid = element_blank(),

      panel.border = element_blank(),

      axis.line = element_line(),

      axis.text.y = element_text(
        size = base.size * 0.8
      )
    )

  # ------------------------------------------------------------
  # Grouping
  # ------------------------------------------------------------
  if ("Theme_wrapped" %in% colnames(sel)) {

    # --------------------------------------------------------
    # Facets
    # --------------------------------------------------------
    if (grouping == "facet") {

      p <- p +

        facet_grid(
          Theme_wrapped ~ .,
          scales = "free_y",
          space = "free_y",
          switch = "y",
          labeller = label_value
        ) +

        theme(

          strip.placement = "inside",

          strip.text.y.left = element_text(
            angle = 0,
            hjust = 0,
            size = base.size * 0.75,
            margin = margin(
              r = 3,
              l = 3
            )
          ),

          strip.background = element_rect(
            fill = "grey92",
            colour = NA
          ),

          strip.switch.pad.grid = unit(
            0.05,
            "cm"
          ),

          panel.spacing.y = unit(
            0.05,
            "lines"
          )
        )

    # --------------------------------------------------------
    # Brackets
    # --------------------------------------------------------
    } else if (grouping == "bracket") {

      sel <- sel %>%
        dplyr::mutate(
          ypos = as.numeric(
            forcats::fct_reorder(
              term_wrapped,
              NES
            )
          )
        )

      theme_groups <- sel %>%
        dplyr::group_by(
          Theme_wrapped
        ) %>%
        dplyr::summarise(
          ymin = min(ypos),
          ymax = max(ypos),
          ymid = mean(ypos),
          .groups = "drop"
        )

      p <- p +

        coord_flip() +

        ggforce::geom_mark_bracket(
          data = theme_groups,

          aes(
            xmin = ymin,
            xmax = ymax,
            y.position = max(x_limits) * 1.05,
            label = Theme_wrapped
          ),

          inherit.aes = FALSE,

          label.hjust = 0,
          label.vjust = 0,
          label.angle = 0
        )
    }
  }

  return(p)
}
