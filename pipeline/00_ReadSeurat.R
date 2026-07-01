#!/usr/bin/env Rscript
# -------------------------------------------------------------------

# Load required packages

# All packages should be installed via environment.yml

# -------------------------------------------------------------------

.required_packages <- c(
"Seurat",
"yaml",
"ggplot2",
"patchwork",
"dplyr",
"here",
"doParallel",
"parallel",
"data.table",
"Matrix",
"SoupX",
"stringr",
"naturalsort",
"scales",
"openxlsx",
"cowplot",
"ggthemes",
"SingleCellExperiment",
"scDblFinder",
"dittoSeq",
"clustifyr",
"scds",
"UCell"
)

# Check that all required packages are installed
.missing_packages <- .required_packages[
!sapply(.required_packages, requireNamespace, quietly = TRUE)
]

if (length(.missing_packages) > 0) {
stop(
"Missing required packages:\n",
paste(.missing_packages, collapse = ", "),
"\n\nPlease update your Conda/Mamba environment:\n",
"mamba env update -f environment.yml",
call. = FALSE
)
}

# Load libraries

invisible(lapply(
.required_packages,
library,
character.only = TRUE
))

# Print session information for reproducibility
message("R version: ", R.version.string)
message("Loaded ", length(.required_packages), " packages.")
message("\nLoaded package versions:")
for (pkg in .required_packages) {
  message(sprintf(" %-30s %s",
  pkg,
  as.character(packageVersion(pkg))))
}
message("==============================================\n")

# -------------------------------------------------------------------

# Set up parallel processing on cluster

# -------------------------------------------------------------------

# Set up parallel processing on cluster 
if (Sys.info()["nodename"] == "ribnode020") {
  ncores = 35
} else {
  ncores = 5
}

# Read in configs and helper functions
manifest_path <- here("config/configs.yml")
manifest <- yaml.load_file(manifest_path)
source(here("scripts/preprocess.R"))
source(here("scripts/helper.R"))
source(here("scripts/plots.R"))
batch_name <- "Batch1"
# Define output directories
seurat.path <- manifest$input$cellranger
seurat.output <- manifest$output$processed_seurat
qc_dir <- make_outdir("results/qc")
seurat_dir <- make_outdir(seurat.output)

# -------------------------------------------------------------------

# Check Cell Ranger output directories and process metrics

# -------------------------------------------------------------------
cellranger_paths <- list(
  Batch1 = seurat.path 
  # Add more batches like: Batch2 = manifest$cellranger$batch2
)
# metrics_all <- lapply(names(cellranger_paths), function(batch) {
#   process_cellranger_metrics(cellranger_paths[[batch]], batch_name = batch)
# })
# combined_metrics <- bind_rows(metrics_all)

# p_metrics <- generate_cellranger_metrics_plot(
#   combined_metrics = combined_metrics,
#   output_path = qc_dir,
#   file_name = "metrics_comparison_across_samples_batch1.png"
# )

######## Read cellranger outputs and create Seurat object
for (batch_name in names(cellranger_paths)) {
  seurat.path <- cellranger_paths[[batch_name]]
  
  message("Processing batch: ", batch_name)

  cellranger.dirs <- list.dirs(path = seurat.path, full.names = TRUE, recursive = FALSE)
  cellranger.samples <- basename(cellranger.dirs)
  dirs <- file.path(cellranger.dirs, "count/sample_filtered_feature_bc_matrix/")
  names(dirs) <- cellranger.samples

  seurat.l <- parallel::mclapply(names(dirs), function(sample) {
  # seurat.l <- lapply(names(dirs), function(sample) {
    dir <- dirs[[sample]]
    id <- sample
    message("  Sample: ", id)
    fltrd_counts <- Read10X(dir)
    seu_obj_filtered <- CreateSeuratObject(counts = fltrd_counts, project = batch_name)
    seu_obj_filtered <- RenameCells(seu_obj_filtered, new.names = gsub("multi_", "", colnames(seu_obj_filtered)))
    seu_obj_filtered@meta.data$orig.ident <- gsub("multi_", "", id)

    return(seu_obj_filtered)
  }, mc.cores = length(dirs))

  names(seurat.l) <- names(dirs)

  # Merge all samples in the batch into one Seurat object
  se.meta <- if (length(seurat.l) == 1) {
    seurat.l[[1]]
  } else {
    merge(
      seurat.l[[1]], y = seurat.l[2:length(seurat.l)],
      add.cell.ids = names(seurat.l),
      project = batch_name
    )
  }

  se.meta@meta.data$orig.ident <- factor(se.meta@meta.data$orig.ident)
  se.meta[["RNA"]] <- JoinLayers(se.meta[["RNA"]])

  # Save the merged Seurat object
  output.file <- file.path(seurat_dir, paste0("00_seurat_merged_", (batch_name), ".Rds"))
  saveRDS(se.meta, file = output.file)

  message("Saved Seurat object for ", batch_name, " to: ", output.file)
}
