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
set.seed(2026)  # For reproducibility
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

######## Calculate QC metrics and filter cells
se.meta <- readRDS(file.path(seurat_dir, paste0("00_seurat_merged_", (batch_name), ".Rds")))
# Add QC metrics
se.meta[["Perc_of_mito_genes"]] <- PercentageFeatureSet(se.meta, pattern = "^MT-")
se.meta[["Perc_of_ribosomal_genes"]] <- PercentageFeatureSet(se.meta, pattern = "^RPL|^RPS")
se.meta$log10GenesPerUMI <- log10(se.meta$nFeature_RNA) / log10(se.meta$nCount_RNA)

# Generate QC violin plots
generate_qc_plots(se.meta, batch_name = batch_name, manifest = manifest, output_path = qc_dir)

# Track cell counts before filtering
cell.track <- count_cells_per_sample(c(se.meta))

qc <- manifest$qc_cutoffs[[batch_name]]
nFeature_low_cutoff  <- qc$nFeature_low
nFeature_high_cutoff <- qc$nFeature_high
nCount_low_cutoff    <- qc$nCount_low
nCount_high_cutoff   <- qc$nCount_high
mt_cutoff            <- qc$mt_cutoff
complx_cutoff        <- qc$complx_cutoff
# Label cells based on QC cutoffs
se.meta <- label_cells_rm(se.meta, nFeature_low_cutoff, nFeature_high_cutoff, nCount_low_cutoff, nCount_high_cutoff, mt_cutoff, complx_cutoff)

# Count after filtering
cell.track <- count_cells_per_sample(c(subset(se.meta, subset = KEEP_CELL == TRUE)), cell.track, "low_quality")

# Apply QC filtering
se.meta <- subset(se.meta, subset = KEEP_CELL == TRUE)
se.meta@meta.data <- droplevels(se.meta@meta.data)
se.meta$KEEP_CELL <- NULL

# Split for doublet calling
seurat.l <- Split_Object(se.meta, split.by = "orig.ident")

# Run scDblFinder and scds in parallel
seurat.l <- parallel::mclapply(seurat.l, perform_scDblFinder, mc.cores = length(seurat.l))
seurat.l <- parallel::mclapply(seurat.l, scds_doublets, mc.cores = length(seurat.l))

# Add consensus doublet labels
for (i in seq_along(seurat.l)) {
  seurat.l[[i]]$DOUBLETS_CONSENSUS_Strict  <- seurat.l[[i]]$scDblFinder_class == "doublet" & seurat.l[[i]]$hybrid_call == TRUE
  seurat.l[[i]]$DOUBLETS_CONSENSUS_Conserv <- seurat.l[[i]]$scDblFinder_class == "doublet" | seurat.l[[i]]$hybrid_call == TRUE
}

# Merge back to single object
se.meta <- if (length(seurat.l) == 1) {
  seurat.l[[1]]
} else {
  merge(seurat.l[[1]], y = seurat.l[2:length(seurat.l)], project = batch_name)
}
se.meta@meta.data$orig.ident <- factor(se.meta@meta.data$orig.ident)
se.meta[["RNA"]] <- JoinLayers(se.meta[["RNA"]])

# Subset doublets
se.meta_conservative <- subset(se.meta, subset = DOUBLETS_CONSENSUS_Conserv == FALSE)
se.meta_conservative@meta.data <- droplevels(se.meta_conservative@meta.data)
se.meta_conservative$DOUBLETS_CONSENSUS_Conserv <- NULL
cell.track <- count_cells_per_sample(c(se.meta_conservative), cell.track, "doublets_conserv")

se.meta_strict <- subset(se.meta, subset = DOUBLETS_CONSENSUS_Strict == FALSE)
se.meta_strict@meta.data <- droplevels(se.meta_strict@meta.data)
se.meta_strict$DOUBLETS_CONSENSUS_Strict <- NULL
cell.track <- count_cells_per_sample(c(se.meta_strict), cell.track, "doublets_strict")

# Normalize filtered object (conservative)
se.meta <- NormalizeData(se.meta_conservative, assay = "RNA", normalization.method = "LogNormalize")
se.meta$DOUBLETS_CONSENSUS_Strict <- NULL

# Save final cell tracking table
celltrack_file <- file.path(qc_dir, paste0("stats_celltrack_", (batch_name), ".Rds"))
saveRDS(cell.track, file = celltrack_file)

# Save final filtered Seurat object
output.file <- file.path(manifest$output$processed_seurat, paste0("01_seurat_filtered_", (batch_name), ".Rds"))
saveRDS(se.meta, file = output.file)

message("Finished batch: ", batch_name)
