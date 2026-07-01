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
"UCell",
"scGate"
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
integration_dir <- make_outdir("results/integration")
seurat_dir <- make_outdir(seurat.output)

se.meta <- readRDS(file.path(seurat_dir, paste0("02_annotated_", batch_name, ".Rds")))

##### Integrate samples 
se.meta = CellCycleScoring(
  se.meta, s.features = s.genes, g2m.features = g2m.genes,
  assay = 'RNA', search = TRUE
)
Idents(se.meta) <- "Final_Annotation"
nk_cells <- subset(se.meta, idents = c("NK"))
nk_cells$orig.ident <- droplevels(factor(nk_cells$orig.ident))
table(nk_cells$orig.ident, nk_cells$Final_Annotation)

dims.use.rna = 15
nk_cells = integration(
  obj = nk_cells,
  no.ftrs = 2000,
  threads = 20,
  .nbr.dims = dims.use.rna,
  cc.regr = T, 
  run.integration = T, # harmony instead of pca (?)
  harmony.group.vars = c("orig.ident"),
)

col1 <- DimPlot(
  nk_cells,
  group.by = "orig.ident",
  reduction = "umap",
  pt.size = 0.05
) + coord_fixed()
col2 <- DimPlot(
  nk_cells,
  group.by = "RNA_snn_res.0.3",
  reduction = "umap",
  pt.size = 0.05
) + coord_fixed()

combined_plot <- col1 + col2 + plot_layout(ncol = 2)
combined_plot

ggsave(
    filename = file.path(integration_dir, "NK_cells_UMAP.png"),
    plot = combined_plot,
    width = 180,
    height = 100,
    bg = "white",
    dpi = 300,
    units = "mm"
)

seu.obj.path <- paste0(manifest$output$processed_seurat, "03_integrated_nk_", batch_name, ".Rds")

saveRDS(nk_cells, file = seu.obj.path)

message("Finished batch: ", batch_name)
