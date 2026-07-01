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

# -------------------------------------------------------------------
# Set up parallel processing on cluster
# -------------------------------------------------------------------

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
annotation_dir <- make_outdir("results/annotation")
seurat_dir <- make_outdir(seurat.output)

se.meta <- readRDS(file.path(seurat_dir, paste0("01_seurat_filtered_", batch_name, ".Rds")))

scGate_models_DB = get_scGateDB("data/metadata/scGateDB")
Idents(se.meta) = "orig.ident"

se.meta.l = SplitObject(se.meta, split.by = "orig.ident")
se.meta.l <- mclapply(names(se.meta.l), function(name) {
  obj <- se.meta.l[[name]]
  # Standard Seurat preprocessing
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  return(obj)
}, mc.cores = 8)

names(se.meta.l) <- sapply(se.meta.l, function(obj) unique(obj$orig.ident))

scgate_nk <- scGate_models_DB$human$PBMC[c("NK")]
scgate_myloid <- scGate_models_DB$human$generic[c("Myeloid")]

se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = scgate_nk)
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = scgate_myloid)
# Extract metadata
meta_df <- se.meta[[]]

meta_df$scGate_Annotation <- case_when(
  meta_df$SCGATE__NK == "Pure" ~ "NK",
  meta_df$SCGATE__Myeloid == "Pure" ~ "AML",
  TRUE ~ "Conflict"
)
se.meta <- AddMetaData(se.meta, metadata = meta_df$scGate_Annotation, col.name = "scGate_Annotation")

# Refine scGate classification, NK cell with high AML-like signature, and AML cells with high NK cell like signature
nk_genes <- c("NKG7","PRF1","GNLY","KLRD1","GZMB","GZMH")
se.meta <- AddModuleScore(se.meta, list(nk_genes), name = "NKscore")
aml_genes <- c("CD33","IL3RA","CSF3R","FCGR1A","FLT3","LYZ")
se.meta <- AddModuleScore(se.meta, list(aml_genes), name = "AMLscore")

thr_nk <- quantile(
  se.meta$AMLscore1[se.meta$scGate_Annotation == "NK"], 
  0.99, na.rm = TRUE
)
thr_aml <- quantile(
  se.meta$NKscore1[se.meta$scGate_Annotation == "AML"], 
  0.99, na.rm = TRUE
)

meta_df <- se.meta[[]]
meta_df$scGate_Annotation <- as.character(meta_df$scGate_Annotation)

meta_df$Final_Annotation <- case_when(
  meta_df$scGate_Annotation == "NK"  & meta_df$AMLscore1 > thr_aml  ~ "NK_conflict",
  meta_df$scGate_Annotation == "AML" & meta_df$NKscore1  > thr_nk  ~ "AML_conflict",
  TRUE ~ meta_df$scGate_Annotation
)
se.meta <- AddMetaData(se.meta, metadata = meta_df$Final_Annotation, col.name = "Final_Annotation")

p <- FeatureScatter(
  se.meta,
  feature1 = "NKscore1",
  feature2 = "AMLscore1",
  group.by = "Final_Annotation"
)

# add cutoffs
p <- p + 
  geom_hline(yintercept = thr_aml, linetype = "dashed", color = "red") +
  geom_vline(xintercept = thr_nk, linetype = "dashed", color = "blue") +
  annotate("text", x = max(se.meta$NKscore1), y = thr_aml, label = paste0("AML cutoff: ", round(thr_aml, 2)), hjust = 1, vjust = -0.5, color="red") +
  annotate("text", x = thr_nk, y = max(se.meta$AMLscore1), label = paste0("NK cutoff: ", round(thr_nk, 2)), vjust = 1, hjust = -0.1, color="blue")
p
ggsave(
    filename = file.path(annotation_dir, "NKscore_vs_AMLscore_scatter.png"),
    plot = p,
    width = 300,
    height = 150,
    bg = "white",
    dpi = 300,
    units = "mm"
)

rm(se.meta.l)
gc()

# Save final filtered Seurat object
output.file <- file.path(manifest$output$processed_seurat, paste0("02_annotated_", (batch_name), ".Rds"))
saveRDS(se.meta, file = output.file)

se.meta.l = SplitObject(se.meta, split.by = "orig.ident")
# Define the object names you want to keep
keep_samples <- c("Prod_CD123-OCIAML2", "SmS_CD123-OCIAML2", "SmS_CD33-OCIAML2")

# Filter se.meta.l to only include those objects
se.meta.l <- se.meta.l[names(se.meta.l) %in% keep_samples]

# Optional: subset each of them to Final_Annotation == "NK"
se.meta.l <- lapply(se.meta.l, function(obj) {
  subset(obj, subset = Final_Annotation == "NK")
})

car_genes <- list(
  "Prod_CD123-OCIAML2" = c("CAR-CD123"), 
  "SmS_CD123-OCIAML2" = c("CAR-CD123"), 
  "SmS_CD33-OCIAML2" = c("CAR-CD33")
)

car_expression_by_flow <- list(
  "Prod_CD123-OCIAML2" = c("40.0"),
  "SmS_CD123-OCIAML2" = c("60.0"), 
  "SmS_CD33-OCIAML2" = c("20.0")
  # "SondentestCLEC12A_CD19" = c("50.0",
  #  "50.0")
)

results <- list()

# Loop through each sample and its corresponding genes
for (sample in names(car_genes)) {
  # Get Seurat object
  obj <- se.meta.l[[sample]]
  
  # Get list of genes to evaluate
  genes <- car_genes[[sample]]
  
  # Total number of cells in the sample
  total_cells <- ncol(obj)
  
  # Loop through each gene separately
  for (gene in genes) {
    if (gene %in% rownames(obj)) {
      # Get expression vector (raw data)
      expr_values <- Seurat::GetAssayData(obj, layer = "counts")[gene, ]
      
      # Count cells with expression > 0
      expressed_cells <- sum(expr_values > 0)
      
      # Calculate percentage
      expressed_percent <- (expressed_cells / total_cells) * 100
      
      # Save result
      results[[length(results) + 1]] <- data.frame(
        sample = sample,
        gene = gene,
        expressing_cells = expressed_cells,
        total_cells = total_cells,
        percent_expressing = expressed_percent
      )
    } else {
      warning(paste("Gene", gene, "not found in sample", sample))
    }
  }
}

# Combine all rows into one data.frame
result_df <- do.call(rbind, results)
result_df_unique <- result_df[!duplicated(result_df), ]
# Ensure the sample order matches car_genes
# Step 1: Prepare scRNA result
result_df_unique <- result_df %>%
  distinct(sample, gene, .keep_all = TRUE) %>%
  mutate(source = "scRNA")

# Step 2: Build flow data frame
flow_list <- list()
for (sample in names(car_expression_by_flow)) {
  genes <- car_genes[[sample]]
  flow_percents <- as.numeric(car_expression_by_flow[[sample]])
  
  # Match gene to corresponding flow value
  for (i in seq_along(genes)) {
    flow_list[[length(flow_list) + 1]] <- data.frame(
      sample = sample,
      gene = genes[i],
      expressing_cells = NA,   # not applicable
      total_cells = NA,        # not applicable
      percent_expressing = flow_percents[i],
      source = "Flow"
    )
  }
}

# Combine flow and scRNA data
flow_df <- do.call(rbind, flow_list)
combined_df <- bind_rows(result_df_unique, flow_df)

# Step 3: Generate consistent sample-gene labels and ordering
combined_df$sample_gene <- paste(combined_df$sample, combined_df$gene, sep = " - ")

ordered_combos <- unlist(lapply(names(car_genes), function(sample) {
  paste(sample, car_genes[[sample]], sep = " - ")
}))

combined_df$sample_gene <- factor(combined_df$sample_gene, levels = ordered_combos)
combined_df$sample <- factor(combined_df$sample, levels = names(car_genes))

# Step 4: Plot
flow_car_expression <- ggplot(combined_df, aes(x = sample_gene, y = percent_expressing, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    title = "CAR+ NK cells - scRNA vs. Flow",
    x = "Sample - CAR Gene",
    y = "Percent CAR Expressing Cells",
    fill = "Data Source"
  ) +
  theme_minimal(base_size = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(flow_car_expression)

ggsave(
    filename = file.path(annotation_dir, "Flow_vs_Sc.png"),
    plot = flow_car_expression,
    width = 180,
    height = 180,
    bg = "white",
    dpi = 300,
    units = "mm"
)

message("Finished batch: ", batch_name)
