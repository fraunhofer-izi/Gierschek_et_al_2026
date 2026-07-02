# 2026_Gierschek_et_al_CD123_CD33_CAR_NK

Code to reproduce the analysis and figures from the publication "Comparison with CD33-CAR-NK cells reveals large-scale CD123-CAR-NK cells as highly promising strategy to combat AML" (manuscript not yet released)

## Repository Structure

```text
.
├── assets/         # Figures, logos, documentation assets
├── config/         # Configuration files
├── data/           # Symlinks to input .pxl files
├── environment.yml # Conda/Mamba environment specification
├── environment_frozen.yml # Specific frozen Conda/Mamba environment specification
├── Makefile        # Project setup utilities
├── notebooks/      # Analysis notebooks
├── packages/       # Custom R packages and scripts
├── pipeline/       # Pipeline submodules
├── README.md       # README 
├── results/        # Symlinks to analysis output directories
└── scripts/        # Utility scripts
```

# Reproducibility

This repository supports two reproducible workflows depending on your objective.

## Option 1. Generate Publication Figures 

If your goal is to reproduce the figures presented in the manuscript, you can start directly from the processed integrated Seurat object.

This workflow minimizes variability introduced during preprocessing and integration, such as stochastic initialization, graph-based clustering, and differences between software versions or computing environments.


### 1. Clone the repository

```bash
git clone <repository-url>
cd Gierschek_et_al_2026
```

### 2. Create the environment

```bash
mamba env create -f environment.yml
mamba activate 2026_Gierschek_et_al_CD123_CD33_CAR_NK
```

### 3. Download the processed Seurat object

Download the processed integrated Seurat object from the GEO repository and place it in:

```text
results/processed_seurat/
```

### 4. Generate the figures

```bash
make figures
```

This workflow reproduces the publication figures directly from the processed Seurat object used for the manuscript.

## Option 2. Full Analysis Pipeline (from raw data)

Use this workflow to reproduce the complete single-cell analysis starting from the processed CellRanger output.

### 1. Clone the repository

```bash
git clone <repository-url>
cd Gierschek_et_al_2026
```

### 2. Create the environment

```bash
mamba env create -f environment.yml
mamba activate 2026_Gierschek_et_al_CD123_CD33_CAR_NK
```

### 3. Download CellRanger output

Download the processed CellRanger output from the GEO repository and organize each sample as:

```text
data/per_sample_outs/NTNK-OCIAML2/count/sample_raw_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

### 4. Run the complete analysis

```bash
make all
```

This workflow reproduces all preprocessing steps, including quality control, annotation, integration, and generation of the processed Seurat object.

---

The pipeline consists of the following stages:

| Step                    | Script                 | Description                       |
| ----------------------- | ---------------------- | --------------------------------- |
| Read CellRanger output | `00_ReadSeurat.R` | Create Seurat object       |
| Quality control         | `01_QC.R`               | Cell and gene filtering           |
| Annotation              | `02_Anno.R`          | Annotation of cells  |
| Integration              | `03_Integration.R`       | Integration of NK cells across samples |

Intermediate results are stored in:

```text
results/
├── qc/
├── annotation/
├── integration/
├── processed_seurat/
```

---

### 5. Publication Figures

Publication figures are generated using:

```bash
make figures
```

Results are stored in:

```text
results/
├── publication_figures
```

Quarto version: 1.8.26

# Configuration Files

Configs (path to data, results, and qc filters) are stored in:

```text
config/
├── configs.yml
```

# Citation

If you use this repository, please cite:

> Gierscheck et al. – *In-depth comparative characterization identifies scalable CD123-CAR-NK cells as promising strategy to combat AML*

---

# Contact

For questions regarding the analysis pipeline, please open an issue in the repository.

---

# License

Copyright 2025 Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V.

Licensed under the GPL-3.0. You may obtain a copy of the License in the LICENSE file.