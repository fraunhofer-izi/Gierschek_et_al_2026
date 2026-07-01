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

# Reproducability 

## 1. Clone Repository

```
git clone <repository-url>
cd Gierschek_et_al_2026
```

Requirements:
- conda/miniforge installed
- mamba installed

## 2. Setup

```bash
mamba env create -f environment.yml
mamba activate 2026_Gierschek_et_al_CD123_CD33_CAR_NK
```

## 3. Download data 

```
mkdir -p data/per_sample_outs
cd data/per_sample_outs
```

Download data from GEO repository.
After extraction, organize the processed CellRanger output for each sample as:

```
data/per_sample_outs/NTNK-OCIAML2/count/sample_raw_feature_bc_matrix/barcodes.tsv.gz
data/per_sample_outs/NTNK-OCIAML2/count/sample_raw_feature_bc_matrix/features.tsv
data/per_sample_outs/NTNK-OCIAML2/count/sample_raw_feature_bc_matrix/matrix.mtx
```

## 4. Run Pipeline


The full analysis pipeline can be executed using:

```bash
make all
```

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

# 4. Publication Figures

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