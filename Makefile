.PHONY: all 01_read 02_qc 03_anno 04_integration figures

all: 04_integration

04_integration: 03_anno
	Rscript pipeline/03_Integration.R

03_anno: 02_qc
	Rscript pipeline/02_Anno.R

02_qc: 01_read
	Rscript pipeline/01_QC.R

01_read:
	Rscript pipeline/00_ReadSeurat.R

figures:
	quarto render notebooks/Figures.qmd --execute-debug