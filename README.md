# Project 1 — PBC blood transcriptomics

## Biological question

Which gene-expression features differ between patients with primary biliary cholangitis (PBC) and healthy controls in whole blood?

## Dataset

- **GEO accession:** GSE119600
- **Organism:** *Homo sapiens*
- **Material:** whole-blood RNA
- **Technology:** Illumina HumanHT-12 v4 expression beadchip (GPL10558)
- **Main comparison:** 90 PBC samples vs 47 healthy controls
- **Other cohorts in the original GEO series:** primary sclerosing cholangitis (PSC), Crohn's disease (CD), and ulcerative colitis (UC). They are excluded from the primary comparison.
- **Preprocessing in GEO:** quantile-normalized signal intensities (lumi workflow reported in sample records). The R script inspects the expression scale and applies log2 transformation when needed; it does not quantile-normalize the matrix a second time.

Official sources:
- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119600
- GEO platform: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558
- Example GEO sample processing record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3378164

## Analysis workflow

1. Download the GEO Series Matrix with `GEOquery`.
2. Inspect expression matrix and sample metadata.
3. Identify PBC and healthy-control samples using the GEO sample labels.
4. Align metadata rows to expression-matrix columns by GEO sample accession.
5. Inspect expression-value scale; apply log2 transformation if needed.
6. Perform basic quality control (sample distribution plot and PCA).
7. Fit a differential-expression model with `limma`.
8. Correct p-values for multiple testing using Benjamini–Hochberg FDR.
9. Annotate Illumina probe IDs using the GEO platform record.
10. Visualize results (volcano plot and heatmap) and optionally run GO Biological Process enrichment.

## Statistical design

The primary contrast is **PBC minus Control**. Positive log2 fold changes indicate higher expression in the PBC group; negative values indicate lower expression in PBC.

The script reports:
- all probes ranked by differential-expression evidence;
- probes with FDR < 0.05;
- a stricter exploratory set with FDR < 0.05 and absolute log2 fold change >= 1;
- PCA and heatmap visualizations;
- optional GO enrichment if the relevant Bioconductor packages and probe annotation are available.

## What this project can and cannot demonstrate

It can demonstrate familiarity with a standard bulk-transcriptomics workflow, public-dataset retrieval, sample metadata handling, quality control, differential-expression analysis, multiple-testing correction, visualization, and biological interpretation.

It **cannot** establish that the resulting genes are validated diagnostic biomarkers, that they are causal in PBC, or that the model can distinguish PBC from AIH in clinical practice. The primary comparison is PBC vs healthy controls, not PBC vs AIH or PBC-AIH variant syndrome. The data come from whole blood, not liver tissue.

## Files produced by the R script

The script creates a `results/` folder containing:

- full GEO sample metadata;
- selected PBC/control metadata;
- an RDS expression matrix;
- differential-expression tables;
- probe annotation (if retrieved successfully);
- PCA coordinates;
- a concise analysis summary;
- PDF figures for sample distributions, PCA, volcano plot, and heatmap;
- optional GO enrichment table and dotplot.

## How to run

1. Open RStudio.
2. Open `01_PBC_blood_transcriptomics.R`.
3. Set the working directory to this project folder (in RStudio: Session → Set Working Directory → To Source File Location).
4. Run the whole script, or use Source.
5. Wait for the GEO download and package setup to finish.
6. Review the outputs in `results/`.

The script is designed to stop if it cannot confidently identify exactly 90 PBC samples and 47 controls. If that happens, inspect `results/GSE119600_all_sample_metadata.csv` before changing any group labels.
