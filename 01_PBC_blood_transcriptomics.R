###############################################################################
# PROJECT 1 - PBC BLOOD TRANSCRIPTOMICS
# Dataset: GSE119600 (NCBI Gene Expression Omnibus)
# Question: Which blood gene-expression features differ between PBC patients
#           and healthy controls?
#
# IMPORTANT:
# - The GEO series contains several diseases. This script selects only PBC
#   samples and healthy controls, based on the labels in the GEO sample titles.
# - GEO reports quantile-normalized Illumina expression intensities.
#   The script checks their scale and applies log2 transformation if needed.
# - This is an exploratory differential-expression analysis, not a validated
#   diagnostic classifier or a clinical biomarker study.
###############################################################################

# ---------------------------- 0. SETUP ---------------------------------------

# Run this script from its project folder in RStudio.
# Output files will be saved in ./results and ./results/figures.

if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)

# Install required packages only when missing.
cran_packages <- c("ggplot2", "pheatmap")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c("GEOquery", "limma")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

# ----------------------- 1. DOWNLOAD AND IMPORT -----------------------------

# GEOquery downloads the processed Series Matrix file from NCBI GEO.
# The first object is the ExpressionSet for the platform used in this series.
message("Downloading GSE119600 from NCBI GEO...")
gse_list <- getGEO("GSE119600", GSEMatrix = TRUE)

if (length(gse_list) < 1) {
  stop("No ExpressionSet was returned for GSE119600.")
}

eset <- gse_list[[1]]
expr_all <- exprs(eset)
metadata_all <- pData(eset)

message("Expression matrix dimensions: ", nrow(expr_all), " probes x ",
        ncol(expr_all), " samples")
message("Metadata dimensions: ", nrow(metadata_all), " samples x ",
        ncol(metadata_all), " variables")

# ------------------------- 2. SAMPLE ANNOTATION ------------------------------

# In this GEO series, the sample titles contain the group suffix:
#   *_PBC = primary biliary cholangitis
#   *_K   = healthy controls
# Other groups (PSC, CD, UC) are not used in the primary comparison.
if (!"title" %in% colnames(metadata_all)) {
  stop("The GEO metadata do not contain a 'title' column. Inspect the saved\n",
       "metadata file and adapt the group-classification section.")
}

sample_title <- trimws(as.character(metadata_all$title))
group <- rep(NA_character_, length(sample_title))
group[grepl("PBC\\s*$", sample_title, ignore.case = TRUE)] <- "PBC"
group[grepl("(^|[_[:space:]])K\\s*$", sample_title, ignore.case = TRUE)] <- "Control"
metadata_all$group <- group

# Save all metadata, including the group label derived from each GEO title.
write.csv(metadata_all,
          file = "results/GSE119600_all_sample_metadata.csv",
          row.names = TRUE)

cat("\nGroup labels across the complete GEO series:\n")
print(table(metadata_all$group, useNA = "ifany"))

# Verify that the expected cohort sizes have been identified.
# This validation prevents accidental inclusion of PSC/IBD or mislabelled samples.
expected_pbc <- 90
expected_controls <- 47
observed_pbc <- sum(metadata_all$group == "PBC", na.rm = TRUE)
observed_controls <- sum(metadata_all$group == "Control", na.rm = TRUE)

if (observed_pbc != expected_pbc || observed_controls != expected_controls) {
  stop(
    "The sample-title parsing did not identify the expected cohort sizes.\n",
    "Expected: 90 PBC and 47 controls. Found: ", observed_pbc,
    " PBC and ", observed_controls, " controls.\n",
    "Open results/GSE119600_all_sample_metadata.csv, inspect the 'title'\n",
    "column, and adapt the group labels before proceeding."
  )
}

# Align metadata rows to expression-matrix columns using GEO sample IDs.
if (is.null(colnames(expr_all)) || is.null(rownames(metadata_all))) {
  stop("Sample identifiers are missing from the expression matrix or metadata.")
}
metadata_index <- match(colnames(expr_all), rownames(metadata_all))
if (anyNA(metadata_index)) {
  stop("Could not align all expression columns to metadata rows.")
}
metadata_all <- metadata_all[metadata_index, , drop = FALSE]

# Keep only PBC and healthy control samples.
keep_samples <- metadata_all$group %in% c("Control", "PBC")
expr <- expr_all[, keep_samples, drop = FALSE]
metadata <- metadata_all[keep_samples, , drop = FALSE]
metadata$group <- factor(metadata$group, levels = c("Control", "PBC"))

# Final alignment check: each expression column must correspond to the same
# sample as the metadata row in that position.
stopifnot(identical(colnames(expr), rownames(metadata)))
stopifnot(ncol(expr) == nrow(metadata))

cat("\nSelected cohort:\n")
print(table(metadata$group))
write.csv(metadata,
          file = "results/selected_sample_metadata.csv",
          row.names = TRUE)

# --------------------- 3. EXPRESSION SCALE / PREPARATION ---------------------

# GEO sample records describe the submitted values as quantile-normalized
# signal intensities. The values in this dataset are expected to be on an
# intensity scale rather than log2 scale. Inspect the distribution first.
cat("\nExpression-value quantiles before log check:\n")
print(quantile(as.numeric(expr),
               probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
               na.rm = TRUE))

if (any(!is.finite(expr), na.rm = TRUE)) {
  stop("The expression matrix contains non-finite values. Inspect before analysis.")
}
if (min(expr, na.rm = TRUE) <= 0) {
  stop("Expression values include zero or negative values; the log2 step needs review.")
}

# A 99th percentile much higher than 30 indicates unlogged intensity values.
# Log2 is applied once; no second quantile-normalization is performed.
if (as.numeric(quantile(expr, 0.99, na.rm = TRUE)) > 30) {
  message("Applying log2 transformation to quantile-normalized intensities.")
  expr <- log2(expr)
} else {
  message("Values appear to be on a log-like scale; no additional log2 applied.")
}

# Remove probes with no usable measurements.
valid_probes <- rowSums(is.finite(expr)) == ncol(expr)
expr <- expr[valid_probes, , drop = FALSE]

saveRDS(expr, file = "results/expression_matrix_PBC_vs_control.rds")

# ---------------------------- 4. QUALITY CONTROL -----------------------------

# 4A. Distribution of expression values per sample.
pdf("results/figures/01_sample_expression_boxplot.pdf", width = 15, height = 7)
boxplot(expr,
        outline = FALSE,
        las = 2,
        cex.axis = 0.45,
        main = "Expression distribution by sample",
        ylab = "Expression (log2 scale)",
        xlab = "Samples")
dev.off()

# 4B. Principal component analysis (PCA) to visualize global sample variation.
pca <- prcomp(t(expr), center = TRUE, scale. = FALSE)
variance_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = metadata$group
)
write.csv(pca_df, "results/PCA_coordinates.csv", row.names = FALSE)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(
    title = "PCA of PBC and healthy-control blood samples",
    x = paste0("PC1 (", round(variance_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme_classic(base_size = 12)

ggsave("results/figures/02_PCA_PBC_vs_control.pdf",
       plot = p_pca, width = 7, height = 5)

# ----------------------- 5. DIFFERENTIAL EXPRESSION --------------------------

# limma fits a linear model for each probe and uses empirical Bayes moderation.
# The contrast is PBC minus Control, so positive logFC means higher expression
# in PBC, and negative logFC means lower expression in PBC.
design <- model.matrix(~ 0 + group, data = metadata)
colnames(design) <- levels(metadata$group)

fit <- lmFit(expr, design)
contrast_matrix <- makeContrasts(PBC - Control, levels = design)
fit_contrasted <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrasted)

deg_table <- topTable(
  fit_ebayes,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
deg_table$probe_id <- rownames(deg_table)

# Multiple-testing correction: Benjamini-Hochberg false-discovery rate (FDR).
# A |log2FC| >= 1 threshold corresponds to at least a two-fold expression change.
deg_table$significance <- "Not significant"
deg_table$significance[deg_table$adj.P.Val < 0.05 & deg_table$logFC >= 1] <- "Higher in PBC"
deg_table$significance[deg_table$adj.P.Val < 0.05 & deg_table$logFC <= -1] <- "Lower in PBC"

deg_table <- deg_table[, c("probe_id", setdiff(colnames(deg_table), "probe_id"))]
write.csv(deg_table,
          file = "results/differential_expression_all_probes.csv",
          row.names = FALSE)

significant_table <- subset(deg_table, adj.P.Val < 0.05 & abs(logFC) >= 1)
write.csv(significant_table,
          file = "results/differentially_expressed_probes_FDR05_logFC1.csv",
          row.names = FALSE)

# Volcano plot.
deg_table$minus_log10_FDR <- -log10(pmax(deg_table$adj.P.Val, .Machine$double.xmin))
p_volcano <- ggplot(deg_table,
                    aes(x = logFC, y = minus_log10_FDR, color = significance)) +
  geom_point(alpha = 0.65, size = 1.2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", size = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.35) +
  labs(
    title = "Differential expression: PBC vs healthy controls",
    x = "log2 fold change (PBC - Control)",
    y = "-log10 adjusted p-value",
    color = "Result"
  ) +
  theme_classic(base_size = 12)

ggsave("results/figures/03_volcano_plot.pdf",
       plot = p_volcano, width = 7, height = 5)

# ---------------------------- 6. PROBE ANNOTATION ----------------------------

# Illumina arrays measure probes. The platform table links probe IDs to genes.
# We retrieve the official GEO platform annotation and look for a gene-symbol
# column. If it is not available, the probe-level analysis remains valid and
# the annotation step can be reviewed separately.
message("Retrieving Illumina HumanHT-12 v4 platform annotation...")
annotation_ok <- FALSE
try({
  gpl <- getGEO("GPL10558")
  gpl_table <- Table(gpl)
  id_candidates <- grep("^ID$", names(gpl_table), ignore.case = TRUE, value = TRUE)
  symbol_candidates <- grep("symbol", names(gpl_table), ignore.case = TRUE, value = TRUE)

  if (length(id_candidates) >= 1 && length(symbol_candidates) >= 1) {
    id_col <- id_candidates[1]
    symbol_col <- symbol_candidates[1]
    annotation <- data.frame(
      probe_id = as.character(gpl_table[[id_col]]),
      gene_symbol = as.character(gpl_table[[symbol_col]]),
      stringsAsFactors = FALSE
    )
    annotation$gene_symbol <- trimws(annotation$gene_symbol)
    annotation$gene_symbol[annotation$gene_symbol %in% c("", "---", "NA")] <- NA
    annotation$gene_symbol <- sub("\\s*///.*$", "", annotation$gene_symbol)

    deg_table$gene_symbol <- annotation$gene_symbol[
      match(deg_table$probe_id, annotation$probe_id)
    ]
    write.csv(annotation, "results/GPL10558_probe_annotation.csv", row.names = FALSE)
    write.csv(deg_table,
              "results/differential_expression_annotated.csv",
              row.names = FALSE)
    write.csv(subset(deg_table, adj.P.Val < 0.05 & abs(logFC) >= 1),
              "results/differentially_expressed_genes_FDR05_logFC1.csv",
              row.names = FALSE)
    annotation_ok <- TRUE
  }
}, silent = TRUE)

if (!annotation_ok) {
  message("Automatic gene-symbol annotation was not completed. Probe-level results were saved.")
}

# ------------------------------- 7. HEATMAP ----------------------------------

# Display the top 30 probes by adjusted p-value. If fewer than 30 probes are
# significant, this remains an exploratory visualization of the top-ranked set.
top_probe_ids <- head(deg_table$probe_id, min(30, nrow(deg_table)))
heatmap_matrix <- expr[top_probe_ids, , drop = FALSE]
heatmap_matrix_z <- t(scale(t(heatmap_matrix)))
heatmap_matrix_z[!is.finite(heatmap_matrix_z)] <- 0

sample_annotation <- data.frame(Group = metadata$group)
rownames(sample_annotation) <- rownames(metadata)

pdf("results/figures/04_top30_probe_heatmap.pdf", width = 10, height = 9)
pheatmap(
  heatmap_matrix_z,
  annotation_col = sample_annotation,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Top 30 probes ranked by adjusted p-value"
)
dev.off()

# ------------------------- 8. OPTIONAL GO ENRICHMENT -------------------------

# This section is optional because it requires additional Bioconductor packages.
# It runs only if clusterProfiler and org.Hs.eg.db are already installed and
# a gene-symbol annotation was obtained above.
if (annotation_ok &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {

  significant_symbols <- unique(na.omit(
    deg_table$gene_symbol[deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) >= 1]
  ))

  if (length(significant_symbols) >= 10) {
    suppressPackageStartupMessages({
      library(clusterProfiler)
      library(org.Hs.eg.db)
    })

    enrichment <- tryCatch({
      ids <- bitr(
        significant_symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      )
      enrichGO(
        gene = unique(ids$ENTREZID),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
      )
    }, error = function(e) {
      message("GO enrichment could not be completed: ", conditionMessage(e))
      NULL
    })

    if (!is.null(enrichment) && nrow(as.data.frame(enrichment)) > 0) {
      write.csv(as.data.frame(enrichment),
                "results/GO_Biological_Process_enrichment.csv",
                row.names = FALSE)
      pdf("results/figures/05_GO_enrichment_dotplot.pdf", width = 9, height = 6)
      print(dotplot(enrichment, showCategory = 15) +
              ggtitle("GO Biological Process enrichment"))
      dev.off()
    }
  } else {
    message("Fewer than 10 significant annotated genes: GO enrichment skipped.")
  }
} else {
  message("Optional GO enrichment skipped (annotation or optional packages missing).")
}

# ---------------------------- 9. PROJECT SUMMARY -----------------------------

n_sig_fdr <- sum(deg_table$adj.P.Val < 0.05, na.rm = TRUE)
n_sig_strict <- nrow(significant_table)

summary_lines <- c(
  "PROJECT 1 - PBC BLOOD TRANSCRIPTOMICS",
  "Dataset: GSE119600",
  "Platform: Illumina HumanHT-12 v4 expression beadchip (GPL10558)",
  "Tissue: whole blood",
  paste("PBC samples:", sum(metadata$group == "PBC")),
  paste("Healthy controls:", sum(metadata$group == "Control")),
  paste("Probes analysed:", nrow(expr)),
  paste("Probes with FDR < 0.05:", n_sig_fdr),
  paste("Probes with FDR < 0.05 and |log2FC| >= 1:", n_sig_strict),
  "Contrast: PBC minus Control (positive logFC = higher expression in PBC)",
  "Statistical method: limma with empirical Bayes moderation",
  "Multiple-testing correction: Benjamini-Hochberg",
  "",
  "Interpretation note:",
  "This is an exploratory comparison of blood transcriptomes. It does not",
  "demonstrate diagnostic performance, causality, or clinical utility."
)
writeLines(summary_lines, "results/analysis_summary.txt")

cat("\nAnalysis finished. Review the files in the results/ folder.\n")
cat("Key outputs: PCA, volcano plot, heatmap, differential-expression tables,\n")
cat("sample metadata, and analysis_summary.txt.\n")
