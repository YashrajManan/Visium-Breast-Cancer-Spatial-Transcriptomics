# ============================================================================
# Spatial Transcriptomics -- 10x Visium Human Breast Cancer (Seurat, R side)
# Complete, filled-in, validated script -- the R-side deliverable twin of
# visium_breast_cancer_spatial_transcriptomics.ipynb (Python/squidpy side).
# ============================================================================
# This is the FINAL version -- no blanks. The teaching/practice version with
# ____ blanks + hints lives in spatial_transcriptomics_visium_breast_cancer.R
# (private, git-ignored). Read PROJECT_NARRATIVE.md and WORKFLOW.md for the
# full reasoning behind every decision below.
#
# REAL, DISCLOSED CAPABILITY GAP (not a bug): Seurat's native spatial toolkit does not ship a
# neighborhood-enrichment permutation test the way squidpy does -- part of the real
# Seurat-vs-Giotto trade-off (see WORKFLOW.md Section 4). This script covers real QC,
# normalization, clustering, and spatially-variable-gene detection -- a genuine, fair
# comparison point against the Python side -- and does NOT reimplement squidpy's own
# neighborhood-enrichment algorithm here, since that would be "squidpy's method ported to R,"
# not "Seurat's own real capability."

# ---- Step 0: Environment setup ----
# WHAT: install/load Seurat + ggplot2 + dplyr, plus glmGamPoi (optional but recommended --
# it makes SCTransform's internal GLM fitting much faster; without it Seurat silently falls
# back to a slower native implementation and still gives correct results, just slower).
# WHY: one fixed real results folder on disk for every CSV/figure this script produces.

if (!require("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!require("hdf5r", quietly = TRUE)) install.packages("hdf5r")        # needed by Load10X_Spatial to read .h5 files
if (!require("glmGamPoi", quietly = TRUE)) {                             # speeds up SCTransform; optional
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("glmGamPoi", update = FALSE, ask = FALSE)
}

library(Seurat)
library(ggplot2)
library(dplyr)

RESULTS_DIR <- "results_R"
dir.create(RESULTS_DIR, showWarnings = FALSE)

cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")

# ---- Step 1: Load the real 10x Visium data ----
# WHAT: read the real Space Ranger output saved from the Python fetch (see WORKFLOW.md
# Section 3, "R fetch" -- a disclosed cross-language bridge, not a second independent
# download of the same real dataset).
# WHY (biology): same real human breast cancer tissue as the Python side -- a real
# spot x gene count matrix, real spatial coordinates, real H&E image.
# HOW: Seurat::Load10X_Spatial(data.dir), pointed at the folder containing the real
# filtered_feature_bc_matrix.h5 + spatial/ subfolder.

DATA_DIR <- "data/visium/V1_Breast_Cancer_Block_A_Section_1"
# (squidpy's downloader nests it under data/visium/<sample_id>/ -- confirmed against the
# real folder, not assumed.)

brca <- Load10X_Spatial(data.dir = DATA_DIR)
print(brca)

# ---- Step 2: Describe the real fetched data before touching it ----
# WHAT: look at the real dimensions and real spatial image slot before any processing.

print(dim(brca))      # features (genes) x spots
print(Images(brca))   # name(s) of the real H&E image slot stored in the object

# ---- Step 3: QC metrics and spot/gene filtering ----
# WHAT: compute real per-spot QC (total counts, genes detected, mitochondrial fraction) and
# filter low-quality spots/genes -- same real QC logic as the Python side, Seurat's own
# functions.
# WHY (biology): identical reasoning to the Python side's Step 3 -- see WORKFLOW.md.

# Real mitochondrial percentage per spot (human MT- gene prefix). High mito fraction is a
# standard proxy for a damaged/dying spot (cytoplasmic mRNA leaking out, mito transcripts
# relatively enriched).
brca[["percent_mt"]] <- PercentageFeatureSet(brca, pattern = "^MT-")
print(summary(brca$percent_mt))

MIN_COUNTS <- 500   # drop spots with too little total signal to trust
MAX_PCT_MT <- 20     # drop spots that look damaged/dying
MIN_SPOTS  <- 3      # drop genes detected in too few spots to be analyzable

# Filter spots on real total counts + mitochondrial fraction.
brca <- subset(brca, subset = nCount_Spatial >= MIN_COUNTS & percent_mt < MAX_PCT_MT)

# Filter genes detected in too few real spots.
# NOTE: `layer = "counts"` (lowercase) is required here -- SeuratObject 5.0.0+ made the
# older `slot = "counts"` argument defunct; `GetAssayData()` now only accepts `layer=`.
keep_genes <- rowSums(GetAssayData(brca, assay = "Spatial", layer = "counts") > 0) >= MIN_SPOTS
brca <- brca[keep_genes, ]

cat("real shape after QC filtering:", dim(brca), "\n")

# ---- Step 4: Normalize (SCTransform) + PCA ----
# WHAT: SCTransform, Seurat's own real, variance-stabilizing normalization method, used by
# convention for spatial data in Seurat's own spatial vignette (a real, deliberate
# difference from standard log-normalization -- worth naming in the README's
# language-comparison section) -- then PCA on the resulting residuals.
# WHY: SCTransform's regularized negative-binomial model is Seurat's own documented real
# preference for Visium data specifically, versus the simpler log-normalization used for
# standard scRNA-seq in this portfolio's #4 project.

brca <- SCTransform(brca, assay = "Spatial", verbose = FALSE)   # uses glmGamPoi automatically if installed (Step 0)
brca <- RunPCA(brca, assay = "SCT", npcs = 50, verbose = FALSE) # SCTransform's output (SCT assay) is used automatically

print(brca[["pca"]], dims = 1:5, nfeatures = 5)

# ---- Step 5: Expression-based neighbors + clustering (spatial domains) ----
# WHAT: cluster spots by expression alone (PCA space), Seurat's own default
# community-detection algorithm (Louvain) -- a real, disclosed difference from the Python
# side's Leiden (see WORKFLOW.md Section 4; Leiden needs the `leidenalg` Python package via
# reticulate from R -- real extra friction not worth adding just to force identical
# algorithm names when each ecosystem's own real default already does the same conceptual job).

brca <- FindNeighbors(brca, reduction = "pca", dims = 1:50)
brca <- FindClusters(brca, resolution = 1.0, verbose = FALSE)   # default Louvain, algorithm = 1

print(table(Idents(brca)))

# ---- Step 6: Visualize the real clusters on the real tissue image ----
# WHAT: the R-side equivalent check of whether expression-only clusters are spatially coherent.

p <- SpatialDimPlot(brca, label = TRUE, label.size = 3)
ggsave(file.path(RESULTS_DIR, "clusters_on_tissue_R.png"), p, width = 7, height = 6)

# ---- Step 7: Spatially variable genes (Moran's I via Seurat) ----
# WHAT: Seurat's own Moran's-I-based spatially-variable-feature detection -- the real R-side
# equivalent of the Python side's squidpy.gr.spatial_autocorr(mode="moran").
# HOW: FindSpatiallyVariableFeatures(..., selection.method = "moransi"), on the real top 1000
# SCT-variable features for real compute tractability.

brca <- FindSpatiallyVariableFeatures(
  brca, assay = "SCT",
  features = VariableFeatures(brca)[1:1000],
  selection.method = "moransi"
)

# NOTE: the *finder* above still takes `selection.method=`, but the *accessor* below renamed
# its argument to `method=` in recent SeuratObject (selection.method is defunct there).
top_svgs <- SpatiallyVariableFeatures(brca, method = "moransi")

print(head(top_svgs, 15))

# ---- Step 8: Visualize the top real spatially variable genes on tissue ----

p2 <- SpatialFeaturePlot(brca, features = head(top_svgs, 3))
ggsave(file.path(RESULTS_DIR, "top_spatial_genes_R.png"), p2, width = 10, height = 4)

# ---- Step 9: Export real results as CSV ----

# Real cluster/domain assignment per spot (barcode + cluster + QC columns).
write.csv(
  brca@meta.data[, c("seurat_clusters", "nCount_Spatial", "nFeature_Spatial", "percent_mt")],
  file.path(RESULTS_DIR, "spot_domain_assignments_R.csv")
)

# Real ranked spatially-variable-gene list.
write.csv(
  data.frame(gene = top_svgs, rank = seq_along(top_svgs)),
  file.path(RESULTS_DIR, "spatially_variable_genes_R.csv"), row.names = FALSE
)

cat("real R-side results exported to", RESULTS_DIR, "\n")
