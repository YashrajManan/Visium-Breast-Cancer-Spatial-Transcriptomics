# Visium Breast Cancer Spatial Transcriptomics

Spatial transcriptomics analysis of a real 10x Visium human breast cancer section — mapping gene expression
back onto real tissue coordinates to find spatially coherent expression domains, genes whose expression is
non-randomly organized across physical space, and which tissue domains sit adjacent to each other more than
chance would predict.

## Aim

Find out whether a real breast tumor section's transcriptional heterogeneity is spatially organized, and if
so, which genes drive that organization and which tissue domains are physically adjacent to which others —
questions dissociated single-cell RNA-seq cannot answer, since dissociation destroys spatial position before
sequencing.

## Objective

Fetch a real public 10x Visium dataset; run standard QC and expression-based clustering to define candidate
spatial domains; formally test which genes are spatially organized (Moran's I) and which domains are
spatially adjacent more than random chance predicts (a permutation-based neighborhood enrichment test);
visualize both directly on the real tissue image.

## Data fetch

`V1_Breast_Cancer_Block_A_Section_1` — 10x Genomics' real, public Visium Spatial Gene Expression dataset,
fresh-frozen human invasive ductal carcinoma tissue, released under CC BY 4.0. Fetched live and
programmatically via `squidpy.datasets.visium(sample_id="V1_Breast_Cancer_Block_A_Section_1")`, which
downloads the real Space Ranger output (filtered feature-barcode matrix, spot coordinates, scale factors,
and the real H&E tissue image) directly from 10x Genomics' dataset server.

## Data describe

3798 real spots under tissue, 36601 genes measured per spot before filtering. Mean 21,815 UMIs and 5,622
genes detected per spot; mitochondrial fraction ranged 1.1%–14.4% across spots (mean 4.0%), well under the
20% quality threshold used — every one of the 3798 spots passed QC. After filtering out genes detected in
fewer than 3 spots, 22,240 genes remained for downstream analysis.

## Methods / Workflow — what we did

Standard scRNA-seq preprocessing (total-count normalization, log-transform, top-2000 highly variable genes,
50-component PCA), then two separate graphs built on the same data: an expression-based k-NN graph (PCA
space) used for Leiden clustering into candidate spatial domains, and a real spatial neighbors graph (from
real spot pixel coordinates, Visium's hexagonal grid geometry) used for every genuinely spatial statistic.
Clustering was deliberately run on expression alone, spatially blind, so that neighborhood enrichment
afterward is an independent spatial test rather than something built into the clustering objective. Moran's
I (a classical spatial-autocorrelation statistic) was computed per gene on the real spatial graph to
formally rank spatially variable genes. Neighborhood enrichment (a 1000-permutation test) measured whether
each pair of Leiden domains is spatially adjacent more or less than random domain-label shuffles on the real
graph would predict.

**R/Seurat side (cross-language comparison).** The same real Space Ranger output (saved once from the
Python fetch) was independently re-analyzed in R: `SCTransform` normalization (Seurat's own documented
preference for Visium data, in place of simple log-normalization) followed by PCA, then Seurat's default
Louvain clustering (in place of Leiden — Leiden isn't bundled with Seurat and pulling it in via `reticulate`
would just add friction without changing the conceptual comparison), and `FindSpatiallyVariableFeatures(...,
selection.method = "moransi")` as the R-side equivalent of squidpy's Moran's I. Seurat has no native
neighborhood-enrichment permutation test (a real, disclosed capability gap — see `WORKFLOW.md`), so that
statistic is Python-only.

## Results

12 real Leiden clusters were found (sizes 100–544 spots, summing exactly to all 3798 spots). The top
Moran's-I-ranked spatially variable genes were MALAT1 (I = 0.751), followed by a run of immunoglobulin
genes — IGHG3 (0.700), IGLC2 (0.691), IGKC (0.682), IGHG4 (0.628), IGHG1 (0.623) — then CPB1 (0.672), APOE
(0.598), C3 (0.581), and CXCL14 (0.561), all statistically significant (FDR-corrected p < 0.02). The
neighborhood-enrichment z-score matrix (12×12) was strongly diagonal-dominant — every domain's strongest
enrichment was with itself, meaning spots from the same expression-defined cluster are real, physical
neighbors far more often than random chance predicts. Full tables: `results/spatially_variable_genes.csv`,
`results/spot_domain_assignments.csv`, `results/nhood_enrichment_zscores.csv`. Figures:
`figures/clusters_on_tissue.png`, `figures/top_spatial_genes.png`, `figures/nhood_enrichment.png`.

**R/Seurat side.** Louvain clustering on the same real, QC-filtered 3798 spots found 20 clusters (sizes
33–520 spots, also summing exactly to 3798) — more, smaller clusters than Python's Leiden run, an expected
real difference between the two default community-detection algorithms at the same resolution rather than a
disagreement about the underlying biology. The top Moran's-I-ranked spatially variable genes (via Seurat's
`FindSpatiallyVariableFeatures`, on the top 1000 SCT-variable genes) were CRISP3, CXCL14, TTLL12, KRT8,
SLITRK6, GFRA1, MGP, S100A11, CCND1, and KRT18. Four genes — CRISP3, CXCL14, CPB1, and SLITRK6 — appear in
*both* languages' top-15 spatially variable gene lists, a genuine cross-tool agreement signal given the two
pipelines differ in normalization (log vs. SCTransform) and gene universe (all genes vs. top-1000
SCT-variable genes). Full tables: `results_R/spatially_variable_genes_R.csv`,
`results_R/spot_domain_assignments_R.csv`. Figures: `results_R/clusters_on_tissue_R.png`,
`results_R/top_spatial_genes_R.png`.

## Biology interpretation of results

The diagonal-dominant enrichment matrix is itself the headline finding: this tissue's expression-only
clustering recovered domains that are also spatially coherent, meaning the dominant axis of transcriptional
variation in this section tracks real physical tissue organization rather than cutting across it randomly —
consistent with how real tumor biology works, where cells sharing an expression state (tumor, stroma, immune
infiltrate) also tend to occupy the same physical region. The specific top spatially variable genes are
biologically coherent, not arbitrary: a cluster of immunoglobulin genes (IGHG3, IGLC2, IGKC, IGHG4, IGHG1)
ranking among the most spatially organized genes in the whole dataset is a real, recognizable signature of
plasma-cell/immune infiltrate concentrated in specific tissue regions rather than spread evenly throughout —
exactly the kind of spatially-localized immune response that dissociated scRNA-seq would average away. MALAT1
topping the list is a common, expected finding (it's a near-ubiquitously highly expressed transcript whose
raw signal often tracks tissue density/quality as much as biology, worth treating as a partial technical
confound rather than pure biological signal). CPB1, APOE, and C3 are all real markers plausible for breast
epithelial and stromal/complement activity respectively, consistent with the expected cell-type composition
of invasive ductal carcinoma tissue.

## Learning through project

The real, non-trivial bug in this build wasn't in the analysis logic — it was in where results landed.
`squidpy`/`scanpy`'s plotting functions default to writing `save=` figures into `./figures/` (scanpy's own
`sc.settings.figdir` default), not into whatever output folder a script happens to `os.makedirs()` itself;
the CSV exports (plain `pandas.to_csv` calls) were unaffected since they don't go through that setting at
all. Caught by checking the real local filesystem after a run rather than assuming a `save="name.png"`
argument writes wherever a results directory variable happens to be set — fixed by explicitly setting
`sc.settings.figdir` once, in the same cell the results directory is created. A small thing, but a real
example of why checking what a library's default behavior actually does beats assuming a parameter name
implies a destination.

## Caveats

A single Visium capture area, one patient, one tissue section — real, genuine spatial biology, but not a
multi-patient cohort, so findings here describe this specific section's organization, not a general claim
about breast cancer spatial architecture. Each Visium spot captures RNA from several overlapping real cells,
not one cell, so "domain" and "spatially variable gene" describe local-tissue-neighborhood patterns, not
single-cell states. MALAT1's top ranking is flagged above as a partial technical-signal confound worth
weighing against the more clearly cell-type-specific genes below it.
