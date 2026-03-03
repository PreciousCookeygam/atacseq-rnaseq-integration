# =============================================================================
# analysis/02_pi_scores.R — π-score concordance plots & concordance stats
# =============================================================================
#
# Produces:
#   1. π-score scatter — All TCGA vs P0
#   2. Concordance stats — All TCGA
#   3. π-score scatter — LumP vs P0
#   4. Concordance stats — LumP
#   5. Exports concordant gene lists
#
# Run:
#   source(here::here("analysis", "02_pi_scores.R"))
# =============================================================================

# --- Load modules -------------------------------------------------------------
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "load_data.R"))
source(here::here("R", "pi_score.R"))
source(here::here("R", "concordance.R"))

# --- Ensure output directory exists -------------------------------------------
ensure_output_dir(OUTPUT_DIR)

# --- Load data (once) ---------------------------------------------------------
raw_tcga <- load_annotated_results(FILE_TCGA)
raw_lump <- load_annotated_results(FILE_LUMP)

combined_tcga <- prepare_combined_data(raw_tcga)
combined_lump <- prepare_combined_data(raw_lump)

# =============================================================================
# 1) π-score plot — All TCGA vs P0
# =============================================================================
pi_tcga <- compute_pi_scores(combined_tcga)

plot_pi_scores(
  df          = pi_tcga,
  title       = expression("π-Score Concordance: RNA-seq vs ATAC-seq (All TCGA)"),
  xy_limit    = 7,
  output_file = "pi_score_allTCGA_vs_P0.png"
)

# =============================================================================
# 2) Concordance statistics — All TCGA
# =============================================================================
conc_tcga <- compute_concordance(combined_tcga, label = "All TCGA")

# =============================================================================
# 3) π-score plot — LumP vs P0
# =============================================================================
pi_lump <- compute_pi_scores(combined_lump)

plot_pi_scores(
  df          = pi_lump,
  title       = expression("π-Score Concordance: RNA-seq vs ATAC-seq (LumP)"),
  xy_limit    = 6,
  output_file = "pi_score_LumP_vs_P0.png"
)

# =============================================================================
# 4) Concordance statistics — LumP
# =============================================================================
conc_lump <- compute_concordance(combined_lump, label = "LumP")

# =============================================================================
# 5) Compare concordant gene sets & export
# =============================================================================
genes_tcga <- conc_tcga$concordant_genes
genes_lump <- conc_lump$concordant_genes

cat(sprintf(
  "\nConcordant genes: TCGA=%d, LumP=%d, shared=%d, identical=%s\n",
  length(genes_tcga),
  length(genes_lump),
  length(intersect(genes_tcga, genes_lump)),
  setequal(genes_tcga, genes_lump)
))

# Export TCGA concordant gene list
write.table(
  genes_tcga,
  file.path(OUTPUT_DIR, "concordant_core_genes_tcga.tsv"),
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

message("\n[02_pi_scores] π-score plots and concordance analysis complete.")
