# =============================================================================
# analysis/01_volcanos.R — Differential accessibility & expression volcanos
# =============================================================================
#
# Produces four volcano plots:
#   1. ATAC-seq — All TCGA vs P0
#   2. ATAC-seq — LumP vs P0
#   3. RNA-seq  — All TCGA vs P0
#   4. RNA-seq  — LumP vs P0
#
# Each file is loaded once (via the caching layer) and passed to the shared
# `make_volcano()` function from R/volcano.R.
#
# Run:
#   source(here::here("analysis", "01_volcanos.R"))
# =============================================================================

# --- Load modules -------------------------------------------------------------
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "load_data.R"))
source(here::here("R", "volcano.R"))

# --- Ensure output directory exists -------------------------------------------
ensure_output_dir(OUTPUT_DIR)

# --- Load data (once) ---------------------------------------------------------
raw_tcga <- load_annotated_results(FILE_TCGA)
raw_lump <- load_annotated_results(FILE_LUMP)

# =============================================================================
# 1) ATAC-seq Volcano — All TCGA vs P0
# =============================================================================
atac_tcga <- prepare_atac_data(raw_tcga)

make_volcano(
  df            = atac_tcga,
  fc_col        = "log2FoldChange",
  pval_col      = "padj",
  title         = "Differential Accessibility of (All TCGA) vs Normal Bladder",
  colors        = VOLCANO_COLORS_ATAC,
  x_limits      = c(-10, 10),
  y_limits      = c(0, 10),
  right_label   = "TCGA",
  left_label    = "P0",
  right_color   = "red",
  left_color    = "blue",
  label_x_pos   = c(-8, 8),
  label_y_pos   = 0.5,
  output_file   = "volcano_ATAC_allTCGA_vs_P0.png"
)

# =============================================================================
# 2) ATAC-seq Volcano — LumP vs P0
# =============================================================================
atac_lump <- prepare_atac_data(raw_lump)

make_volcano(
  df            = atac_lump,
  fc_col        = "log2FoldChange",
  pval_col      = "padj",
  title         = "Differential Accessibility of (LumP) vs Normal Bladder",
  colors        = VOLCANO_COLORS_ATAC,
  x_limits      = c(-10, 10),
  y_limits      = c(0, 10),
  right_label   = "LumP",
  left_label    = "P0",
  right_color   = "red",
  left_color    = "blue",
  label_x_pos   = c(-8, 8),
  label_y_pos   = 0.5,
  output_file   = "volcano_ATAC_LumP_vs_P0.png"
)

# =============================================================================
# 3) RNA-seq Volcano — All TCGA vs P0
# =============================================================================
rna_tcga <- prepare_rna_data(raw_tcga)

make_volcano(
  df            = rna_tcga,
  fc_col        = "RNAseq_Log2FC",
  pval_col      = "RNAseq_qval",
  title         = "Differential Gene Expression of (All TCGA) vs Normal Bladder",
  colors        = VOLCANO_COLORS_RNA,
  top_n         = TOP_N_GENES,
  point_size    = 1.5,
  point_alpha   = 0.7,
  x_limits      = NULL,
  y_limits      = NULL,
  right_label   = "TCGA",
  left_label    = "PO",
  right_color   = "#b2182b",
  left_color    = "#2166ac",
  label_x_pos   = c(-7, 8),
  label_y_pos   = "data_min",
  base_size     = 14,
  output_file   = "volcano_RNA_allTCGA_vs_P0.png"
)

# =============================================================================
# 4) RNA-seq Volcano — LumP vs P0 (with jitter)
# =============================================================================
rna_lump <- prepare_rna_data(raw_lump)

make_volcano(
  df            = rna_lump,
  fc_col        = "RNAseq_Log2FC",
  pval_col      = "RNAseq_qval",
  title         = "Differential Gene Expression of (LumP) vs Normal Bladder",
  colors        = VOLCANO_COLORS_RNA,
  top_n         = TOP_N_GENES,
  point_size    = 1.5,
  point_alpha   = 0.7,
  x_limits      = NULL,
  y_limits      = NULL,
  right_label   = "LumP",
  left_label    = "PO",
  right_color   = "#b2182b",
  left_color    = "#2166ac",
  label_x_pos   = c(-7, 8),
  label_y_pos   = "data_min",
  jitter        = TRUE,
  jitter_seed   = 123,
  base_size     = 14,
  output_file   = "volcano_RNA_LumP_vs_P0.png"
)

message("\n[01_volcanos] All four volcano plots complete.")
