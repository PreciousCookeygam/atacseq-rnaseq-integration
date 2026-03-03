# =============================================================================
# config.R — Centralised project configuration
# =============================================================================
#
# All analysis thresholds, colour palettes, and shared plot settings live here.
# Import this file at the top of every analysis script so constants are never
# hard-coded inline.
#
# Usage:
#   source(here::here("R", "config.R"))
# =============================================================================

library(here)

# ---------------------------------------------------------------------------
# Data file paths (relative to project root via `here`)
# ---------------------------------------------------------------------------

DATA_DIR <- here::here("data")

FILE_TCGA <- here::here("data", "ATACseq_annotated_results_TCGA-BLCA_vs_P0-Bladder.tsv")
FILE_LUMP <- here::here("data", "ATACseq_annotated_results_TCGA-BLCA-LumP_vs_P0-Bladder.tsv")
FILE_LUMP_VS_BASQ <- here::here("data", "ATACseq_annotated_results_TCGA-BLCA-LumP_vs_TCGA-BLCA-BaSq.tsv")
FILE_MERGED <- here::here("data", "Merged_ATAC_RNA_results.tsv")
FILE_CONCORDANT_GENES <- here::here("data", "concordant_core_152_genes.tsv")
FILE_RNA_RANK <- here::here("data", "rna_seq_rank.rnk")

# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------

OUTPUT_DIR <- here::here("output")

# ---------------------------------------------------------------------------
# Statistical thresholds
# ---------------------------------------------------------------------------

FC_THRESHOLD   <- 1      # |log2 fold-change| cut-off
PADJ_THRESHOLD <- 0.05   # Adjusted p-value / q-value cut-off

# ---------------------------------------------------------------------------
# Gene labelling
# ---------------------------------------------------------------------------

TOP_N_GENES     <- 10   # Top genes to label on volcano plots
TOP_N_PI_GENES  <- 20   # Top genes to label on pi-score plots

# ---------------------------------------------------------------------------
# Colour palettes
# ---------------------------------------------------------------------------

#' ATAC-seq volcano colours
VOLCANO_COLORS_ATAC <- c(
  "Up"      = "red",
  "Down"    = "blue",
  "Not Sig" = "grey"
)

#' RNA-seq volcano colours
VOLCANO_COLORS_RNA <- c(
  "Up"      = "#b2182b",
  "Down"    = "#2166ac",
  "Not Sig" = "lightgrey"
)

#' Pi-score significance groups
PI_SCORE_COLORS <- c(
  "Both" = "black",
  "RNA"  = "red",
  "ATAC" = "blue",
  "None" = "lightgrey"
)

#' Venn diagram fill
VENN_FILL <- c("#87CEEB", "#FF00FF")
VENN_ALPHA <- c(0.6, 0.6)

# ---------------------------------------------------------------------------
# Publication-quality base theme
# ---------------------------------------------------------------------------

#' Consistent ggplot2 theme for all plots
#'
#' @param base_size Numeric. Base font size (default 13).
#' @return A ggplot2 theme object.
theme_publication <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid     = ggplot2::element_blank(),
      panel.border   = ggplot2::element_blank(),
      axis.line.x    = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y    = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks     = ggplot2::element_line(color = "black"),
      axis.text      = ggplot2::element_text(size = base_size - 1),
      axis.title     = ggplot2::element_text(size = base_size + 1, face = "bold"),
      plot.title     = ggplot2::element_text(size = base_size + 1, hjust = 0.5,
                                             face = "bold"),
      legend.title   = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.margin    = ggplot2::margin(10, 40, 10, 10)
    )
}
