# =============================================================================
# concordance.R — Cross-modality concordance analysis
# =============================================================================
#
# Replaces the two near-identical concordance blocks (TCGA / LumP) that used
# inconsistent variable names (N_A/X_A/Y_A vs M_A/U_A/V_A).
#
# Usage:
#   source(here::here("R", "concordance.R"))
#   stats <- compute_concordance(df, label = "All TCGA")
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

#' Compute RNA-seq / ATAC-seq concordance statistics
#'
#' A gene is "concordant" when it is significant in **both** modalities
#' (|log2FC| ≥ threshold & adjusted p < threshold) **and** the fold-change
#' direction is the same.
#'
#' Two denominators are reported:
#' \describe{
#'   \item{Denominator A}{All genes with complete RNA + ATAC values.}
#'   \item{Denominator B}{Only genes significant in at least one modality.}
#' }
#'
#' @param df    Data frame produced by \code{\link{prepare_combined_data}}.
#'              Must contain `RNAseq_Log2FC`, `RNAseq_qval`,
#'              `log2FoldChange`, `padj`.
#' @param label Character. Label for the comparison (e.g. `"All TCGA"`,
#'              `"LumP"`). Used in console output.
#' @param fc_threshold   Numeric. |FC| cut-off.
#' @param padj_threshold Numeric. Adjusted p-value cut-off.
#' @return A named list with elements:
#'   \describe{
#'     \item{merged}{The annotated data frame with significance columns.}
#'     \item{concordant_genes}{Character vector of concordant gene symbols.}
#'     \item{n_total}{Denominator A count.}
#'     \item{n_concordant}{Numerator for denominator A.}
#'     \item{pct_concordant}{Percentage (denominator A).}
#'     \item{n_sig_any}{Denominator B count.}
#'     \item{n_concordant_sig}{Numerator for denominator B.}
#'     \item{pct_concordant_sig}{Percentage (denominator B).}
#'   }
compute_concordance <- function(df,
                                label = "Comparison",
                                fc_threshold = FC_THRESHOLD,
                                padj_threshold = PADJ_THRESHOLD) {

  merged <- df %>%
    dplyr::filter(!is.na(RNAseq_Log2FC), !is.na(RNAseq_qval)) %>%
    dplyr::mutate(
      RNAseq_qval  = as.numeric(RNAseq_qval),
      padj          = as.numeric(padj),
      pi_score_RNA  = sign(RNAseq_Log2FC) * -log10(RNAseq_qval),
      pi_score_ATAC = sign(log2FoldChange) * -log10(padj),
      RNAseq_Significance  = (abs(RNAseq_Log2FC)  >= fc_threshold) &
                              (RNAseq_qval < padj_threshold),
      ATACseq_Significance = (abs(log2FoldChange) >= fc_threshold) &
                              (padj < padj_threshold),
      same_dir   = sign(RNAseq_Log2FC) == sign(log2FoldChange),
      concordant = RNAseq_Significance & ATACseq_Significance & same_dir
    )

  # Denominator A: all genes with complete data
  n_total       <- nrow(merged)
  n_concordant  <- sum(merged$concordant, na.rm = TRUE)
  pct_total     <- round(100 * n_concordant / n_total, 1)

  # Denominator B: genes significant in ≥ 1 modality
  sig_any           <- merged %>%
    dplyr::filter(RNAseq_Significance | ATACseq_Significance)
  n_sig_any         <- nrow(sig_any)
  n_concordant_sig  <- sum(sig_any$concordant, na.rm = TRUE)
  pct_sig           <- round(100 * n_concordant_sig / n_sig_any, 1)

  concordant_genes <- merged %>%
    dplyr::filter(concordant) %>%
    dplyr::pull(Gene_Name) %>%
    unique()

  # Console summary
  cat(sprintf(
    "\n=== Concordance: %s ===\n  Denom A (complete data): %d / %d (%.1f%%)\n  Denom B (sig in >= 1):   %d / %d (%.1f%%)\n  Concordant genes: %d\n",
    label, n_concordant, n_total, pct_total,
    n_concordant_sig, n_sig_any, pct_sig,
    length(concordant_genes)
  ))

  invisible(list(
    merged             = merged,
    concordant_genes   = concordant_genes,
    n_total            = n_total,
    n_concordant       = n_concordant,
    pct_concordant     = pct_total,
    n_sig_any          = n_sig_any,
    n_concordant_sig   = n_concordant_sig,
    pct_concordant_sig = pct_sig
  ))
}
