# =============================================================================
# pi_score.R — π-score computation and concordance scatter plots
# =============================================================================
#
# Replaces the two duplicated π-score blocks (All TCGA and LumP) from the
# original script.  The π-score is defined as:
#   π = sign(FC) × –log₁₀(q-value)
#
# Usage:
#   source(here::here("R", "pi_score.R"))
#   scored <- compute_pi_scores(df)
#   p      <- plot_pi_scores(scored, title = "All TCGA vs P0")
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

#' Compute π-scores and significance categories
#'
#' Expects a data frame with columns: `Gene_Name`, `log2FoldChange`, `padj`,
#' `RNAseq_Log2FC`, `RNAseq_qval` (as returned by
#' \code{\link{prepare_combined_data}}).
#'
#' @param df             Data frame (one row per gene).
#' @param fc_threshold   Numeric. |FC| cut-off for significance flags.
#' @param padj_threshold Numeric. Adjusted p-value cut-off.
#' @return The input data frame augmented with `pi_score_RNA`,
#'   `pi_score_ATAC`, `RNAseq_Significance`, `ATACseq_Significance`, and
#'   `color_group`.
compute_pi_scores <- function(df,
                              fc_threshold = FC_THRESHOLD,
                              padj_threshold = PADJ_THRESHOLD) {
  df %>%
    dplyr::mutate(
      RNAseq_qval  = as.numeric(RNAseq_qval),
      padj          = as.numeric(padj),
      pi_score_RNA  = sign(RNAseq_Log2FC) * -log10(RNAseq_qval),
      pi_score_ATAC = sign(log2FoldChange) * -log10(padj),
      RNAseq_Significance  = (abs(RNAseq_Log2FC) >= fc_threshold) &
                              (RNAseq_qval < padj_threshold),
      ATACseq_Significance = (abs(log2FoldChange) >= fc_threshold) &
                              (padj < padj_threshold),
      color_group = dplyr::case_when(
        RNAseq_Significance &  ATACseq_Significance ~ "Both",
        RNAseq_Significance & !ATACseq_Significance ~ "RNA",
        !RNAseq_Significance & ATACseq_Significance ~ "ATAC",
        TRUE ~ "None"
      )
    )
}


#' Draw a π-score concordance scatter plot
#'
#' @param df       Data frame produced by \code{\link{compute_pi_scores}}.
#' @param title    Character. Plot title.
#' @param top_n       Integer. Number of "Both" genes to label (default
#'                    \code{TOP_N_PI_GENES}).
#' @param xy_limit    Numeric. Symmetric axis range (default 7).
#' @param colors      Named character vector of four colours.
#' @param base_size   Numeric. Theme base font size.
#' @param output_file Character or NULL. If provided, saves to
#'                    \code{file.path(outdir_path, output_file)}.
#' @param outdir_path Character. Output directory (default \code{OUTPUT_DIR}).
#' @param width       Numeric. Saved plot width in inches (default 10).
#' @param height      Numeric. Saved plot height in inches (default 7).
#' @param dpi         Numeric. Resolution (default 300).
#' @return A \code{ggplot} object (also printed and optionally saved).
plot_pi_scores <- function(df,
                           title       = expression("π-Score Concordance: RNA-seq vs ATAC-seq"),
                           top_n       = TOP_N_PI_GENES,
                           xy_limit    = 7,
                           colors      = PI_SCORE_COLORS,
                           base_size   = 13,
                           output_file = NULL,
                           outdir_path = OUTPUT_DIR,
                           width       = 10,
                           height      = 7,
                           dpi         = 300) {

  top_genes <- df %>%
    dplyr::filter(color_group == "Both") %>%
    dplyr::arrange(dplyr::desc(abs(pi_score_RNA) + abs(pi_score_ATAC))) %>%
    dplyr::slice_head(n = top_n)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = pi_score_RNA,
                                         y = pi_score_ATAC,
                                         color = color_group)) +
    ggplot2::geom_point(alpha = 0.7, size = 2.5) +
    ggplot2::scale_color_manual(values = colors, name = "Significance") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggrepel::geom_text_repel(
      data          = top_genes,
      ggplot2::aes(label = Gene_Name),
      size          = 4,
      max.overlaps  = 100,
      box.padding   = 0.3,
      point.padding = 0.3,
      color         = "black"
    ) +
    ggplot2::labs(
      title = title,
      x     = expression("RNA-seq π Score  (sign(FC) · –log"[10]*"q)"),
      y     = expression("ATAC-seq π Score  (sign(FC) · –log"[10]*"q)")
    ) +
    ggplot2::coord_cartesian(xlim = c(-xy_limit, xy_limit),
                             ylim = c(-xy_limit, xy_limit)) +
    theme_publication(base_size = base_size)

  # --- save to file (if requested) -------------------------------------------
  if (!is.null(output_file)) {
    out_path <- file.path(outdir_path, output_file)
    ggplot2::ggsave(out_path, p, width = width, height = height, dpi = dpi)
    message(sprintf("[pi_score] Saved: %s", out_path))
  }

  print(p)
  invisible(p)
}
