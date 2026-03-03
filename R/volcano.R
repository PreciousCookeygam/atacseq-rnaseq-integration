# =============================================================================
# volcano.R — Volcano plot construction
# =============================================================================
#
# Replaces the four near-identical volcano blocks from the original script with
# a single parameterised function.  Supports both ATAC-seq and RNA-seq styles
# through arguments rather than code duplication.
#
# Usage:
#   source(here::here("R", "volcano.R"))
#   p <- make_volcano(df, fc_col = "log2FoldChange", pval_col = "padj",
#                     title = "ATAC: All TCGA vs P0")
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

#' Classify genes as Up / Down / Not Sig
#'
#' Adds `regulation` and `neg_log10_pval` columns to the data frame.
#'
#' @param df           Data frame with at least the columns named by
#'                     \code{fc_col} and \code{pval_col}.
#' @param fc_col       Character. Column name for log2 fold-change.
#' @param pval_col     Character. Column name for adjusted p-value / q-value.
#' @param fc_threshold Numeric. Absolute fold-change threshold (default from
#'                     \code{FC_THRESHOLD}).
#' @param padj_threshold Numeric. Significance threshold (default from
#'                       \code{PADJ_THRESHOLD}).
#' @return The input data frame with two new columns appended.
classify_regulation <- function(df,
                                fc_col = "log2FoldChange",
                                pval_col = "padj",
                                fc_threshold = FC_THRESHOLD,
                                padj_threshold = PADJ_THRESHOLD) {
  df %>%
    dplyr::mutate(
      regulation = dplyr::case_when(
        .data[[fc_col]] >=  fc_threshold & .data[[pval_col]] < padj_threshold ~ "Up",
        .data[[fc_col]] <= -fc_threshold & .data[[pval_col]] < padj_threshold ~ "Down",
        TRUE ~ "Not Sig"
      ),
      neg_log10_pval = -log10(.data[[pval_col]])
    )
}


#' Build a publication-quality volcano plot
#'
#' A single function that replaces the four copy-pasted volcano blocks from the
#' original analysis.
#'
#' @param df             Data frame (one row per gene, deduplicated).
#' @param fc_col         Character. Column with log2 fold-change values.
#' @param pval_col       Character. Column with adjusted p-value / q-value.
#' @param gene_col       Character. Column with gene symbols (default
#'                       `"Gene_Name"`).
#' @param title          Character. Plot title.
#' @param fc_threshold   Numeric. |FC| cut-off.
#' @param padj_threshold Numeric. Significance cut-off.
#' @param colors         Named character vector of three colours
#'                       (`Up`, `Down`, `Not Sig`).
#' @param top_n          Integer. Number of top Up and Down genes to label.
#' @param point_size     Numeric. Size of geom_point.
#' @param point_alpha    Numeric. Alpha transparency.
#' @param label_size     Numeric. Font size for gene labels.
#' @param x_limits       Numeric vector of length 2, or NULL for auto.
#' @param y_limits       Numeric vector of length 2, or NULL for auto.
#' @param left_label     Character or NULL. Text annotation on the left side.
#' @param right_label    Character or NULL. Text annotation on the right side.
#' @param left_color     Colour for \code{left_label}.
#' @param right_color    Colour for \code{right_label}.
#' @param label_x_pos    Numeric vector of length 2 (left x, right x).
#' @param label_y_pos    Numeric or `"data_min"`. Y position for side labels.
#' @param jitter         Logical. If TRUE, add reproducible jitter to the FC
#'                       and -log10(p) axes (used for the LumP RNA volcano).
#' @param jitter_seed    Integer. Seed for reproducibility when jittering.
#' @param base_size      Numeric. Base font size for \code{theme_publication}.
#' @param output_file    Character or NULL. If provided, the plot is saved to
#'                       \code{file.path(outdir_path, output_file)}.
#' @param outdir_path    Character. Output directory (default \code{OUTPUT_DIR}).
#' @param width          Numeric. Saved plot width in inches (default 10).
#' @param height         Numeric. Saved plot height in inches (default 7).
#' @param dpi            Numeric. Resolution for PNG output (default 300).
#' @return A \code{ggplot} object (also printed and optionally saved).
make_volcano <- function(df,
                         fc_col         = "log2FoldChange",
                         pval_col       = "padj",
                         gene_col       = "Gene_Name",
                         title          = "Volcano Plot",
                         fc_threshold   = FC_THRESHOLD,
                         padj_threshold = PADJ_THRESHOLD,
                         colors         = VOLCANO_COLORS_ATAC,
                         top_n          = TOP_N_GENES,
                         point_size     = 1.8,
                         point_alpha    = 0.8,
                         label_size     = 4,
                         x_limits       = c(-10, 10),
                         y_limits       = c(0, 10),
                         left_label     = NULL,
                         right_label    = NULL,
                         left_color     = "blue",
                         right_color    = "red",
                         label_x_pos    = c(-8, 8),
                         label_y_pos    = 0.5,
                         jitter         = FALSE,
                         jitter_seed    = 123,
                         base_size      = 13,
                         output_file    = NULL,
                         outdir_path    = OUTPUT_DIR,
                         width          = 10,
                         height         = 7,
                         dpi            = 300) {

  # --- classify regulation ---------------------------------------------------
  df <- classify_regulation(df, fc_col, pval_col, fc_threshold, padj_threshold)

  # --- optional jitter --------------------------------------------------------
  if (jitter) {
    set.seed(jitter_seed)
    df <- df %>%
      dplyr::mutate(
        .jitter_fc = .data[[fc_col]] + stats::runif(dplyr::n(), -0.2, 0.2),
        .jitter_pv = neg_log10_pval  + stats::runif(dplyr::n(), -0.4, 0.4)
      )
    x_aes <- ".jitter_fc"
    y_aes <- ".jitter_pv"
  } else {
    x_aes <- fc_col
    y_aes <- "neg_log10_pval"
  }

  # --- select top genes for labelling -----------------------------------------
  top_up <- df %>%
    dplyr::filter(regulation == "Up") %>%
    dplyr::arrange(.data[[pval_col]]) %>%
    dplyr::slice_head(n = top_n)

  top_down <- df %>%
    dplyr::filter(regulation == "Down") %>%
    dplyr::arrange(.data[[pval_col]]) %>%
    dplyr::slice_head(n = top_n)

  annotated_genes <- dplyr::bind_rows(top_up, top_down)

  n_up   <- sum(df$regulation == "Up")
  n_down <- sum(df$regulation == "Down")
  message(sprintf("[volcano] %s  |  Up: %d  Down: %d", title, n_up, n_down))

  # --- build the plot ---------------------------------------------------------
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_aes]],
                                         y = .data[[y_aes]],
                                         color = regulation)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggrepel::geom_text_repel(
      data         = annotated_genes,
      ggplot2::aes(label = .data[[gene_col]]),
      size         = label_size,
      box.padding  = 0.4,
      point.padding = 0.3,
      segment.color = "grey50",
      max.overlaps  = 100,
      show.legend   = FALSE
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_vline(
      xintercept = c(-fc_threshold, fc_threshold),
      linetype = "dashed", color = "grey50"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_threshold),
      linetype = "dashed", color = "grey50"
    ) +
    ggplot2::labs(
      title = title,
      x     = expression(log[2]("Fold Change")),
      y     = expression(-log[10]("adj. P-value")),
      color = "Regulation"
    ) +
    theme_publication(base_size = base_size)

  # --- axis limits (if specified) ---------------------------------------------
  if (!is.null(x_limits)) p <- p + ggplot2::xlim(x_limits)
  if (!is.null(y_limits)) p <- p + ggplot2::ylim(y_limits)

  # --- side annotations -------------------------------------------------------
  if (!is.null(right_label)) {
    y_pos <- if (identical(label_y_pos, "data_min")) {
      min(df$neg_log10_pval, na.rm = TRUE)
    } else {
      label_y_pos
    }
    p <- p +
      ggplot2::annotate("text", x = label_x_pos[2], y = y_pos,
                        label = right_label, size = 3,
                        fontface = "bold", color = right_color)
  }
  if (!is.null(left_label)) {
    y_pos <- if (identical(label_y_pos, "data_min")) {
      min(df$neg_log10_pval, na.rm = TRUE)
    } else {
      label_y_pos
    }
    p <- p +
      ggplot2::annotate("text", x = label_x_pos[1], y = y_pos,
                        label = left_label, size = 3,
                        fontface = "bold", color = left_color)
  }

  # --- save to file (if requested) -------------------------------------------
  if (!is.null(output_file)) {
    out_path <- file.path(outdir_path, output_file)
    ggplot2::ggsave(out_path, p, width = width, height = height, dpi = dpi)
    message(sprintf("[volcano] Saved: %s", out_path))
  }

  print(p)
  invisible(p)
}
