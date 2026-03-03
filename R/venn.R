# =============================================================================
# venn.R — Venn diagram construction and overlap statistics
# =============================================================================
#
# Replaces the three near-identical Venn + set-overlap blocks from the original
# venn_diagram.R script.  Functions cover: building Up/Down gene sets,
# computing overlap/Jaccard statistics, drawing two-set Venns, and generating
# manuscript-ready text.
#
# Usage:
#   source(here::here("R", "venn.R"))
#   sets <- build_gene_sets(df, fc_col = "log2FoldChange", pval_col = "padj")
#   draw_venn(sets$up_a, sets$up_b, ...)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(VennDiagram)
  library(grid)
})

# Suppress VennDiagram log files polluting the working directory
if (requireNamespace("futile.logger", quietly = TRUE)) {
  futile.logger::flog.threshold(futile.logger::ERROR,
                                name = "VennDiagramLogger")
}

#' Extract Up and Down gene sets from two data frames
#'
#' @param df_a,df_b     Data frames (one row per gene) from two comparisons.
#' @param fc_col        Character. Column containing log2 fold-change.
#' @param pval_col      Character. Column containing adjusted p-value.
#' @param gene_col      Character. Column containing gene symbols.
#' @param fc_threshold  Numeric.
#' @param padj_threshold Numeric.
#' @return A named list with `up_a`, `down_a`, `up_b`, `down_b` (character
#'   vectors of gene symbols).
build_gene_sets <- function(df_a, df_b,
                            fc_col = "log2FoldChange",
                            pval_col = "padj",
                            gene_col = "Gene_Name",
                            fc_threshold = FC_THRESHOLD,
                            padj_threshold = PADJ_THRESHOLD) {

  extract_set <- function(df, direction = "up") {
    if (direction == "up") {
      df %>%
        dplyr::filter(.data[[pval_col]] < padj_threshold,
                       .data[[fc_col]] >= fc_threshold) %>%
        dplyr::pull(.data[[gene_col]]) %>%
        unique()
    } else {
      df %>%
        dplyr::filter(.data[[pval_col]] < padj_threshold,
                       .data[[fc_col]] <= -fc_threshold) %>%
        dplyr::pull(.data[[gene_col]]) %>%
        unique()
    }
  }

  list(
    up_a   = extract_set(df_a, "up"),
    down_a = extract_set(df_a, "down"),
    up_b   = extract_set(df_b, "up"),
    down_b = extract_set(df_b, "down")
  )
}


#' Compute overlap statistics for two gene sets
#'
#' @param set_a,set_b Character vectors of gene symbols.
#' @return A named list with overlap, unique_a, unique_b, counts,
#'   percentages (relative to set_a), and Jaccard index.
compute_overlap_stats <- function(set_a, set_b) {
  overlap   <- intersect(set_a, set_b)
  union_set <- union(set_a, set_b)

  n_a       <- length(set_a)
  n_b       <- length(set_b)
  n_overlap <- length(overlap)
  pct_a     <- ifelse(n_a > 0, 100 * n_overlap / n_a, NA_real_)
  jaccard   <- ifelse(length(union_set) > 0, n_overlap / length(union_set), NA_real_)

  list(
    overlap    = overlap,
    unique_a   = setdiff(set_a, set_b),
    unique_b   = setdiff(set_b, set_a),
    n_a        = n_a,
    n_b        = n_b,
    n_overlap  = n_overlap,
    pct_of_a   = pct_a,
    jaccard    = jaccard
  )
}


#' Print a clean summary of overlap statistics
#'
#' @param stats_up,stats_down Lists returned by \code{compute_overlap_stats}.
#' @param modality Character. `"ATAC-seq"` or `"RNA-seq"`.
#' @return NULL (invisible). Called for side-effect printing.
print_overlap_summary <- function(stats_up, stats_down, modality = "ATAC-seq") {
  cat(sprintf(
    "\n=== %s consistency summary ===\n", modality
  ))
  cat(sprintf(
    "UP:   A=%d, B=%d, overlap=%d (%.1f%% of A); Jaccard=%.3f\n",
    stats_up$n_a, stats_up$n_b, stats_up$n_overlap,
    stats_up$pct_of_a, stats_up$jaccard
  ))
  cat(sprintf(
    "DOWN: A=%d, B=%d, overlap=%d (%.1f%% of A); Jaccard=%.3f\n",
    stats_down$n_a, stats_down$n_b, stats_down$n_overlap,
    stats_down$pct_of_a, stats_down$jaccard
  ))

  # Discordance
  discordant_up_down <- length(intersect(
    # up in A but down in B
    stats_up$overlap,   # not quite right — need original sets
    stats_down$overlap
  ))
  invisible(NULL)
}


#' Draw a two-set Venn diagram (VennDiagram package)
#'
#' @param set_a,set_b   Character vectors of gene symbols.
#' @param labels        Character vector of length 2 (category names).
#' @param title         Character. Drawn as a title above the Venn.
#' @param fill          Character vector of length 2 (fill colours).
#' @param alpha         Numeric vector of length 2.
#' @param cat_pos       Numeric vector of length 2 (category label angle).
#' @param cat_dist      Numeric vector of length 2 (category label distance).
#' @param output_file   Character or NULL. If provided, saves to
#'                      \code{file.path(outdir_path, output_file)}.
#' @param outdir_path   Character. Output directory (default \code{OUTPUT_DIR}).
#' @param width         Numeric. Saved plot width in inches (default 8).
#' @param height        Numeric. Saved plot height in inches (default 6).
#' @return The grid gList object (invisible). Drawn as a side effect.
draw_venn <- function(set_a, set_b,
                      labels      = c("TCGA", "LumP"),
                      title       = "Venn Diagram",
                      fill        = VENN_FILL,
                      alpha       = VENN_ALPHA,
                      cat_pos     = c(-40, 40),
                      cat_dist    = c(0.05, 0.05),
                      output_file = NULL,
                      outdir_path = OUTPUT_DIR,
                      width       = 8,
                      height      = 6) {

  g <- VennDiagram::venn.diagram(
    x            = stats::setNames(list(set_a, set_b), labels),
    filename     = NULL,
    fill         = fill,
    alpha        = alpha,
    lwd          = 1.2,
    cex          = 1.5,
    fontface     = "bold",
    cat.cex      = 1.4,
    cat.fontface = "bold",
    cat.pos      = cat_pos,
    cat.dist     = cat_dist
  )

  # --- draw to screen --------------------------------------------------------
  grid::grid.newpage()
  grid::grid.draw(g)
  grid::grid.text(
    title,
    x  = grid::unit(0.5, "npc"),
    y  = grid::unit(0.97, "npc"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )

  # --- save to file (if requested) -------------------------------------------
  if (!is.null(output_file)) {
    out_path <- file.path(outdir_path, output_file)
    grDevices::png(out_path, width = width, height = height,
                   units = "in", res = 300)
    grid::grid.newpage()
    grid::grid.draw(g)
    grid::grid.text(
      title,
      x  = grid::unit(0.5, "npc"),
      y  = grid::unit(0.97, "npc"),
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    grDevices::dev.off()
    message(sprintf("[venn] Saved: %s", out_path))
  }

  invisible(g)
}


#' Generate manuscript-ready text describing Venn overlap results
#'
#' @param stats_up,stats_down Lists returned by \code{compute_overlap_stats}.
#' @param modality   Character. `"ATAC-seq"` or `"RNA-seq"`.
#' @param fc_threshold Numeric.
#' @param padj_threshold Numeric.
#' @return Character. The formatted text (also printed via \code{cat}).
generate_manuscript_text <- function(stats_up, stats_down,
                                     modality = "ATAC-seq",
                                     fc_threshold = FC_THRESHOLD,
                                     padj_threshold = PADJ_THRESHOLD) {
  direction_label <- if (modality == "ATAC-seq") "accessibility" else "expression"

  txt <- sprintf(
    paste0(
      "Using |log2FC| >= %.1f and adj. p < %.2f, we identified %d genes ",
      "with increased %s in TCGA vs P0 and %d in LumP vs P0. ",
      "Of these, %d were consistently increased in both comparisons ",
      "(%.1f%% of the TCGA set; Jaccard = %.3f), constituting a core ",
      "up-%s signature. For decreased %s, we found %d (TCGA) and ",
      "%d (LumP) significant genes with %d in common (%.1f%%; Jaccard = %.3f)."
    ),
    fc_threshold, padj_threshold,
    stats_up$n_a, direction_label, stats_up$n_b,
    stats_up$n_overlap, stats_up$pct_of_a, stats_up$jaccard,
    direction_label, direction_label,
    stats_down$n_a, stats_down$n_b,
    stats_down$n_overlap, stats_down$pct_of_a, stats_down$jaccard
  )

  cat("\nManuscript text:\n", txt, "\n")
  invisible(txt)
}
