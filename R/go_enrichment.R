# =============================================================================
# go_enrichment.R — Gene Ontology enrichment and standardised dotplots
# =============================================================================
#
# Cleans up and de-duplicates the GO enrichment logic from the original
# venn_diagram.R.  The enrichment function is called exactly once per gene
# set (not twice as in the original), and the standardised dotplot builder
# is extracted into its own reusable function.
#
# Usage:
#   source(here::here("R", "go_enrichment.R"))
#   ego <- go_enrich_and_plot(gene_symbols, "ATAC_core_UP",
#                             universe_symbols = bg, outdir_path = OUTPUT_DIR)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(scales)
})

#' Run GO enrichment and export dotplot + barplot PNGs
#'
#' Maps gene symbols to Entrez IDs, runs \code{clusterProfiler::enrichGO},
#' and saves a dotplot and barplot as PNG files in \code{outdir_path}.
#'
#' @param gene_symbols     Character vector of HGNC gene symbols.
#' @param list_label        Character. Label used in file names and plot titles.
#' @param universe_symbols  Character vector of background gene symbols, or NULL.
#' @param ont               Character. GO ontology: `"BP"`, `"MF"`, or `"CC"`.
#' @param p_cutoff          Numeric. P-value cut-off for enrichGO.
#' @param q_cutoff          Numeric. Q-value cut-off for enrichGO.
#' @param outdir_path       Character. Directory to write output files.
#' @return An \code{enrichResult} object (invisible), or NULL on failure.
go_enrich_and_plot <- function(gene_symbols,
                               list_label,
                               universe_symbols = NULL,
                               ont              = "BP",
                               p_cutoff         = 0.05,
                               q_cutoff         = 0.10,
                               outdir_path      = OUTPUT_DIR) {

  # --- Map SYMBOL → ENTREZ (dedupe by SYMBOL) --------------------------------
  mapped <- suppressMessages(
    clusterProfiler::bitr(gene_symbols,
                          fromType = "SYMBOL",
                          toType   = "ENTREZID",
                          OrgDb    = org.Hs.eg.db)
  )
  if (is.null(mapped) || nrow(mapped) == 0) {
    message(sprintf("[GO] %s — no symbols mapped to Entrez. Skipping.", list_label))
    return(invisible(NULL))
  }
  mapped <- dplyr::distinct(mapped, SYMBOL, .keep_all = TRUE)

  # --- Optional background universe ------------------------------------------
  bg_entrez <- NULL
  if (!is.null(universe_symbols)) {
    bg_map <- suppressMessages(
      clusterProfiler::bitr(universe_symbols,
                            fromType = "SYMBOL",
                            toType   = "ENTREZID",
                            OrgDb    = org.Hs.eg.db)
    )
    if (!is.null(bg_map) && nrow(bg_map) > 0) {
      bg_map    <- dplyr::distinct(bg_map, SYMBOL, .keep_all = TRUE)
      bg_entrez <- unique(bg_map$ENTREZID)
    }
  }

  # --- Enrichment -------------------------------------------------------------
  ego <- tryCatch({
    clusterProfiler::enrichGO(
      gene          = unique(mapped$ENTREZID),
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      universe      = bg_entrez,
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff,
      readable      = TRUE
    )
  }, error = function(e) {
    message(sprintf("[GO] %s — enrichGO error: %s", list_label, e$message))
    return(NULL)
  })

  if (is.null(ego)) return(invisible(NULL))

  tbl <- as.data.frame(ego)
  if (nrow(tbl) == 0) {
    message(sprintf("[GO] %s — no enriched GO terms at current cut-offs.", list_label))
    return(invisible(ego))
  }

  # --- Save results table -----------------------------------------------------
  readr::write_tsv(
    tbl,
    file.path(outdir_path, sprintf("%s_GO_%s_results.tsv", list_label, ont))
  )

  # --- Dotplot ----------------------------------------------------------------
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  p_dot <- enrichplot::dotplot(ego, showCategory = 20,
                               title = sprintf("%s — GO:%s", list_label, ont))
  ggplot2::ggsave(
    file.path(outdir_path, sprintf("%s_GO_%s_dotplot_%s.png", list_label, ont, stamp)),
    p_dot, width = 8, height = 6, dpi = 300
  )

  # --- Barplot ----------------------------------------------------------------
  p_bar <- barplot(ego, showCategory = 20,
                   title = sprintf("%s — GO:%s", list_label, ont))
  ggplot2::ggsave(
    file.path(outdir_path, sprintf("%s_GO_%s_barplot_%s.png", list_label, ont, stamp)),
    p_bar, width = 8, height = 6, dpi = 300
  )

  message(sprintf("[GO] %s — %d enriched terms. Saved to %s/", list_label,
                  nrow(tbl), outdir_path))
  invisible(ego)
}


# ---------------------------------------------------------------------------
# Standardised dotplot helpers
# ---------------------------------------------------------------------------

#' Compute shared scales across multiple enrichResult objects
#'
#' Determines the global ranges for GeneRatio, Count, and p.adjust so that
#' multiple dotplots can share identical axes and colour bars.
#'
#' @param ego_list A list of \code{enrichResult} objects (NULLs are tolerated).
#' @return A named list: `max_ratio`, `max_count`, `padj_range`, `x_breaks`.
compute_shared_scales <- function(ego_list) {

  tbls <- lapply(ego_list, function(e) {
    if (is.null(e)) NULL else as.data.frame(e)
  })

  max_ratio <- max(vapply(tbls, function(t) {
    if (is.null(t)) 0 else max(ratio_to_numeric(t$GeneRatio), na.rm = TRUE)
  }, numeric(1)), na.rm = TRUE)


  max_count <- max(vapply(tbls, function(t) {
    if (is.null(t)) 0 else max(t$Count, na.rm = TRUE)
  }, numeric(1)), na.rm = TRUE)

  padj_rng <- range(unlist(lapply(tbls, function(t) {
    if (is.null(t)) NA_real_ else t$p.adjust
  })), na.rm = TRUE)

  x_breaks <- seq(0, max_ratio, length.out = 4)

  list(
    max_ratio  = max_ratio,
    max_count  = max_count,
    padj_range = padj_rng,
    x_breaks   = x_breaks
  )
}


#' Generate -log10(FDR) colour-bar breaks
#'
#' @param padj_range Numeric vector of length 2 (min, max p.adjust).
#' @param n          Number of breaks.
#' @return Numeric vector of p-values to use as colour breaks.
neglog10_breaks <- function(padj_range, n = 5) {
  lo <- floor(-log10(max(padj_range)))
  hi <- ceiling(-log10(min(padj_range)))
  b  <- seq(lo, hi, length.out = n)
  10^(-b)
}


#' Build a standardised dotplot with shared axes
#'
#' @param ego          An \code{enrichResult} object.
#' @param title        Character. Plot title.
#' @param shared       List returned by \code{compute_shared_scales}.
#' @param show         Integer. Number of categories to show.
#' @param wrap         Integer. Label wrap width.
#' @param outdir_path  Character. Output directory for the saved PNG.
#' @param filename     Character. Output filename (without path).
#' @return A \code{ggplot} object (invisible, also printed).
make_standardised_dotplot <- function(ego,
                                      title,
                                      shared,
                                      show        = 12,
                                      wrap        = 45,
                                      outdir_path = OUTPUT_DIR,
                                      filename    = NULL) {

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message(sprintf("[dotplot] %s — no data, skipping.", title))
    return(invisible(NULL))
  }

  col_breaks <- neglog10_breaks(shared$padj_range, n = 5)
  col_labels <- function(x) sprintf("%.1f", -log10(x))

  p <- enrichplot::dotplot(ego, showCategory = show, label_format = wrap) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(6, 12, 6, 32, "pt")
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, shared$max_ratio),
      breaks = shared$x_breaks,
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::scale_size_continuous(
      limits = c(0, shared$max_count),
      range  = c(3, 10)
    ) +
    ggplot2::scale_colour_gradient(
      limits = shared$padj_range,
      breaks = col_breaks,
      labels = col_labels,
      trans  = "reverse",
      low    = "red",
      high   = "blue",
      oob    = scales::squish
    ) +
    ggplot2::guides(
      size   = ggplot2::guide_legend(title = "Gene count"),
      colour = ggplot2::guide_colourbar(title = expression(-log[10] * "(FDR)"))
    )

  print(p)

  if (!is.null(filename)) {
    ggplot2::ggsave(
      file.path(outdir_path, filename), p,
      width = 12, height = 8, dpi = 300
    )
  }

  invisible(p)
}
