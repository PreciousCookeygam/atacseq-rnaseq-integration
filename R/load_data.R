# =============================================================================
# load_data.R — Data loading and preparation
# =============================================================================
#
# Provides a single-load caching layer so each large TSV is read at most once
# per session, and helper functions that extract ATAC-only, RNA-only, or
# combined columns with deduplication by gene name.
#
# Usage:
#   source(here::here("R", "load_data.R"))
#   df_tcga <- load_annotated_results(FILE_TCGA)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ---------------------------------------------------------------------------
# Internal cache environment (avoids polluting the global namespace)
# ---------------------------------------------------------------------------

.data_cache <- new.env(parent = emptyenv())

#' Load an annotated results TSV with caching
#'
#' Reads the file once and stores it in an internal environment. Subsequent
#' calls with the same path return the cached copy without touching disk.
#'
#' @param path Character. Absolute or project-relative path to a TSV file.
#' @return A tibble with all columns from the file.
load_annotated_results <- function(path) {
  key <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (exists(key, envir = .data_cache)) {
    return(get(key, envir = .data_cache))
  }
  message(sprintf("[load_data] Reading: %s", basename(path)))
  df <- readr::read_tsv(path, show_col_types = FALSE)
  assign(key, df, envir = .data_cache)
  df
}

#' Clear the data cache (useful during interactive development)
#'
#' @return NULL (invisible).
clear_data_cache <- function() {
  rm(list = ls(envir = .data_cache), envir = .data_cache)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Column-selection helpers (with deduplication by gene)
# ---------------------------------------------------------------------------

#' Prepare ATAC-seq columns for analysis
#'
#' Selects `Gene_Name`, `log2FoldChange`, `padj`; removes rows with any NA
#' in those columns; and keeps the first occurrence per gene.
#'
#' @param df A tibble loaded by \code{\link{load_annotated_results}}.
#' @return A deduplicated tibble with three columns.
prepare_atac_data <- function(df) {
  df %>%
    dplyr::select(Gene_Name, log2FoldChange, padj) %>%
    dplyr::filter(!is.na(Gene_Name), !is.na(log2FoldChange), !is.na(padj)) %>%
    dplyr::distinct(Gene_Name, .keep_all = TRUE)
}

#' Prepare RNA-seq columns for analysis
#'
#' Selects `Gene_Name`, `RNAseq_Log2FC`, `RNAseq_qval`; removes incomplete
#' rows; and keeps the first occurrence per gene.
#'
#' @param df A tibble loaded by \code{\link{load_annotated_results}}.
#' @return A deduplicated tibble with three columns.
prepare_rna_data <- function(df) {
  df %>%
    dplyr::select(Gene_Name, RNAseq_Log2FC, RNAseq_qval) %>%
    dplyr::filter(!is.na(Gene_Name), !is.na(RNAseq_Log2FC), !is.na(RNAseq_qval)) %>%
    dplyr::distinct(Gene_Name, .keep_all = TRUE)
}

#' Prepare combined ATAC + RNA columns for pi-score / concordance analysis
#'
#' Selects both ATAC (`log2FoldChange`, `padj`) and RNA (`RNAseq_Log2FC`,
#' `RNAseq_qval`) columns together with `Gene_Name`, removes rows missing
#' the ATAC values (RNA NAs are kept for flexibility), and deduplicates.
#'
#' @param df A tibble loaded by \code{\link{load_annotated_results}}.
#' @return A deduplicated tibble with five columns.
prepare_combined_data <- function(df) {
  df %>%
    dplyr::select(Gene_Name, log2FoldChange, padj, RNAseq_Log2FC, RNAseq_qval) %>%
    dplyr::filter(!is.na(Gene_Name), !is.na(log2FoldChange), !is.na(padj)) %>%
    dplyr::distinct(Gene_Name, .keep_all = TRUE)
}
