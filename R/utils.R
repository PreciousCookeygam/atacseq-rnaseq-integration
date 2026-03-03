# =============================================================================
# utils.R — Shared utility helpers
# =============================================================================
#
# Small, general-purpose helper functions used across multiple modules.
#
# Usage:
#   source(here::here("R", "utils.R"))
# =============================================================================

#' Ensure an output directory exists and is writable
#'
#' Creates the directory (recursively) if it does not yet exist, then performs
#' a write test to verify the path is usable.
#'
#' @param path Character. Path to the output directory.
#' @return The normalised absolute path (character), invisibly.
ensure_output_dir <- function(path) {
  if (file.exists(path) && !dir.exists(path)) {
    stop(
      sprintf("Path exists but is a FILE, not a folder: %s. Delete or rename it.", path),
      call. = FALSE
    )
  }

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  # Quick write test

  tmp <- file.path(path, paste0(".write_test_", Sys.getpid(), ".tmp"))
  ok  <- tryCatch(file.create(tmp), error = function(e) FALSE, warning = function(w) FALSE)
  if (!isTRUE(ok)) {
    stop(
      sprintf("Cannot write to '%s' (cwd: %s). Check permissions.", path, getwd()),
      call. = FALSE
    )
  }
  unlink(tmp, force = TRUE)

  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}


#' Generate a timestamped filename
#'
#' Appends a `YYYYMMDD-HHMMSS` timestamp before the file extension to avoid
#' overwriting previous outputs.
#'
#' @param prefix Character. Base name (e.g. `"ATAC_core_UP_GO_BP_dotplot"`).
#' @param ext    Character. File extension including the dot (default `".png"`).
#' @return A character string like `"ATAC_core_UP_GO_BP_dotplot_20260303-143021.png"`.
timestamped_filename <- function(prefix, ext = ".png") {
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  paste0(prefix, "_", stamp, ext)
}


#' Convert a GeneRatio string (e.g. "12/300") to a numeric proportion
#'
#' Used internally when standardising enrichment dotplot scales.
#'
#' @param x Character vector of ratio strings, or numeric.
#' @return Numeric vector.
ratio_to_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  vapply(
    strsplit(as.character(x), "/"),
    function(z) as.numeric(z[1]) / as.numeric(z[2]),
    numeric(1)
  )
}


#' Source all R module files from the R/ directory
#'
#' Convenience function to load every module in the correct order.
#'
#' @return NULL (invisible). Called for side effects.
source_all_modules <- function() {
  modules <- c("config.R", "utils.R", "load_data.R", "volcano.R",
                "pi_score.R", "concordance.R", "venn.R", "go_enrichment.R")
  for (m in modules) {
    source(here::here("R", m), local = FALSE)
  }
  invisible(NULL)
}
