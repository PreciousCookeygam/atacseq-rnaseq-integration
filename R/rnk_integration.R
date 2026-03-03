# =============================================================================
# rnk_integration.R — Ensure Python-generated rank files are available
# =============================================================================

#' Ensure required .rnk files exist (generate via Python if needed)
#'
#' Checks for `data/rna_seq_rank.rnk` and `data/atac_seq_rank.rnk`. If either
#' is missing and `allow_generate = TRUE`, calls the Python prerank script to
#' generate them from merged data.
#'
#' @param allow_generate Logical. If TRUE, attempt Python generation when
#'   missing files are detected.
#' @param python_exe Character. Python executable path/command.
#' @param python_script Character. Path to run_gsea_prerank.py.
#' @return A list with file paths, existence flags, generation status, and
#'   log messages.
ensure_rnk_files <- function(allow_generate = TRUE,
                             python_exe = PYTHON_EXE,
                             python_script = PYTHON_GSEA_SCRIPT) {
  result <- list(
    rna_path = FILE_RNA_RANK,
    atac_path = FILE_ATAC_RANK,
    rna_exists = file.exists(FILE_RNA_RANK),
    atac_exists = file.exists(FILE_ATAC_RANK),
    generated = FALSE,
    messages = character()
  )

  add_msg <- function(msg) {
    result$messages <<- c(result$messages, msg)
  }

  if (result$rna_exists && result$atac_exists) {
    add_msg("Both .rnk files already exist.")
    return(result)
  }

  missing_now <- c(
    if (!result$rna_exists) FILE_RNA_RANK,
    if (!result$atac_exists) FILE_ATAC_RANK
  )

  add_msg(sprintf("Missing .rnk files: %s", paste(missing_now, collapse = ", ")))

  if (!allow_generate) {
    stop(
      paste(
        "Required .rnk files are missing and generation is disabled.",
        "Set `allow_generate = TRUE` or place the files in data/:",
        paste(missing_now, collapse = ", "),
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  if (!file.exists(python_script)) {
    stop(
      paste(
        "Python prerank script not found:",
        python_script,
        "Check `PYTHON_GSEA_SCRIPT` in R/config.R.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  py_dir <- normalizePath(here::here("python"), winslash = "/", mustWork = FALSE)
  py_sep <- if (.Platform$OS.type == "windows") ";" else ":"
  py_old <- Sys.getenv("PYTHONPATH", unset = "")
  py_new <- if (nzchar(py_old)) paste(py_dir, py_old, sep = py_sep) else py_dir

  add_msg(sprintf("Generating .rnk files via Python: %s", python_exe))

  cmd_out <- tryCatch(
    system2(
      command = python_exe,
      args = c(shQuote(python_script), "--score-mode", "pi_score"),
      stdout = TRUE,
      stderr = TRUE,
      env = c(sprintf("PYTHONPATH=%s", py_new))
    ),
    error = function(e) {
      structure(
        sprintf("system2 failed: %s", conditionMessage(e)),
        status = 1
      )
    }
  )

  status <- attr(cmd_out, "status")
  if (is.null(status)) status <- 0L

  result$generated <- TRUE
  result$rna_exists <- file.exists(FILE_RNA_RANK)
  result$atac_exists <- file.exists(FILE_ATAC_RANK)

  if (length(cmd_out)) {
    add_msg(paste(cmd_out, collapse = "\n"))
  }

  if (status != 0L || !result$rna_exists || !result$atac_exists) {
    missing_after <- c(
      if (!result$rna_exists) FILE_RNA_RANK,
      if (!result$atac_exists) FILE_ATAC_RANK
    )

    stop(
      paste(
        "Failed to ensure required .rnk files.",
        sprintf("Python exit status: %s", status),
        sprintf("Python executable: %s", python_exe),
        sprintf("Script: %s", python_script),
        sprintf(
          "Command: %s %s --score-mode pi_score",
          python_exe,
          shQuote(python_script)
        ),
        sprintf(
          "Still missing: %s",
          if (length(missing_after)) paste(missing_after, collapse = ", ") else "none"
        ),
        "Action: verify Python deps in python/requirements.txt and set PYTHON_EXE env var.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  add_msg("Successfully ensured both .rnk files.")
  result
}
