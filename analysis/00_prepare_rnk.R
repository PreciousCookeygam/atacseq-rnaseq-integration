# =============================================================================
# analysis/00_prepare_rnk.R — Ensure required GSEA rank files are available
# =============================================================================

library(here)

# --- Load modules -------------------------------------------------------------
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "rnk_integration.R"))

status <- ensure_rnk_files(allow_generate = TRUE)

cat(sprintf("RNA rank:  %s [%s]\n", status$rna_path,
            ifelse(status$rna_exists, "OK", "MISSING")))
cat(sprintf("ATAC rank: %s [%s]\n", status$atac_path,
            ifelse(status$atac_exists, "OK", "MISSING")))

if (length(status$messages)) {
  cat(paste(status$messages, collapse = "\n"), "\n")
}

message("[00_prepare_rnk] Rank-file preparation complete.")
