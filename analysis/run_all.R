# =============================================================================
# analysis/run_all.R — Master orchestration script
# =============================================================================
#
# Sources every analysis step in order.  Run this to reproduce the full
# analysis from a clean R session:
#
#   Rscript analysis/run_all.R
#
# Or interactively:
#   source(here::here("analysis", "run_all.R"))
# =============================================================================

library(here)

cat("=================================================================\n")
cat(" ATAC-seq / RNA-seq Integrative Analysis Pipeline\n")
cat("=================================================================\n\n")

t0 <- Sys.time()

# Step 1 — Volcano plots (ATAC + RNA, TCGA + LumP)
cat("[1/4] Generating volcano plots ...\n")
source(here::here("analysis", "01_volcanos.R"), local = new.env())

# Step 2 — π-score concordance & statistics
cat("\n[2/4] Computing π-scores and concordance ...\n")
source(here::here("analysis", "02_pi_scores.R"), local = new.env())

# Step 3 — Venn diagrams & core gene-set exports
cat("\n[3/4] Building Venn diagrams ...\n")
source(here::here("analysis", "03_venns.R"), local = new.env())

# Step 4 — GO enrichment & standardised dotplots
cat("\n[4/4] Running GO enrichment ...\n")
source(here::here("analysis", "04_go_enrichment.R"), local = new.env())

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)

cat("\n=================================================================\n")
cat(sprintf(" Pipeline complete in %s min.\n", elapsed))
cat(sprintf(" Outputs saved to: %s\n", here::here("output")))
cat("=================================================================\n")
