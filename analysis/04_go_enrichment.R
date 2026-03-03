# =============================================================================
# analysis/04_go_enrichment.R — GO enrichment for core gene sets
# =============================================================================
#
# Runs GO:BP enrichment for the four core gene sets (ATAC Up/Down, RNA Up/Down)
# from the Venn overlap analysis.  Then builds standardised dotplots with
# shared axes / colour bars for direct visual comparison.
#
# Prerequisites:
#   - 03_venns.R must have been run first (or objects are built here).
#
# Run:
#   source(here::here("analysis", "04_go_enrichment.R"))
# =============================================================================

# --- Load modules -------------------------------------------------------------
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "load_data.R"))
source(here::here("R", "venn.R"))
source(here::here("R", "go_enrichment.R"))

# --- Close stale graphics devices (Windows PNG lock workaround) ---------------
graphics.off()

# --- Ensure output directory exists -------------------------------------------
outdir <- ensure_output_dir(OUTPUT_DIR)

# --- Load data (once, from cache) ---------------------------------------------
raw_tcga <- load_annotated_results(FILE_TCGA)
raw_lump <- load_annotated_results(FILE_LUMP)

# =============================================================================
# 1) Build core gene sets (mirrors 03_venns.R logic)
# =============================================================================
atac_tcga <- prepare_atac_data(raw_tcga)
atac_lump <- prepare_atac_data(raw_lump)
rna_tcga  <- prepare_rna_data(raw_tcga)
rna_lump  <- prepare_rna_data(raw_lump)

atac_sets <- build_gene_sets(atac_tcga, atac_lump,
                             fc_col = "log2FoldChange", pval_col = "padj")
rna_sets  <- build_gene_sets(rna_tcga, rna_lump,
                             fc_col = "RNAseq_Log2FC", pval_col = "RNAseq_qval")

overlap_up_atac   <- intersect(atac_sets$up_a,   atac_sets$up_b)
overlap_down_atac <- intersect(atac_sets$down_a, atac_sets$down_b)
overlap_up_rna    <- intersect(rna_sets$up_a,    rna_sets$up_b)
overlap_down_rna  <- intersect(rna_sets$down_a,  rna_sets$down_b)

# Background universes (per modality)
bg_atac <- unique(c(atac_tcga$Gene_Name, atac_lump$Gene_Name))
bg_rna  <- unique(c(rna_tcga$Gene_Name,  rna_lump$Gene_Name))

# =============================================================================
# 2) Run GO enrichment (each gene set called exactly once)
# =============================================================================
ego_atac_up   <- go_enrich_and_plot(overlap_up_atac,   "ATAC_core_UP",
                                    universe_symbols = bg_atac, outdir_path = outdir)
ego_atac_down <- go_enrich_and_plot(overlap_down_atac, "ATAC_core_DOWN",
                                    universe_symbols = bg_atac, outdir_path = outdir)
ego_rna_up    <- go_enrich_and_plot(overlap_up_rna,    "RNA_core_UP",
                                    universe_symbols = bg_rna, outdir_path = outdir)
ego_rna_down  <- go_enrich_and_plot(overlap_down_rna,  "RNA_core_DOWN",
                                    universe_symbols = bg_rna, outdir_path = outdir)

cat(sprintf("\n[GO] All enrichment outputs written to: %s\n", outdir))

# =============================================================================
# 3) Standardised dotplots (shared scales for publication)
# =============================================================================
ego_list <- list(ego_atac_up, ego_atac_down, ego_rna_up, ego_rna_down)
shared   <- compute_shared_scales(ego_list)

make_standardised_dotplot(
  ego_atac_up, "ATAC core UP \u2014 GO:BP", shared,
  outdir_path = outdir, filename = "ATAC_core_UP_GO_BP_std.png"
)

make_standardised_dotplot(
  ego_atac_down, "ATAC core DOWN \u2014 GO:BP", shared,
  outdir_path = outdir, filename = "ATAC_core_DOWN_GO_BP_std.png"
)

make_standardised_dotplot(
  ego_rna_up, "RNA core UP \u2014 GO:BP", shared,
  outdir_path = outdir, filename = "RNA_core_UP_GO_BP_std.png"
)

make_standardised_dotplot(
  ego_rna_down, "RNA core DOWN \u2014 GO:BP", shared,
  outdir_path = outdir, filename = "RNA_core_DOWN_GO_BP_std.png"
)

message("\n[04_go_enrichment] GO enrichment and standardised dotplots complete.")
