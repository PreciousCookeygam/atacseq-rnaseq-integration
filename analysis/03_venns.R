# =============================================================================
# analysis/03_venns.R — Venn diagrams for ATAC-seq and RNA-seq consistency
# =============================================================================
#
# Produces:
#   1. ATAC Up Venn   (TCGA vs LumP overlap)
#   2. ATAC Down Venn (TCGA vs LumP overlap)
#   3. RNA Up Venn    (TCGA vs LumP overlap)
#   4. RNA Down Venn  (TCGA vs LumP overlap)
#   5. Manuscript-ready text snippets
#   6. Core gene-list TSV exports
#
# Run:
#   source(here::here("analysis", "03_venns.R"))
# =============================================================================

# --- Load modules -------------------------------------------------------------
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "load_data.R"))
source(here::here("R", "venn.R"))

# --- Ensure output directory exists -------------------------------------------
ensure_output_dir(OUTPUT_DIR)

# --- Load data (once) ---------------------------------------------------------
raw_tcga <- load_annotated_results(FILE_TCGA)
raw_lump <- load_annotated_results(FILE_LUMP)

# =============================================================================
# ATAC-seq Venns
# =============================================================================
atac_tcga <- prepare_atac_data(raw_tcga)
atac_lump <- prepare_atac_data(raw_lump)

atac_sets <- build_gene_sets(
  atac_tcga, atac_lump,
  fc_col   = "log2FoldChange",
  pval_col = "padj"
)

# Overlap stats
atac_up_stats   <- compute_overlap_stats(atac_sets$up_a,   atac_sets$up_b)
atac_down_stats <- compute_overlap_stats(atac_sets$down_a, atac_sets$down_b)

print_overlap_summary(atac_up_stats, atac_down_stats, modality = "ATAC-seq")

# Venn: ATAC Up
draw_venn(
  atac_sets$up_a, atac_sets$up_b,
  labels      = c("TCGA Up", "LumP Up"),
  title       = "ATAC: Up genes \u2014 TCGA vs LumP overlap",
  output_file = "venn_ATAC_Up_TCGA_vs_LumP.png"
)

# Venn: ATAC Down
draw_venn(
  atac_sets$down_a, atac_sets$down_b,
  labels      = c("TCGA Down", "LumP Down"),
  title       = "ATAC: Down genes \u2014 TCGA vs LumP overlap",
  cat_pos     = c(-20, 20),
  cat_dist    = c(0.03, 0.03),
  output_file = "venn_ATAC_Down_TCGA_vs_LumP.png"
)

# Export core ATAC lists
readr::write_tsv(
  data.frame(Gene = atac_up_stats$overlap),
  file.path(OUTPUT_DIR, "ATAC_core_UP_TCGA_and_LumP.tsv")
)
readr::write_tsv(
  data.frame(Gene = atac_down_stats$overlap),
  file.path(OUTPUT_DIR, "ATAC_core_DOWN_TCGA_and_LumP.tsv")
)

# Manuscript text
generate_manuscript_text(atac_up_stats, atac_down_stats, modality = "ATAC-seq")

# =============================================================================
# RNA-seq Venns
# =============================================================================
rna_tcga <- prepare_rna_data(raw_tcga)
rna_lump <- prepare_rna_data(raw_lump)

rna_sets <- build_gene_sets(
  rna_tcga, rna_lump,
  fc_col   = "RNAseq_Log2FC",
  pval_col = "RNAseq_qval"
)

# Overlap stats
rna_up_stats   <- compute_overlap_stats(rna_sets$up_a,   rna_sets$up_b)
rna_down_stats <- compute_overlap_stats(rna_sets$down_a, rna_sets$down_b)

print_overlap_summary(rna_up_stats, rna_down_stats, modality = "RNA-seq")

# Venn: RNA Up
draw_venn(
  rna_sets$up_a, rna_sets$up_b,
  labels      = c("TCGA Up", "LumP Up"),
  title       = "RNA: Up genes \u2014 TCGA vs LumP overlap",
  output_file = "venn_RNA_Up_TCGA_vs_LumP.png"
)

# Venn: RNA Down
draw_venn(
  rna_sets$down_a, rna_sets$down_b,
  labels      = c("TCGA Down", "LumP Down"),
  title       = "RNA: Down genes \u2014 TCGA vs LumP overlap",
  output_file = "venn_RNA_Down_TCGA_vs_LumP.png"
)

# Export core RNA lists
readr::write_tsv(
  data.frame(Gene = rna_up_stats$overlap),
  file.path(OUTPUT_DIR, "RNA_core_UP_TCGA_and_LumP.tsv")
)
readr::write_tsv(
  data.frame(Gene = rna_down_stats$overlap),
  file.path(OUTPUT_DIR, "RNA_core_DOWN_TCGA_and_LumP.tsv")
)

# Manuscript text
generate_manuscript_text(rna_up_stats, rna_down_stats, modality = "RNA-seq")

message("\n[03_venns] All Venn diagrams and core lists complete.")
