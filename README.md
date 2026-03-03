# ATAC-seq / RNA-seq Integrative Analysis — TCGA-BLCA

Integrative analysis of ATAC-seq chromatin accessibility and RNA-seq gene
expression data comparing **TCGA-BLCA** bladder cancer cohorts (All TCGA and
Luminal Papillary subtype) against **normal bladder (P0)** controls.

## Project Structure

```
├── R/                        # Reusable function modules (source these)
│   ├── config.R              # Thresholds, colour palettes, paths
│   ├── utils.R               # General-purpose helpers
│   ├── rnk_integration.R     # Ensure/generate required .rnk files
│   ├── load_data.R           # Cached data loading & column selection
│   ├── volcano.R             # Volcano plot construction
│   ├── pi_score.R            # π-score computation & scatter plots
│   ├── concordance.R         # Cross-modality concordance statistics
│   ├── venn.R                # Venn diagrams & overlap statistics
│   └── go_enrichment.R       # GO enrichment & standardised dotplots
│
├── analysis/                 # Thin orchestration scripts (run these)
│   ├── 00_prepare_rnk.R      # Step 0: ensure/generate ATAC/RNA .rnk files
│   ├── 01_volcanos.R         # ATAC + RNA volcano plots (4 plots)
│   ├── 02_pi_scores.R        # π-score plots + concordance stats
│   ├── 03_venns.R            # Venn diagrams + core gene-set exports
│   ├── 04_go_enrichment.R    # GO:BP enrichment + publication dotplots
│   └── run_all.R             # Master script — runs steps 00–04
│
├── python/                   # Python package for .rnk generation / GSEApy
│   ├── gsea/
│   │   ├── config.py
│   │   ├── io.py
│   │   └── prerank.py
│   ├── scripts/
│   │   └── run_gsea_prerank.py
│   └── requirements.txt
│
├── data/                     # Input data files (not version-controlled)
│   ├── ATACseq_annotated_results_TCGA-BLCA_vs_P0-Bladder.tsv
│   ├── ATACseq_annotated_results_TCGA-BLCA-LumP_vs_P0-Bladder.tsv
│   ├── ATACseq_annotated_results_TCGA-BLCA-LumP_vs_TCGA-BLCA-BaSq.tsv
│   ├── Merged_ATAC_RNA_results.tsv
│   ├── concordant_core_152_genes.tsv
│   ├── rna_seq_rank.rnk
│   └── atac_seq_rank.rnk
│
├── output/                   # Generated plots and tables (gitignored)
│
├── DESCRIPTION               # Project metadata & dependency list
├── atacseq-rnaseq-integration.Rproj  # RStudio project file
├── .gitignore
└── README.md                 # This file
```

## Quick Start

### 1. Install dependencies

**CRAN packages:**

```r
install.packages(c("dplyr", "readr", "ggplot2", "ggrepel",
                    "VennDiagram", "scales", "here"))
```

**Bioconductor packages:**

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
```

**Python packages (for automatic `.rnk` generation):**

```bash
pip install -r python/requirements.txt
```

### 2. Configure Python executable (optional)

By default, R uses `python`. You can override this with environment variable
`PYTHON_EXE`.

**Windows (PowerShell):**

```powershell
$env:PYTHON_EXE="C:\\path\\to\\python.exe"
```

### 3. Place data files

Copy or move the input TSV files into the `data/` directory. If `.rnk` files
are missing, Step 0 auto-generates them from `Merged_ATAC_RNA_results.tsv`.
Expected names are listed in `R/config.R`.

### 4. Run the full pipeline

From the project root:

```r
source("analysis/run_all.R")
```

Or run individual steps:

```r
library(here)
source(here::here("analysis", "00_prepare_rnk.R"))
source(here::here("analysis", "01_volcanos.R"))
source(here::here("analysis", "02_pi_scores.R"))
source(here::here("analysis", "03_venns.R"))
source(here::here("analysis", "04_go_enrichment.R"))
```

All outputs (PNG plots, TSV tables) are written to `output/`.

## Input Data Schema

Each annotated results TSV contains one row per **peak–gene pair** with
~25 columns.  The key columns used by this pipeline:

| Column             | Description                          |
|--------------------|--------------------------------------|
| `Gene_Name`        | HGNC gene symbol                     |
| `log2FoldChange`   | ATAC-seq log2 fold-change            |
| `padj`             | ATAC-seq BH-adjusted p-value         |
| `RNAseq_Log2FC`    | RNA-seq log2 fold-change             |
| `RNAseq_qval`      | RNA-seq q-value (FDR)                |

## Analysis Workflow

| Step | Script              | What it produces                                      |
|------|---------------------|-------------------------------------------------------|
| 0    | `00_prepare_rnk.R`  | Validates rank files; auto-generates missing `.rnk`   |
| 1    | `01_volcanos.R`     | 4 volcano plots (ATAC/RNA × TCGA/LumP)               |
| 2    | `02_pi_scores.R`    | 2 π-score scatter plots + concordance stats + exports |
| 3    | `03_venns.R`        | 4 Venn diagrams + 4 core gene TSVs + manuscript text  |
| 4    | `04_go_enrichment.R`| GO:BP enrichment tables, dotplots, and barplots       |

## Example Output Figures

These figures are intended to answer three questions:

1. **What changes are significant within each cohort comparison?** (Volcano plots)
2. **How consistent are ATAC and RNA effects for the same genes?** (π-score plots)
3. **Which dysregulated genes are shared between cohorts, and what biology do they represent?** (Venn + GO plots)

### Volcano plots

Each volcano plot shows differential signal for one modality and one contrast:

- **x-axis**: log2 fold-change ($\log_2FC$), where right = upregulated and left = downregulated.
- **y-axis**: $-\log_{10}(p_{adj})$, where higher points are more statistically significant.
- Vertical and horizontal threshold lines correspond to the pipeline cut-offs ($|\log_2FC| \ge 1$ and $p_{adj} < 0.05$).

Interpretation:

- Points in the **upper-right** quadrant are significantly upregulated.
- Points in the **upper-left** quadrant are significantly downregulated.
- Comparing TCGA vs LumP volcanoes highlights subtype-specific versus broadly shared effects.

![ATAC all TCGA volcano](docs/images/volcano_ATAC_allTCGA_vs_P0.png)
![ATAC LumP volcano](docs/images/volcano_ATAC_LumP_vs_P0.png)
![RNA all TCGA volcano](docs/images/volcano_RNA_allTCGA_vs_P0.png)
![RNA LumP volcano](docs/images/volcano_RNA_LumP_vs_P0.png)

### π-score concordance

π-score integrates effect size and significance to compare ATAC and RNA changes per gene.

- **x-axis**: ATAC π-score
- **y-axis**: RNA π-score
- Each point represents a gene present in both modalities.

Interpretation:

- Genes near the **positive diagonal** show concordant activation (both ATAC and RNA increased).
- Genes near the **negative diagonal** show concordant repression (both decreased).
- Opposite-sign quadrants indicate discordant regulation (accessibility and expression moving in different directions).

![PI score all TCGA](docs/images/pi_score_allTCGA_vs_P0.png)
![PI score LumP](docs/images/pi_score_LumP_vs_P0.png)

### Venn overlap plots

Venn diagrams compare significant gene sets between cohorts (All TCGA vs LumP), separately for upregulated and downregulated genes in each modality.

Interpretation:

- **Intersection region** = shared core dysregulated genes across cohorts.
- **Non-overlapping regions** = cohort-specific dysregulation.
- These intersections are exported as core gene lists and used for downstream enrichment.

![ATAC up Venn](docs/images/venn_ATAC_Up_TCGA_vs_LumP.png)
![ATAC down Venn](docs/images/venn_ATAC_Down_TCGA_vs_LumP.png)
![RNA up Venn](docs/images/venn_RNA_Up_TCGA_vs_LumP.png)
![RNA down Venn](docs/images/venn_RNA_Down_TCGA_vs_LumP.png)

### Standardised GO dotplots

GO:BP dotplots summarize functional enrichment of the core overlapping gene sets.

- **y-axis**: enriched biological process term.
- **x-axis**: standardized gene ratio (comparable across contrasts).
- **dot size**: number of genes mapped to the term.
- **dot color**: enrichment significance (more intense color = lower adjusted p-value).

Interpretation:

- Large, strongly colored dots to the right indicate robust and substantial pathway enrichment.
- Comparing ATAC and RNA panels reveals whether chromatin-level and transcription-level signals converge on similar biology.

![ATAC core up GO](docs/images/ATAC_core_UP_GO_BP_std.png)
![ATAC core down GO](docs/images/ATAC_core_DOWN_GO_BP_std.png)
![RNA core up GO](docs/images/RNA_core_UP_GO_BP_std.png)
![RNA core down GO](docs/images/RNA_core_DOWN_GO_BP_std.png)

## Thresholds

All statistical thresholds are centralised in [`R/config.R`](R/config.R):

- **|log2FC| ≥ 1** — fold-change cut-off
- **adjusted p-value < 0.05** — significance cut-off

## Original Scripts

The original monolithic scripts (`rscript_atac_rna_plot.R` and
`venn_diagram.R`) are retained in the project root for reference.  They are
superseded by the modular `R/` + `analysis/` structure.
