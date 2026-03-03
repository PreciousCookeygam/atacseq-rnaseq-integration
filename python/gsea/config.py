from __future__ import annotations

from pathlib import Path

# Repository root (python/gsea/config.py -> python -> repo)
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "output" / "gsea"

# Default .rnk paths (matches R config expectations)
RNA_RNK_PATH = DATA_DIR / "rna_seq_rank.rnk"
ATAC_RNK_PATH = DATA_DIR / "atac_seq_rank.rnk"

# Default inputs
MERGED_ATAC_RNA_PATH = DATA_DIR / "Merged_ATAC_RNA_results.tsv"

# Optional GMT defaults (update if you relocate files)
CGP_GMT_PATH = Path(r"C:\Users\User\Downloads\msigdb_v2025.1.Hs_files_to_download_locally\msigdb_v2025.1.Hs_files_to_download_locally\msigdb_v2025.1.Hs_GMTs\c2.cgp.v2025.1.Hs.symbols.gmt")
HALLMARK_GMT_PATH = Path(r"C:\Users\User\Downloads\msigdb_v2025.1.Hs_files_to_download_locally\msigdb_v2025.1.Hs_files_to_download_locally\msigdb_v2025.1.Hs_GMTs\h.all.v2025.1.Hs.symbols.gmt")

# Default columns (override via CLI if needed)
GENE_COL = "Gene_Name"
RNA_LOG2FC_COL = "RNAseq_Log2FC"
RNA_QVAL_COL = "RNAseq_qval"
ATAC_LOG2FC_COL = "log2FC"
ATAC_PVAL_COL = "pval"
