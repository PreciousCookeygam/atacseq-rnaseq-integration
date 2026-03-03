from __future__ import annotations

import argparse
import multiprocessing
from pathlib import Path

import pandas as pd

from gsea.config import (
    MERGED_ATAC_RNA_PATH,
    GENE_COL,
    RNA_LOG2FC_COL,
    RNA_QVAL_COL,
    ATAC_LOG2FC_COL,
    ATAC_PVAL_COL,
    RNA_RNK_PATH,
    ATAC_RNK_PATH,
    OUTPUT_DIR,
    CGP_GMT_PATH,
    HALLMARK_GMT_PATH,
)
from gsea.io import compute_pi_score, build_rank_series, read_table, write_rnk
from gsea.prerank import run_prerank


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate .rnk files and optionally run GSEApy prerank."
    )
    parser.add_argument(
        "--input",
        default=str(MERGED_ATAC_RNA_PATH),
        help="Path to merged ATAC/RNA TSV file.",
    )
    parser.add_argument("--gene-col", default=GENE_COL)
    parser.add_argument("--rna-log2fc-col", default=RNA_LOG2FC_COL)
    parser.add_argument("--rna-qval-col", default=RNA_QVAL_COL)
    parser.add_argument("--atac-log2fc-col", default=ATAC_LOG2FC_COL)
    parser.add_argument("--atac-pval-col", default=ATAC_PVAL_COL)

    parser.add_argument(
        "--rna-rnk",
        default=str(RNA_RNK_PATH),
        help="Output path for RNA .rnk file.",
    )
    parser.add_argument(
        "--atac-rnk",
        default=str(ATAC_RNK_PATH),
        help="Output path for ATAC .rnk file.",
    )

    parser.add_argument(
        "--score-mode",
        choices=["pi_score", "log2fc"],
        default="pi_score",
        help="Score type for ranking.",
    )

    parser.add_argument(
        "--run-gsea",
        action="store_true",
        help="Run GSEApy prerank after writing .rnk files.",
    )
    parser.add_argument(
        "--use-existing-rnk",
        action="store_true",
        help="Skip TSV parsing and use existing --rna-rnk/--atac-rnk files directly.",
    )
    parser.add_argument(
        "--cgp-gmt",
        default=str(CGP_GMT_PATH),
        help="Path to CGP GMT file.",
    )
    parser.add_argument(
        "--hallmark-gmt",
        default=str(HALLMARK_GMT_PATH),
        help="Path to Hallmark GMT file.",
    )
    parser.add_argument(
        "--outdir",
        default=str(OUTPUT_DIR),
        help="Base output directory for GSEApy results.",
    )
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--permutation-num", type=int, default=1000)

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.use_existing_rnk:
        rna_df = pd.read_csv(args.rna_rnk, sep="\t", header=None, names=["gene", "score"])
        atac_df = pd.read_csv(args.atac_rnk, sep="\t", header=None, names=["gene", "score"])
        rna_rnk = pd.Series(rna_df["score"].values, index=rna_df["gene"].values)
        atac_rnk = pd.Series(atac_df["score"].values, index=atac_df["gene"].values)
        print("Using existing .rnk files:")
        print(f" - {args.rna_rnk}")
        print(f" - {args.atac_rnk}")
    else:
        df = read_table(args.input)
        gene_col = args.gene_col

        if gene_col not in df.columns:
            raise ValueError(f"Missing gene column: {gene_col}")

        # RNA scores
        if args.score_mode == "pi_score":
            if args.rna_log2fc_col not in df.columns or args.rna_qval_col not in df.columns:
                raise ValueError(
                    "Missing RNA columns for pi_score: "
                    f"{args.rna_log2fc_col}, {args.rna_qval_col}"
                )
            df["rna_score"] = compute_pi_score(
                df[args.rna_log2fc_col], df[args.rna_qval_col]
            )
        else:
            if args.rna_log2fc_col not in df.columns:
                raise ValueError(f"Missing RNA log2FC column: {args.rna_log2fc_col}")
            df["rna_score"] = pd.to_numeric(df[args.rna_log2fc_col], errors="coerce")

        # ATAC scores
        if args.score_mode == "pi_score":
            if args.atac_log2fc_col in df.columns and args.atac_pval_col in df.columns:
                df["atac_score"] = compute_pi_score(
                    df[args.atac_log2fc_col], df[args.atac_pval_col]
                )
            else:
                df["atac_score"] = pd.to_numeric(df.get(args.atac_log2fc_col), errors="coerce")
        else:
            df["atac_score"] = pd.to_numeric(df.get(args.atac_log2fc_col), errors="coerce")

        rna_rnk = build_rank_series(df, gene_col, "rna_score")
        atac_rnk = build_rank_series(df, gene_col, "atac_score")

        write_rnk(rna_rnk, args.rna_rnk)
        write_rnk(atac_rnk, args.atac_rnk)

        print("Saved .rnk files:")
        print(f" - {args.rna_rnk}")
        print(f" - {args.atac_rnk}")

    if not args.run_gsea:
        return

    threads = args.threads or max(1, multiprocessing.cpu_count() - 1)
    outdir = Path(args.outdir)

    if args.cgp_gmt:
        run_prerank(
            rna_rnk,
            args.cgp_gmt,
            outdir / "gsea_cgp_RNA_results",
            threads=threads,
            permutation_num=args.permutation_num,
        )
        run_prerank(
            atac_rnk,
            args.cgp_gmt,
            outdir / "gsea_cgp_ATAC_results",
            threads=threads,
            permutation_num=args.permutation_num,
        )

    if args.hallmark_gmt:
        run_prerank(
            rna_rnk,
            args.hallmark_gmt,
            outdir / "gsea_hall_RNA_results",
            threads=threads,
            permutation_num=args.permutation_num,
        )
        run_prerank(
            atac_rnk,
            args.hallmark_gmt,
            outdir / "gsea_hall_ATAC_results",
            threads=threads,
            permutation_num=args.permutation_num,
        )


if __name__ == "__main__":
    main()
