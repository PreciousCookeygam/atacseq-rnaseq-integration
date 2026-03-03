from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def read_table(path: Path | str, sep: str = "\t") -> pd.DataFrame:
    # Read all columns as strings first for robust parsing of mixed-type TSVs.
    # Downstream code explicitly coerces required numeric columns.
    return pd.read_csv(path, sep=sep, dtype=str, low_memory=False)


def compute_pi_score(log2fc: pd.Series, pvals: pd.Series, eps: float = 1e-300) -> pd.Series:
    pvals = pd.to_numeric(pvals, errors="coerce").clip(lower=eps)
    log2fc = pd.to_numeric(log2fc, errors="coerce")
    return log2fc * (-np.log10(pvals))


def build_rank_series(
    df: pd.DataFrame,
    gene_col: str,
    score_col: str,
    drop_duplicates: bool = True,
) -> pd.Series:
    r = df[[gene_col, score_col]].dropna().copy()
    if drop_duplicates:
        r = r.drop_duplicates(subset=[gene_col])
    r = r.rename(columns={gene_col: "gene", score_col: "score"})
    r = r.sort_values("score", ascending=False)
    return pd.Series(r["score"].values, index=r["gene"].values)


def write_rnk(series: pd.Series, path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    series.to_csv(path, sep="\t", header=False)


def coerce_columns(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df
