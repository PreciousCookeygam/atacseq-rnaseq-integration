from __future__ import annotations

from pathlib import Path
from typing import Iterable

import gseapy as gp
import pandas as pd


def run_prerank(
    rnk: pd.Series | pd.DataFrame,
    gene_sets: str | Path | Iterable[str],
    outdir: str | Path,
    min_size: int = 15,
    max_size: int = 500,
    permutation_num: int = 1000,
    seed: int = 42,
    threads: int = 1,
    verbose: bool = True,
):
    outdir = str(outdir)
    return gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        outdir=outdir,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        seed=seed,
        threads=threads,
        verbose=verbose,
    )
