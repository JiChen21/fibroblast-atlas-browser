#!/usr/bin/env python3
"""Create a lightweight .h5ad tailored for this Streamlit browser."""

import argparse
import os
from typing import List

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

DEFAULT_OBS_KEEP = ["condition", "region", "cell_type"]
DEFAULT_OBSM_KEEP = ["X_umap"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Slim an AnnData file for fibroblast-atlas-browser usage.")
    parser.add_argument("--input", required=True, help="Input .h5ad path.")
    parser.add_argument("--output", required=True, help="Output slim .h5ad path.")
    parser.add_argument(
        "--keep-obs",
        nargs="+",
        default=DEFAULT_OBS_KEEP,
        help="obs columns to keep (default: condition region cell_type).",
    )
    parser.add_argument(
        "--keep-var-cols",
        nargs="*",
        default=[],
        help="Optional var columns to keep besides var_names index.",
    )
    parser.add_argument(
        "--x-dtype",
        choices=["float32", "float64"],
        default="float32",
        help="dtype for expression matrix X (default: float32).",
    )
    parser.add_argument(
        "--no-categorical-obs",
        action="store_true",
        help="Do not cast kept obs columns to pandas category.",
    )
    return parser.parse_args()


def existing_columns(columns: List[str], available: pd.Index) -> List[str]:
    return [c for c in columns if c in available]


def main() -> int:
    args = parse_args()
    if not os.path.exists(args.input):
        print(f"ERROR: input not found: {args.input}")
        return 1

    print(f"Loading {args.input}")
    adata = ad.read_h5ad(args.input)

    keep_obs = existing_columns(args.keep_obs, adata.obs.columns)
    if not keep_obs:
        print("ERROR: none of --keep-obs columns exist in adata.obs")
        return 2

    keep_var_cols = existing_columns(args.keep_var_cols, adata.var.columns)

    print("Building slim AnnData...")
    obs = adata.obs.loc[:, keep_obs].copy()
    if not args.no_categorical_obs:
        for col in keep_obs:
            obs[col] = obs[col].astype("category")

    var = adata.var.loc[:, keep_var_cols].copy() if keep_var_cols else pd.DataFrame(index=adata.var_names.copy())
    if var.index.name != adata.var_names.name:
        var.index.name = adata.var_names.name

    x = adata.X
    target_dtype = np.float32 if args.x_dtype == "float32" else np.float64
    if sparse.issparse(x):
        x = x.tocsr().astype(target_dtype)
    else:
        x = np.asarray(x, dtype=target_dtype)

    out = ad.AnnData(X=x, obs=obs, var=var)

    for key in DEFAULT_OBSM_KEEP:
        if key in adata.obsm:
            out.obsm[key] = np.asarray(adata.obsm[key], dtype=np.float32)

    if "X_umap" not in out.obsm:
        print("WARNING: X_umap not found in source; output will not be usable by current app.")

    out.write_h5ad(args.output, compression="gzip")

    in_mb = os.path.getsize(args.input) / (1024 * 1024)
    out_mb = os.path.getsize(args.output) / (1024 * 1024)
    ratio = (out_mb / in_mb) if in_mb > 0 else 0

    print("Done.")
    print(f"Input : {in_mb:.2f} MB")
    print(f"Output: {out_mb:.2f} MB")
    print(f"Ratio : {ratio:.3f}")
    print(f"Kept obs columns: {keep_obs}")
    print(f"Kept var columns: {keep_var_cols if keep_var_cols else '[index only]'}")
    print(f"Kept obsm keys: {[k for k in DEFAULT_OBSM_KEEP if k in out.obsm]}")
    print("Dropped by design: raw, layers, uns, obsp, varm, extra obsm keys.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

