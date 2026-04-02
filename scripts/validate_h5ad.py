#!/usr/bin/env python3
"""Minimal dataset validation script for AnnData portal."""

import argparse
import os
import sys

import anndata as ad

DEFAULT_REQUIRED_OBS = [
    "condition",
    "region",
    "cell_type",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate H5AD compatibility for the Streamlit app")
    parser.add_argument("--path", default=os.getenv("H5AD_PATH", "./data/FBs_adata.h5ad"))
    parser.add_argument("--require-obs", nargs="*", default=DEFAULT_REQUIRED_OBS)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not os.path.exists(args.path):
        print(f"ERROR: H5AD file missing at: {args.path}")
        return 1

    try:
        adata = ad.read_h5ad(args.path)
    except Exception as exc:
        print(f"ERROR: failed to load H5AD: {exc}")
        return 1

    if "X_umap" not in adata.obsm:
        print("ERROR: missing required embedding adata.obsm['X_umap']")
        return 1

    missing_obs = [col for col in args.require_obs if col not in adata.obs.columns]
    if missing_obs:
        print("ERROR: missing required obs columns:", ", ".join(missing_obs))
        return 1

    print("OK: dataset validation passed")
    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")
    print(f"obs columns: {len(adata.obs.columns)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
