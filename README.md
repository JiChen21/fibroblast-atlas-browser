# Fibroblast Atlas Browser

A lightweight, production-friendly Streamlit app for exploring a single-cell AnnData (`.h5ad`) dataset.

## Features

- Precomputed UMAP visualization from `adata.obsm["X_umap"]`.
- UMAP coloring by metadata (`adata.obs`).
- On-demand single-gene expression overlay.
- Structured interface with dedicated Atlas Overview / Metadata Explore / Gene Query / Disease–Subtype Compare modules.
- Contextual sidebar controls are shown only in analysis modules (hidden on Home).
- Gene Query module uses a compact two-row layout to reduce scrolling.
- Gene Query charts enforce fixed display orders for `cell_type` and `condition`.
- Multi-select metadata filters for common fibroblast atlas fields.
- Optional downsampling for responsive large-scale plotting.
- Data dictionary tab (`obs` column names + dtypes).
- Demo/mock mode when the real dataset file is missing.

## Environment variables

- `H5AD_PATH` (default: `./data/FBs_adata.h5ad`)
- `UMAP_MAX_POINTS` (default: `200000`)
- `STRICT_DATA` (default: `false`; when `true`, app fails fast instead of using demo fallback)
- `LOG_LEVEL` (default: `INFO`)
- `HOME_IMAGE_PATH` (optional; if not set, app auto-checks `./assets/home_overview.png` then `/data/chenji/fibroblast-atlas-browser-new/plot/home_page.jpg`)

> Do **not** commit real `.h5ad` files to GitHub.

## Local setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Set dataset path:

```bash
export H5AD_PATH=./data/FBs_adata.h5ad
export UMAP_MAX_POINTS=200000  # optional
```

The app reads dataset path from `H5AD_PATH` (environment variable) and does not expose a path picker in the UI.

Run app:

```bash
streamlit run app.py
```

Run on Linux server (public binding):

```bash
streamlit run app.py --server.address 0.0.0.0 --server.port 8501
```

## Smoke test / validation flow

After installing requirements, run this validation command:

```bash
python scripts/validate_h5ad.py --path "$H5AD_PATH"
```

What it validates:

- H5AD file exists and can be loaded
- `adata.obsm["X_umap"]` exists
- Required `obs` columns exist by default (`condition`, `region`, `cell_type`)

You can override required metadata fields via `--require-obs`, for example:

```bash
python scripts/validate_h5ad.py --path "$H5AD_PATH" --require-obs condition cell_type
```

This is a minimal smoke test to verify dataset compatibility before launching Streamlit.

## Slim large `.h5ad` for this browser

If your source AnnData contains analysis artifacts not needed by this app, create a browser-focused slim file:

```bash
python scripts/slim_h5ad_for_browser.py \
  --input ./data/FBs_adata.h5ad \
  --output ./data/FBs_adata.slim.h5ad \
  --keep-obs condition region cell_type
```

Default keep/drop behavior:

- Keep: `X` (cast to float32), `obs` (`condition`, `region`, `cell_type`), `var_names`, `obsm["X_umap"]`
- Drop: `raw`, `layers`, `uns`, `obsp`, `varm`, extra `obsm` keys

## Tests

```bash
pytest -q
python -m py_compile app.py scripts/validate_h5ad.py core.py
```

## Add an image to Home page

1. Recommended: put your image file at `./assets/home_overview.png`.
2. Or set `HOME_IMAGE_PATH` to any absolute path (for example `/data/chenji/fibroblast-atlas-browser-new/plot/home_page.jpg`).

```bash
export HOME_IMAGE_PATH=./assets/home_overview.png
```

If `HOME_IMAGE_PATH` is not set, the app will automatically check:
- `./assets/home_overview.png`
- `/data/chenji/fibroblast-atlas-browser-new/plot/home_page.jpg`

If one exists, the app shows it on **Home** below the atlas overview text.

## Docker

Build image:

```bash
docker build -t fibroblast-atlas-browser .
```

Run container:

```bash
docker run --rm -p 8501:8501 \
  -e H5AD_PATH=/data/FBs_adata.h5ad \
  -e UMAP_MAX_POINTS=200000 \
  -v $(pwd)/data:/data \
  fibroblast-atlas-browser
```

## Performance notes for ~740k cells

For ~739,910 cells, full interactive plotting can be slow depending on hardware.

- `Auto` mode (default): plots all filtered points up to cap, otherwise stratified downsamples.
- `Downsampled` mode: always uses stratified subset up to cap.
- `Full` mode: plots all filtered cells (best fidelity, slowest).

Tradeoff:

- Downsampling improves UI responsiveness but does not show every point.
- Full mode preserves all points but may lag on commodity hardware.

The app always avoids densifying the entire expression matrix and computes expression only for the selected gene.

## Expression source compatibility

Gene query uses normalized expression from `X` only.

Gene lookup uses `adata.var_names`.

## Expected AnnData structure

- `adata.obsm["X_umap"]` present
- Required metadata in `adata.obs`
- Gene IDs in `adata.var_names` (default lookup)
