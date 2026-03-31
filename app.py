import os
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from scipy import sparse

DEFAULT_H5AD_PATH = "./data/FBs_adata.h5ad"
DEFAULT_MAX_POINTS = int(os.getenv("UMAP_MAX_POINTS", "200000"))
FILTER_COLUMNS = [
    "condition",
    "region",
    "donor",
    "dataset_name",
    "cell_type",
    "cell_type_V2",
    "cell_type_V3",
    "predicted_cell_type",
    "leiden",
    "leiden_default",
]


def build_mock_adata(n_cells: int = 3000, n_genes: int = 80) -> ad.AnnData:
    """Create a lightweight synthetic AnnData so the UI can run without a real file."""
    rng = np.random.default_rng(42)
    gene_names = np.array([f"GENE_{i:03d}" for i in range(n_genes)])

    counts = rng.poisson(lam=0.3, size=(n_cells, n_genes)).astype(np.float32)
    counts[counts < 1] = 0
    x = sparse.csr_matrix(counts)

    obs = pd.DataFrame(
        {
            "batch": rng.choice(["batch1", "batch2"], size=n_cells),
            "cell_type": rng.choice(["FB", "Immune", "Endothelial"], size=n_cells),
            "donor": rng.choice(["D1", "D2", "D3", "D4"], size=n_cells),
            "region": rng.choice(["Skin", "Lung", "Gut"], size=n_cells),
            "dataset_name": rng.choice(["CohortA", "CohortB"], size=n_cells),
            "condition": rng.choice(["Control", "Disease"], size=n_cells),
            "age": rng.integers(20, 80, size=n_cells).astype(str),
            "gender": rng.choice(["F", "M"], size=n_cells),
            "predicted_cell_type": rng.choice(["FB_1", "FB_2", "FB_3"], size=n_cells),
            "condition_v1": rng.choice(["Healthy", "Inflamed"], size=n_cells),
            "cell_type_V2": rng.choice(["FB_A", "FB_B", "FB_C"], size=n_cells),
            "leiden": rng.integers(0, 12, size=n_cells).astype(str),
            "leiden_default": rng.integers(0, 8, size=n_cells).astype(str),
            "predicted_doublet": rng.choice(["False", "True"], size=n_cells, p=[0.97, 0.03]),
            "doublet_score": rng.random(size=n_cells),
            "cell_type_V3": rng.choice(["Fib_1", "Fib_2", "Fib_3", "Fib_4"], size=n_cells),
            "condition_v2": rng.choice(["Homeostasis", "Perturbed"], size=n_cells),
        }
    )

    var = pd.DataFrame(index=gene_names)
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.var_names = gene_names
    adata.obsm["X_umap"] = rng.normal(loc=0, scale=4, size=(n_cells, 2)).astype(np.float32)
    return adata


@st.cache_resource(show_spinner=False)
def load_adata(path: str) -> Tuple[ad.AnnData, bool, str]:
    try:
        if path and os.path.exists(path):
            return ad.read_h5ad(path), False, path
        msg = (
            f"H5AD file not found at '{path}'. Running in DEMO mode with synthetic data. "
            "Set H5AD_PATH to your real dataset path."
        )
        return build_mock_adata(), True, msg
    except Exception as exc:  # pragma: no cover
        return build_mock_adata(), True, f"Failed loading '{path}' ({exc}). Running in DEMO mode."


@st.cache_data(show_spinner=False)
def get_umap_df(adata: ad.AnnData) -> pd.DataFrame:
    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP coordinates were not found in adata.obsm['X_umap'].")

    coords = adata.obsm["X_umap"]
    if coords.shape[1] < 2:
        raise ValueError("UMAP embedding must have at least 2 columns.")

    return pd.DataFrame({"UMAP1": coords[:, 0], "UMAP2": coords[:, 1]}, index=adata.obs_names)


@st.cache_data(show_spinner=False)
def get_filter_options(adata: ad.AnnData, columns: Tuple[str, ...]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for col in columns:
        if col in adata.obs.columns:
            out[col] = sorted(adata.obs[col].astype(str).dropna().unique().tolist())
    return out


def apply_filters(adata: ad.AnnData, selected_filters: Dict[str, List[str]]) -> np.ndarray:
    mask = np.ones(adata.n_obs, dtype=bool)
    for col, chosen in selected_filters.items():
        if chosen and col in adata.obs.columns:
            mask &= adata.obs[col].astype(str).isin(chosen).to_numpy()
    return mask


def choose_plot_indices(
    filtered_indices: np.ndarray,
    view_mode: str,
    max_points: int,
    random_seed: int = 42,
) -> Tuple[np.ndarray, str]:
    n_filtered = filtered_indices.size
    if n_filtered == 0:
        return filtered_indices, "No cells matched the current filters."

    use_downsample = view_mode == "Downsampled"
    if view_mode == "Auto":
        use_downsample = n_filtered > max_points

    if use_downsample and n_filtered > max_points:
        rng = np.random.default_rng(random_seed)
        selected = np.sort(rng.choice(filtered_indices, size=max_points, replace=False))
        return selected, f"Displaying {max_points:,}/{n_filtered:,} filtered cells (downsampled)."

    return filtered_indices, f"Displaying all {n_filtered:,} filtered cells."


def get_expression_source_options(adata: ad.AnnData) -> List[str]:
    opts = ["X"]
    if adata.raw is not None:
        opts.append("raw")
    if adata.layers.keys():
        opts.extend([f"layer:{name}" for name in adata.layers.keys()])
    return opts


def get_expression_source(adata: ad.AnnData, source_key: str) -> Tuple[object, pd.Index, pd.DataFrame]:
    if source_key == "X":
        return adata.X, adata.var_names, adata.var
    if source_key == "raw":
        if adata.raw is None:
            raise ValueError("Expression source 'raw' is not available in this dataset.")
        return adata.raw.X, adata.raw.var_names, adata.raw.var

    if source_key.startswith("layer:"):
        layer_name = source_key.split(":", 1)[1]
        if layer_name not in adata.layers:
            raise ValueError(f"Layer '{layer_name}' is not available in adata.layers.")
        return adata.layers[layer_name], adata.var_names, adata.var

    raise ValueError(f"Unknown expression source '{source_key}'.")


def resolve_gene_index(
    query: str,
    gene_id_field: str,
    gene_names: pd.Index,
    gene_metadata: pd.DataFrame,
) -> Optional[int]:
    if not query:
        return None

    if gene_id_field == "var_names":
        if query in gene_names:
            return int(gene_names.get_loc(query))
        return None

    if gene_id_field not in gene_metadata.columns:
        return None

    matches = np.flatnonzero(gene_metadata[gene_id_field].astype(str).to_numpy() == query)
    if len(matches) == 0:
        return None
    return int(matches[0])


def extract_gene_for_indices(matrix: object, gene_idx: int, obs_indices: np.ndarray) -> np.ndarray:
    """Extract one gene for selected rows only. Never densifies full matrix."""
    values = matrix[obs_indices, gene_idx]
    if sparse.issparse(values):
        return values.toarray().ravel()
    return np.asarray(values).ravel()


def render_umap(df: pd.DataFrame, color: str, title: str, continuous: bool = False) -> None:
    fig = px.scatter(
        df,
        x="UMAP1",
        y="UMAP2",
        color=color,
        title=title,
        render_mode="webgl",
        opacity=0.75,
        color_continuous_scale="Viridis" if continuous else None,
    )
    fig.update_traces(marker={"size": 3})
    fig.update_layout(height=620, legend={"itemsizing": "constant"})
    st.plotly_chart(fig, use_container_width=True)


def main() -> None:
    st.set_page_config(page_title="Single-cell AnnData Explorer", layout="wide")
    st.title("Single-cell AnnData Explorer")

    default_path = os.getenv("H5AD_PATH", DEFAULT_H5AD_PATH)

    with st.sidebar:
        st.header("Controls")
        h5ad_path = st.text_input("Dataset path", value=default_path)

    adata, is_demo, status_msg = load_adata(h5ad_path)
    if is_demo:
        st.warning(f"DEMO MODE: {status_msg}")
    else:
        st.success(f"Loaded dataset from: {status_msg}")

    if "X_umap" not in adata.obsm:
        st.error("Missing required UMAP embedding: adata.obsm['X_umap'].")
        st.stop()

    obs_cols = adata.obs.columns.astype(str).tolist()
    color_candidates = obs_cols if obs_cols else ["None"]
    filter_options = get_filter_options(adata, tuple(FILTER_COLUMNS))

    with st.sidebar:
        color_by = st.selectbox(
            "Metadata color selector",
            options=color_candidates,
            index=color_candidates.index("cell_type") if "cell_type" in color_candidates else 0,
        )

        st.subheader("Filters")
        if st.button("Reset filters", use_container_width=True):
            for col in filter_options:
                st.session_state[f"filter_{col}"] = []

        selected_filters: Dict[str, List[str]] = {}
        for col, options in filter_options.items():
            key = f"filter_{col}"
            if key not in st.session_state:
                st.session_state[key] = []
            selected_filters[col] = st.multiselect(col, options=options, key=key)

        st.subheader("Visualization performance")
        view_mode = st.radio("UMAP view mode", ["Auto", "Downsampled", "Full"], index=0)
        max_points = int(
            st.number_input(
                "Display point cap",
                min_value=10000,
                max_value=max(10000, adata.n_obs),
                value=min(DEFAULT_MAX_POINTS, max(10000, adata.n_obs)),
                step=10000,
                help="Used in Auto/Downsampled mode to keep plotting responsive.",
            )
        )

        st.subheader("Gene query")
        expression_source = st.selectbox("Expression source", options=get_expression_source_options(adata))
        gene_id_options = ["var_names"] + adata.var.columns.astype(str).tolist()
        gene_id_field = st.selectbox(
            "Gene identifier field",
            options=gene_id_options,
            help="Default is adata.var_names. Use a var column if you keep alternate gene IDs.",
        )
        gene_query = st.text_input(
            "Gene search box",
            value="",
            placeholder="e.g., COL1A1",
            help="Gene names are resolved from var_names by default.",
        ).strip()

    filter_mask = apply_filters(adata, selected_filters)
    filtered_indices = np.flatnonzero(filter_mask)
    n_selected = int(filtered_indices.size)

    plotted_indices, view_msg = choose_plot_indices(filtered_indices, view_mode, max_points)

    tabs = st.tabs(["Explorer", "Data dictionary"])
    with tabs[0]:
        st.info(view_msg)
        if view_mode != "Full":
            st.caption(
                "Tip: Full mode keeps all filtered cells and can be slow at this scale. "
                "Downsampled mode improves interaction speed."
            )

        umap_df = get_umap_df(adata)
        plot_df = umap_df.iloc[plotted_indices].copy()
        plot_df[color_by] = adata.obs.iloc[plotted_indices][color_by].astype(str).to_numpy()

        st.subheader("Dataset summary")
        c1, c2, c3 = st.columns(3)
        c1.metric("Cells (total)", f"{adata.n_obs:,}")
        c2.metric("Genes", f"{adata.n_vars:,}")
        c3.metric("Cells after filters", f"{n_selected:,}")
        st.caption("Available metadata columns: " + ", ".join(obs_cols))

        st.subheader("About this dataset")
        st.markdown(
            """
            This portal visualizes one AnnData (`.h5ad`) file with precomputed UMAP coordinates.
            It supports metadata coloring, metadata filters, and on-demand single-gene overlays.
            For large datasets, the app can display a capped random subset of cells for smoother plotting.
            """
        )

        st.subheader(f"UMAP colored by metadata: {color_by}")
        render_umap(plot_df, color=color_by, title=f"UMAP • {color_by}")

        st.subheader("UMAP colored by gene expression")
        if gene_query:
            try:
                expr_matrix, gene_names, gene_metadata = get_expression_source(adata, expression_source)
                gene_idx = resolve_gene_index(gene_query, gene_id_field, gene_names, gene_metadata)
                if gene_idx is None:
                    st.error(
                        f"Gene '{gene_query}' was not found using '{gene_id_field}' in source '{expression_source}'."
                    )
                else:
                    plot_df["gene_expression"] = extract_gene_for_indices(
                        expr_matrix, gene_idx, plotted_indices
                    )
                    render_umap(
                        plot_df,
                        color="gene_expression",
                        title=f"UMAP • {gene_query} expression ({expression_source})",
                        continuous=True,
                    )
            except ValueError as exc:
                st.error(str(exc))
        else:
            st.info("Select a gene in the sidebar to display expression on UMAP.")

    with tabs[1]:
        st.subheader("Data dictionary (obs)")
        data_dict = pd.DataFrame(
            {
                "column": adata.obs.columns.astype(str),
                "dtype": [str(dtype) for dtype in adata.obs.dtypes],
            }
        )
        st.dataframe(data_dict, use_container_width=True)


if __name__ == "__main__":
    main()
