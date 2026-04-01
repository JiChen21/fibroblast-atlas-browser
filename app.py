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
STRICT_DATA_MODE = os.getenv("STRICT_DATA", "false").strip().lower() in {"1", "true", "yes", "on"}
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
CORE_REQUIRED_OBS = ["cell_type", "condition"]
CELL_TYPE_COLOR_MAP = {
    "F1_Basal": "#8dd3c7",
    "F2_PSL": "#ffffb3",
    "F3_Int": "#bebada",
    "F4_Act": "#fb8072",
    "F5_ECM": "#80b1d3",
    "F6_SLIT3+": "#fdb462",
    "F7_Inflam": "#b3de69",
    "F8_ITM2B+": "#fccde5",
    "F9_Mech": "#e08214",
    "F10_CML": "#ccebc5",
    "F11_ECL": "#bc80bd",
}
CONDITION_COLOR_MAP = {
    "CTRL": "#b8e186",
    "HCM": "#fdb462",
    "ACM": "#fccde5",
    "AS": "#fee090",
    "MI": "#d6604d",
    "ICM": "#1d91c0",
    "DCM": "#fb8072",
    "HF": "#c51b7d",
    "COVID19": "#17becf",
}


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


def get_umap_df(adata: ad.AnnData) -> pd.DataFrame:
    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP coordinates were not found in adata.obsm['X_umap'].")

    coords = adata.obsm["X_umap"]
    if coords.shape[1] < 2:
        raise ValueError("UMAP embedding must have at least 2 columns.")

    return pd.DataFrame({"UMAP1": coords[:, 0], "UMAP2": coords[:, 1]}, index=adata.obs_names)


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
    strata_values: Optional[np.ndarray] = None,
    random_seed: int = 42,
) -> Tuple[np.ndarray, str]:
    n_filtered = filtered_indices.size
    if n_filtered == 0:
        return filtered_indices, "No cells matched the current filters."

    use_downsample = view_mode == "Downsampled"
    if view_mode == "Auto":
        use_downsample = n_filtered > max_points

    if use_downsample and n_filtered > max_points:
        if strata_values is not None and strata_values.size == n_filtered:
            selected = stratified_sample_indices(
                filtered_indices=filtered_indices,
                strata_values=strata_values,
                max_points=max_points,
                random_seed=random_seed,
            )
            return selected, f"Displaying {max_points:,}/{n_filtered:,} filtered cells (stratified downsampled)."
        rng = np.random.default_rng(random_seed)
        selected = np.sort(rng.choice(filtered_indices, size=max_points, replace=False))
        return selected, f"Displaying {max_points:,}/{n_filtered:,} filtered cells (downsampled)."

    return filtered_indices, f"Displaying all {n_filtered:,} filtered cells."


def stratified_sample_indices(
    filtered_indices: np.ndarray,
    strata_values: np.ndarray,
    max_points: int,
    random_seed: int = 42,
) -> np.ndarray:
    """Proportionally downsample while preserving group composition when possible."""
    rng = np.random.default_rng(random_seed)
    strata_values = strata_values.astype(str)
    unique_groups, inverse = np.unique(strata_values, return_inverse=True)
    group_positions = [np.flatnonzero(inverse == i) for i in range(len(unique_groups))]
    group_sizes = np.array([len(pos) for pos in group_positions], dtype=int)
    n_total = int(group_sizes.sum())
    if n_total <= max_points:
        return filtered_indices

    exact_alloc = group_sizes / n_total * max_points
    alloc = np.floor(exact_alloc).astype(int)
    alloc = np.minimum(alloc, group_sizes)

    remainder = max_points - int(alloc.sum())
    if remainder > 0:
        frac = exact_alloc - np.floor(exact_alloc)
        order = np.argsort(-frac)
        for idx in order:
            if remainder == 0:
                break
            if alloc[idx] < group_sizes[idx]:
                alloc[idx] += 1
                remainder -= 1

    selected_parts: List[np.ndarray] = []
    for i, pos in enumerate(group_positions):
        take = int(alloc[i])
        if take <= 0:
            continue
        if take >= len(pos):
            picked_pos = pos
        else:
            picked_pos = rng.choice(pos, size=take, replace=False)
        selected_parts.append(filtered_indices[np.sort(picked_pos)])

    if not selected_parts:
        return np.sort(rng.choice(filtered_indices, size=max_points, replace=False))
    return np.sort(np.concatenate(selected_parts))


def validate_core_metadata(adata: ad.AnnData) -> List[str]:
    return [col for col in CORE_REQUIRED_OBS if col not in adata.obs.columns]


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


def get_discrete_color_map(column_name: str) -> Optional[Dict[str, str]]:
    if column_name == "cell_type":
        return CELL_TYPE_COLOR_MAP
    if column_name == "condition":
        return CONDITION_COLOR_MAP
    return None


def render_umap(df: pd.DataFrame, color: str, title: str, continuous: bool = False) -> None:
    discrete_color_map = None if continuous else get_discrete_color_map(color)
    fig = px.scatter(
        df,
        x="UMAP1",
        y="UMAP2",
        color=color,
        title=title,
        render_mode="webgl",
        opacity=0.75,
        color_continuous_scale="Viridis" if continuous else None,
        color_discrete_map=discrete_color_map,
    )
    fig.update_traces(marker={"size": 3})
    fig.update_layout(
        height=820,
        width=820,
        legend={"itemsizing": "constant"},
    )
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    st.plotly_chart(fig, width="content")


def render_violin(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    color: Optional[str] = None,
) -> None:
    color_discrete_map = get_discrete_color_map(color) if color else None
    fig = px.violin(
        df,
        x=x,
        y=y,
        color=color,
        color_discrete_map=color_discrete_map,
        box=True,
        points=False,
        title=title,
    )
    fig.update_layout(height=520)
    st.plotly_chart(fig, width="stretch")


def build_group_expression_stats(df: pd.DataFrame, group_col: str, expr_col: str) -> pd.DataFrame:
    grouped = (
        df.groupby(group_col, observed=False)[expr_col]
        .agg(
            mean_expression="mean",
            median_expression="median",
            pct_expressing=lambda s: float((s > 0).mean() * 100.0),
            n_cells="size",
        )
        .reset_index()
        .sort_values("mean_expression", ascending=False)
    )
    return grouped


def render_dotplot(
    df: pd.DataFrame,
    x: str,
    y: str,
    size: str,
    color: str,
    title: str,
) -> None:
    color_discrete_map = None
    color_continuous_scale = None
    if color in {"cell_type", "condition"}:
        color_discrete_map = get_discrete_color_map(color)
    else:
        color_continuous_scale = "Viridis"

    fig = px.scatter(
        df,
        x=x,
        y=y,
        size=size,
        color=color,
        color_discrete_map=color_discrete_map,
        color_continuous_scale=color_continuous_scale,
        title=title,
        hover_data={"mean_expression": ":.4f", "median_expression": ":.4f", "pct_expressing": ":.2f", "n_cells": True},
    )
    fig.update_layout(height=460)
    st.plotly_chart(fig, width="stretch")


def main() -> None:
    st.set_page_config(page_title="Single-cell AnnData Explorer", layout="wide")
    st.title("Cross-disease Human Heart Fibroblast Atlas")

    default_path = os.getenv("H5AD_PATH", DEFAULT_H5AD_PATH)
    adata, is_demo, status_msg = load_adata(default_path)
    if STRICT_DATA_MODE and is_demo:
        st.error(
            "STRICT_DATA is enabled. Failed to load a valid dataset from H5AD_PATH; "
            "demo fallback is disabled."
        )
        st.caption(status_msg)
        st.stop()

    if is_demo:
        st.warning(f"DEMO MODE: {status_msg}")
    else:
        st.success(f"Loaded dataset from: {status_msg}")

    if "X_umap" not in adata.obsm:
        st.error("Missing required UMAP embedding: adata.obsm['X_umap'].")
        st.stop()

    missing_core_obs = validate_core_metadata(adata)
    if missing_core_obs:
        st.error(
            "Dataset is missing required metadata columns for this app: "
            + ", ".join(missing_core_obs)
            + "."
        )
        st.caption("Please provide these columns in adata.obs before using this portal.")
        st.stop()

    obs_cols = adata.obs.columns.astype(str).tolist()
    color_candidates = obs_cols if obs_cols else ["None"]
    filter_options = get_filter_options(adata, tuple(FILTER_COLUMNS))

    with st.sidebar:
        st.header("Atlas controls")
        st.caption(f"Dataset source: `{default_path}`")
        color_by = st.selectbox(
            "Metadata color selector",
            options=color_candidates,
            index=color_candidates.index("cell_type") if "cell_type" in color_candidates else 0,
        )

        st.subheader("Filters")
        if st.button("Reset filters", width="stretch"):
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
        gene_query = st.text_input(
            "Gene search box",
            value="",
            placeholder="e.g., COL1A1",
            help="Gene names are resolved directly from adata.var_names.",
        ).strip()

    filter_mask = apply_filters(adata, selected_filters)
    filtered_indices = np.flatnonzero(filter_mask)
    n_selected = int(filtered_indices.size)

    strata_values = None
    if filtered_indices.size > 0:
        if "cell_type" in adata.obs.columns:
            strata_values = adata.obs.iloc[filtered_indices]["cell_type"].astype(str).to_numpy()
        elif "condition" in adata.obs.columns:
            strata_values = adata.obs.iloc[filtered_indices]["condition"].astype(str).to_numpy()

    plotted_indices, view_msg = choose_plot_indices(filtered_indices, view_mode, max_points, strata_values=strata_values)

    tabs = st.tabs(["Home", "Metadata Explorer", "Gene Query Module", "Data Dictionary"])
    with tabs[0]:
        st.subheader("Atlas overview")
        st.markdown(
            """
            Cardiac fibroblasts are central regulators of extracellular matrix remodeling and fibrosis across diverse
            heart diseases, yet their heterogeneity and shared disease-associated states remain incompletely defined.

            This atlas curates and integrates 21 publicly available human heart single-cell and single-nucleus
            transcriptomic datasets, comprising approximately 730,000 fibroblasts from control hearts and eight disease
            conditions: hypertrophic cardiomyopathy (HCM), arrhythmogenic cardiomyopathy (ACM), aortic stenosis (AS),
            myocardial infarction (MI), ischemic cardiomyopathy (ICM), dilated cardiomyopathy (DCM), heart failure (HF),
            and COVID-19-associated cardiac injury.

            The portal is designed as an interactive resource for exploring fibroblast subtypes, gene expression
            patterns, and cross-disease fibroblast remodeling in the human heart.
            """
        )
        st.subheader("How to use this portal")
        st.markdown(
            """
            1. Use **Metadata Explorer** to inspect UMAP distributions and subset cells with sidebar filters.  
            2. Use **Gene Query Module** to visualize gene expression on UMAP and compare cell_type/condition-level patterns.  
            3. Use **Data Dictionary** to inspect available metadata columns and dtypes.
            """
        )
        st.info(view_msg)

    with tabs[1]:
        st.subheader("Metadata Explorer")
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

        st.subheader(f"UMAP colored by metadata: {color_by}")
        render_umap(plot_df, color=color_by, title=f"UMAP • {color_by}")

    with tabs[2]:
        st.subheader("Gene Query Module")
        st.caption("Includes: expression UMAP overlay, cell_type violin plot, and condition-level dot/violin plots.")
        st.subheader("UMAP colored by gene expression")
        if gene_query:
            try:
                expr_matrix, gene_names, gene_metadata = get_expression_source(adata, expression_source)
                gene_idx = resolve_gene_index(gene_query, "var_names", gene_names, gene_metadata)
                if gene_idx is None:
                    st.error(
                        f"Gene '{gene_query}' was not found in var_names for source '{expression_source}'."
                    )
                else:
                    umap_df = get_umap_df(adata)
                    plot_df = umap_df.iloc[plotted_indices].copy()
                    plot_df["gene_expression"] = extract_gene_for_indices(
                        expr_matrix, gene_idx, plotted_indices
                    )
                    render_umap(
                        plot_df,
                        color="gene_expression",
                        title=f"UMAP • {gene_query} expression ({expression_source})",
                        continuous=True,
                    )

                    analysis_indices = filtered_indices
                    analysis_df = adata.obs.iloc[analysis_indices][["cell_type", "condition"]].copy()
                    analysis_df["cell_type"] = analysis_df["cell_type"].astype(str)
                    analysis_df["condition"] = analysis_df["condition"].astype(str)
                    analysis_df["gene_expression"] = extract_gene_for_indices(
                        expr_matrix, gene_idx, analysis_indices
                    )

                    st.subheader(f"{gene_query} expression by cell_type (violin)")
                    render_violin(
                        analysis_df,
                        x="cell_type",
                        y="gene_expression",
                        color="cell_type",
                        title=f"{gene_query} expression across cell types",
                    )

                    st.subheader(f"{gene_query} expression in different conditions")
                    c1, c2 = st.columns(2)
                    with c1:
                        cond_stats = build_group_expression_stats(
                            analysis_df, group_col="condition", expr_col="gene_expression"
                        )
                        render_dotplot(
                            cond_stats,
                            x="condition",
                            y="mean_expression",
                            size="pct_expressing",
                            color="condition",
                            title=f"{gene_query} condition dot plot (size=% expressing)",
                        )
                    with c2:
                        render_violin(
                            analysis_df,
                            x="condition",
                            y="gene_expression",
                            color="condition",
                            title=f"{gene_query} condition violin plot",
                        )
            except ValueError as exc:
                st.error(str(exc))
            except KeyError as exc:
                st.error(
                    f"Dataset is missing required metadata column: {exc}. "
                    "Please ensure obs includes both 'cell_type' and 'condition'."
                )
        else:
            st.info("Select a gene in the sidebar to display expression on UMAP.")

    with tabs[3]:
        st.subheader("Data dictionary (obs)")
        data_dict = pd.DataFrame(
            {
                "column": adata.obs.columns.astype(str),
                "dtype": [str(dtype) for dtype in adata.obs.dtypes],
            }
        )
        st.dataframe(data_dict, width="stretch")


if __name__ == "__main__":
    main()
