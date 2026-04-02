import logging
import os
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from scipy import sparse
from core import (
    apply_filters,
    build_condition_subtype_counts,
    choose_plot_indices,
    compute_roe,
    resolve_gene_index,
    roe_symbol,
    validate_core_metadata,
)

DEFAULT_H5AD_PATH = "./data/FBs_adata.h5ad"
DEFAULT_MAX_POINTS = int(os.getenv("UMAP_MAX_POINTS", "200000"))
STRICT_DATA_MODE = os.getenv("STRICT_DATA", "false").strip().lower() in {"1", "true", "yes", "on"}
FILTER_COLUMNS = [
    "condition",
    "region",
    "cell_type"
]
CELL_TYPE_ORDER = [
    "F1_Basal",
    "F2_PSL",
    "F3_Int",
    "F4_Act",
    "F5_ECM",
    "F6_SLIT3+",
    "F7_Inflam",
    "F8_ITM2B+",
    "F9_Mech",
    "F10_CML",
    "F11_ECL",
]
CONDITION_ORDER = ["CTRL", "HCM", "ACM", "AS", "AMI", "ICM", "DCM", "HF", "COVID19"]
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
    "AMI": "#d6604d",
    "ICM": "#1d91c0",
    "DCM": "#fb8072",
    "HF": "#c51b7d",
    "COVID19": "#17becf",
}

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(
        level=os.getenv("LOG_LEVEL", "INFO").upper(),
        format="%(asctime)s %(levelname)s %(name)s - %(message)s",
    )


def resolve_home_image_path() -> Optional[str]:
    configured = os.getenv("HOME_IMAGE_PATH", "").strip()
    candidates = [
        configured,
        "./assets/home_overview.png",
        "/data/chenji/fibroblast-atlas-browser-new/plot/home_page.jpg",
    ]
    for path in candidates:
        if path and os.path.exists(path):
            return path
    return None


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


@st.cache_resource(show_spinner="Loading AnnData dataset...")
def load_adata(path: str) -> Tuple[ad.AnnData, bool, str]:
    try:
        if path and os.path.exists(path):
            logger.info("Loading dataset from %s", path)
            return ad.read_h5ad(path), False, path
        msg = (
            f"H5AD file not found at '{path}'. Running in DEMO mode with synthetic data. "
            "Set H5AD_PATH to your real dataset path."
        )
        logger.warning(msg)
        return build_mock_adata(), True, msg
    except Exception as exc:  # pragma: no cover
        logger.exception("Failed to load dataset from %s", path)
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


def render_umap(
    df: pd.DataFrame,
    color: str,
    title: str,
    continuous: bool = False,
    height: int = 620,
) -> None:
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
        height=height,
        legend={"itemsizing": "constant"},
        template="plotly_white",
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        font={"color": "#111827"},
    )
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    st.plotly_chart(fig, width="stretch", theme=None)


def render_violin(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    color: Optional[str] = None,
    height: int = 420,
    show_legend: bool = True,
) -> None:
    color_discrete_map = get_discrete_color_map(color) if color else None
    fig = px.violin(
        df,
        x=x,
        y=y,
        color=color,
        color_discrete_map=color_discrete_map,
        category_orders={
            "cell_type": CELL_TYPE_ORDER,
            "condition": CONDITION_ORDER,
        },
        box=True,
        points=False,
        title=title,
    )
    fig.update_layout(
        height=height,
        showlegend=show_legend,
        template="plotly_white",
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        font={"color": "#111827"},
    )
    st.plotly_chart(fig, width="stretch", theme=None)


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
    height: int = 380,
    show_legend: bool = True,
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
        category_orders={
            "cell_type": CELL_TYPE_ORDER,
            "condition": CONDITION_ORDER,
        },
        title=title,
        hover_data={"mean_expression": ":.4f", "median_expression": ":.4f", "pct_expressing": ":.2f", "n_cells": True},
    )
    fig.update_layout(
        height=height,
        showlegend=show_legend,
        template="plotly_white",
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        font={"color": "#111827"},
    )
    st.plotly_chart(fig, width="stretch", theme=None)


def render_condition_stacked_bar(
    df: pd.DataFrame,
    condition_col: str,
    subtype_col: str,
    value_col: str,
    title: str,
    height: int = 430,
) -> None:
    fig = px.bar(
        df,
        x=condition_col,
        y=value_col,
        color=subtype_col,
        barmode="stack",
        category_orders={"condition": CONDITION_ORDER, "cell_type": CELL_TYPE_ORDER},
        color_discrete_map=CELL_TYPE_COLOR_MAP if subtype_col == "cell_type" else None,
        title=title,
        hover_data={value_col: ":.2f", "n_observed": True},
    )
    fig.update_yaxes(title="Proportion (%)")
    fig.update_layout(
        height=height,
        template="plotly_white",
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        font={"color": "#111827"},
    )
    st.plotly_chart(fig, width="stretch", theme=None)


def render_roe_heatmap(
    roe_df: pd.DataFrame,
    condition_col: str,
    subtype_col: str,
    title: str,
    condition_order: Optional[List[str]] = None,
    subtype_order: Optional[List[str]] = None,
    clip_val: float = 3.0,
    height: int = 560,
) -> None:
    heatmap_df = roe_df.pivot(index=subtype_col, columns=condition_col, values="Ro_e")
    if subtype_order:
        row_order = [v for v in subtype_order if v in heatmap_df.index]
        row_order.extend([v for v in heatmap_df.index.tolist() if v not in row_order])
        heatmap_df = heatmap_df.reindex(index=row_order)
    if condition_order:
        col_order = [v for v in condition_order if v in heatmap_df.columns]
        col_order.extend([v for v in heatmap_df.columns.tolist() if v not in col_order])
        heatmap_df = heatmap_df.reindex(columns=col_order)

    roe_matrix = heatmap_df.to_numpy(dtype=float)
    log2_mat = np.log2(roe_matrix)

    finite_vals = log2_mat[np.isfinite(log2_mat)]
    max_abs = float(np.max(np.abs(finite_vals))) if finite_vals.size > 0 else clip_val
    if (not np.isfinite(max_abs)) or max_abs <= 0:
        max_abs = clip_val

    log2_mat[np.isinf(log2_mat) & (roe_matrix == 0)] = -max_abs
    plot_mat = np.clip(log2_mat, -clip_val, clip_val)
    text_mat = np.vectorize(roe_symbol)(roe_matrix)

    fig = px.imshow(
        plot_mat,
        x=heatmap_df.columns.tolist(),
        y=heatmap_df.index.tolist(),
        color_continuous_scale=[(0.0, "#2166AC"), (0.5, "white"), (1.0, "#B2182B")],
        zmin=-clip_val,
        zmax=clip_val,
        labels={"color": "log2(Ro/e)", "x": condition_col, "y": subtype_col},
        title=title,
        aspect="auto",
    )
    fig.update_traces(text=text_mat, texttemplate="%{text}", textfont={"size": 11, "color": "black"})
    fig.update_layout(
        height=height,
        template="plotly_white",
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        font={"color": "#111827"},
        coloraxis_colorbar={
            "tickmode": "array",
            "tickvals": [-3, -2, -1, 0, 1, 2, 3],
            "ticktext": ["≤-3", "-2", "-1", "0", "1", "2", "≥3"],
        },
    )
    st.plotly_chart(fig, width="stretch", theme=None)


def apply_global_styles() -> None:
    st.markdown(
        """
        <style>
        :root {
            --background-color: #ffffff;
            --secondary-background-color: #f7f7f7;
            --text-color: #111827;
        }
        .stApp,
        [data-testid="stAppViewContainer"],
        [data-testid="stHeader"],
        [data-testid="stSidebar"] {
            background-color: #ffffff !important;
        }
        .stApp,
        .stApp p,
        .stApp li,
        .stApp label,
        .stApp span,
        .stApp div,
        .stApp h1,
        .stApp h2,
        .stApp h3,
        .stApp h4,
        .stApp h5,
        .stApp h6 {
            color: #111827 !important;
        }
        /* Streamlit/BaseWeb form controls */
        .stSelectbox [data-baseweb="select"] > div,
        .stMultiSelect [data-baseweb="select"] > div,
        .stTextInput input,
        .stNumberInput input,
        .stTextArea textarea,
        [data-baseweb="select"] div[role="combobox"] {
            background-color: #ffffff !important;
            color: #111827 !important;
        }
        .stButton > button,
        .stDownloadButton > button {
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
        }
        .stRadio > div[role="radiogroup"] label {
            color: #111827 !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def main() -> None:
    st.set_page_config(page_title="Single-cell AnnData Explorer", layout="wide")
    apply_global_styles()
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

    if "X_umap" not in adata.obsm:
        st.error("Missing required UMAP embedding: adata.obsm['X_umap'].")
        st.stop()

    missing_core_obs = validate_core_metadata(adata.obs)
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
    module = st.radio(
        "Navigation",
        ["Home", "Metadata Explorer", "Gene Query Module", "Condition / Disease Browser", "Data Dictionary"],
        horizontal=True,
        label_visibility="collapsed",
    )

    color_by = "cell_type" if "cell_type" in color_candidates else color_candidates[0]
    selected_filters: Dict[str, List[str]] = {}
    view_mode = "Auto"
    max_points = min(DEFAULT_MAX_POINTS, max(10000, adata.n_obs))
    expression_source = "X"
    gene_query = ""

    selected_conditions_for_browser: List[str] = []

    if module in {"Metadata Explorer", "Gene Query Module", "Condition / Disease Browser"}:
        with st.sidebar:
            st.header(f"{module} controls")
            color_by = st.selectbox(
                "Metadata color selector",
                options=color_candidates,
                index=color_candidates.index(color_by) if color_by in color_candidates else 0,
            )

            st.subheader("Filters")
            if st.button("Reset filters", width="stretch"):
                for col in filter_options:
                    st.session_state[f"filter_{col}"] = []

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

            if module == "Gene Query Module":
                st.subheader("Gene query")
                expression_source = st.selectbox("Expression source", options=get_expression_source_options(adata))
                gene_query = st.text_input(
                    "Gene search box",
                    value="",
                    placeholder="e.g., COL1A1",
                    help="Gene names are resolved directly from adata.var_names.",
                ).strip()

            if module == "Condition / Disease Browser":
                condition_options_raw = filter_options.get("condition", [])
                condition_options = [c for c in CONDITION_ORDER if c in condition_options_raw]
                condition_options.extend([c for c in condition_options_raw if c not in condition_options])
                selected_conditions_for_browser = st.multiselect(
                    "Choose conditions",
                    options=condition_options,
                    default=condition_options,
                    help="Select one or multiple conditions to compare subtype composition and enrichment.",
                )

    filter_mask = apply_filters(adata.obs, selected_filters)
    filtered_indices = np.flatnonzero(filter_mask)
    n_selected = int(filtered_indices.size)

    strata_values = None
    if filtered_indices.size > 0:
        if "cell_type" in adata.obs.columns:
            strata_values = adata.obs.iloc[filtered_indices]["cell_type"].astype(str).to_numpy()
        elif "condition" in adata.obs.columns:
            strata_values = adata.obs.iloc[filtered_indices]["condition"].astype(str).to_numpy()

    plotted_indices, view_msg = choose_plot_indices(filtered_indices, view_mode, max_points, strata_values=strata_values)
    logger.info(
        "View mode=%s filtered=%d plotted=%d max_points=%d module=%s",
        view_mode,
        n_selected,
        int(plotted_indices.size),
        max_points,
        module,
    )

    if module == "Home":
        st.subheader("Atlas overview")
        st.markdown(
            """
            Cardiac fibroblasts are central regulators of extracellular matrix remodeling and fibrosis across diverse
            heart diseases, yet their heterogeneity and shared disease-associated states remain incompletely defined.

            This atlas curates and integrates 21 publicly available human heart single-cell and single-nucleus
            transcriptomic datasets, comprising approximately 730,000 fibroblasts from control hearts and eight disease
            conditions: hypertrophic cardiomyopathy (HCM), arrhythmogenic cardiomyopathy (ACM), aortic stenosis (AS),
            acute myocardial infarction (AMI), ischemic cardiomyopathy (ICM), dilated cardiomyopathy (DCM), heart failure (HF),
            and COVID-19-associated cardiac injury.

            The portal is designed as an interactive resource for exploring fibroblast subtypes, gene expression
            patterns, and cross-disease fibroblast remodeling in the human heart.
            """
        )
        home_image_path = resolve_home_image_path()
        if home_image_path:
            logger.info("Using homepage image: %s", home_image_path)
            _, center_col, _ = st.columns([1, 6, 1])
            with center_col:
                st.image(
                    home_image_path,
                    caption="Fibroblast atlas overview figure",
                    use_container_width=True,
                )
        else:
            st.caption(
                "Tip: set `HOME_IMAGE_PATH` or place an image at "
                "`./assets/home_overview.png` "
                "(also supports `/data/chenji/fibroblast-atlas-browser-new/plot/home_page.jpg`)."
            )
        st.subheader("How to use this portal")
        st.markdown(
            """
            1. Use **Metadata Explorer** to inspect UMAP distributions and subset cells with sidebar filters.  
            2. Use **Gene Query Module** to visualize gene expression on UMAP and compare cell_type/condition-level patterns.  
            3. Use **Data Dictionary** to inspect available metadata columns and dtypes.
            """
        )
    elif module == "Metadata Explorer":
        st.subheader("Metadata Explorer")
        st.info(view_msg)
        if view_mode != "Full":
            st.caption(
                "Tip: Full mode keeps all filtered cells and can be slow at this scale. "
                "Downsampled mode improves interaction speed."
            )

        with st.spinner("Preparing metadata visualization..."):
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
        render_umap(plot_df, color=color_by, title=f"UMAP • {color_by}", height=640)

    elif module == "Gene Query Module":
        st.subheader("Gene Query Module")
        st.caption("Compact layout: expression UMAP + cell_type violin on top, condition dot/violin below.")
        if gene_query:
            try:
                with st.spinner(f"Querying gene '{gene_query}' and rendering plots..."):
                    expr_matrix, gene_names, gene_metadata = get_expression_source(adata, expression_source)
                    gene_idx = resolve_gene_index(gene_query, "var_names", gene_names, gene_metadata)
                    if gene_idx is None:
                        logger.info("Gene query miss source=%s query=%s", expression_source, gene_query)
                        st.error(
                            f"Gene '{gene_query}' was not found in var_names for source '{expression_source}'."
                        )
                    else:
                        umap_df = get_umap_df(adata)
                        plot_df = umap_df.iloc[plotted_indices].copy()
                        plot_df["gene_expression"] = extract_gene_for_indices(
                            expr_matrix, gene_idx, plotted_indices
                        )
                        analysis_indices = filtered_indices
                        analysis_df = adata.obs.iloc[analysis_indices][["cell_type", "condition"]].copy()
                        analysis_df["cell_type"] = analysis_df["cell_type"].astype(str)
                        analysis_df["condition"] = analysis_df["condition"].astype(str)
                        analysis_df["gene_expression"] = extract_gene_for_indices(
                            expr_matrix, gene_idx, analysis_indices
                        )

                        top_left, top_right = st.columns(2)
                        with top_left:
                            render_umap(
                                plot_df,
                                color="gene_expression",
                                title=f"UMAP • {gene_query} expression ({expression_source})",
                                continuous=True,
                                height=500,
                            )
                        with top_right:
                            render_violin(
                                analysis_df,
                                x="cell_type",
                                y="gene_expression",
                                color="cell_type",
                                title=f"{gene_query} expression across cell types",
                                height=500,
                                show_legend=False,
                            )

                        st.subheader(f"{gene_query} expression in different conditions")
                        bottom_left, bottom_right = st.columns(2)
                        with bottom_left:
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
                                height=360,
                                show_legend=False,
                            )
                        with bottom_right:
                            render_violin(
                                analysis_df,
                                x="condition",
                                y="gene_expression",
                                color="condition",
                                title=f"{gene_query} condition violin plot",
                                height=360,
                                show_legend=False,
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

    elif module == "Data Dictionary":
        st.subheader("Data dictionary (obs)")
        data_dict = pd.DataFrame(
            {
                "column": adata.obs.columns.astype(str),
                "dtype": [str(dtype) for dtype in adata.obs.dtypes],
            }
        )
        st.dataframe(data_dict, width="stretch")

    elif module == "Condition / Disease Browser":
        st.subheader("Condition / Disease Browser")
        st.caption("Compare subtype proportions across selected conditions and inspect log2(Observed/Expected) enrichment.")

        if "condition" not in adata.obs.columns or "cell_type" not in adata.obs.columns:
            st.error("This module requires both 'condition' and 'cell_type' columns in adata.obs.")
            st.stop()

        analysis_obs = adata.obs.iloc[filtered_indices][["condition", "cell_type"]].copy()
        analysis_obs["condition"] = analysis_obs["condition"].astype(str)
        analysis_obs["cell_type"] = analysis_obs["cell_type"].astype(str)

        if analysis_obs.empty:
            st.warning("No cells matched the current filters. Please adjust filters in the sidebar.")
            st.stop()

        selected_conditions = selected_conditions_for_browser or sorted(analysis_obs["condition"].unique().tolist())

        counts = build_condition_subtype_counts(
            analysis_obs,
            condition_col="condition",
            subtype_col="cell_type",
            selected_conditions=selected_conditions,
            condition_order=CONDITION_ORDER,
            subtype_order=CELL_TYPE_ORDER,
        )
        if counts.empty:
            st.warning("No data available for the selected conditions under current filters.")
            st.stop()

        proportion_df = counts.copy()
        totals = proportion_df.groupby("condition", observed=False)["n_observed"].transform("sum")
        proportion_df["proportion_pct"] = np.where(totals > 0, proportion_df["n_observed"] / totals * 100.0, 0.0)

        left_col, right_col = st.columns(2)
        with left_col:
            st.markdown("#### Subtype composition (stacked proportions)")
            render_condition_stacked_bar(
                proportion_df,
                condition_col="condition",
                subtype_col="cell_type",
                value_col="proportion_pct",
                title="Cell subtype proportions by condition",
            )

        roe_df = compute_roe(counts, condition_col="condition", subtype_col="cell_type")
        with right_col:
            st.markdown("#### Ro/e heatmap (log2 scale)")
            render_roe_heatmap(
                roe_df,
                condition_col="condition",
                subtype_col="cell_type",
                title="Fibroblast subtype enrichment (log2[Observed/Expected])",
                condition_order=CONDITION_ORDER,
                subtype_order=CELL_TYPE_ORDER,
            )

        st.caption(
            "Symbol rules: +++ / --- (|log2FC| ≥ 0.58), ++ / -- (≥ 0.32), + / - (≥ 0.10), +/- (near neutral). "
            "Heatmap color scale is clipped to [-3, 3] on log2(Ro/e)."
        )


if __name__ == "__main__":
    main()
