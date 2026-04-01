import numpy as np
import pandas as pd

from core import apply_filters, choose_plot_indices, resolve_gene_index, stratified_sample_indices, validate_core_metadata


def test_apply_filters_combines_columns():
    obs = pd.DataFrame(
        {
            "condition": ["CTRL", "HCM", "CTRL", "AMI"],
            "cell_type": ["F1", "F1", "F2", "F2"],
        }
    )
    mask = apply_filters(obs, {"condition": ["CTRL"], "cell_type": ["F2"]})
    assert mask.tolist() == [False, False, True, False]


def test_stratified_sample_indices_preserves_total():
    filtered_indices = np.arange(100)
    strata = np.array(["A"] * 80 + ["B"] * 20)
    sampled = stratified_sample_indices(filtered_indices, strata, max_points=25, random_seed=1)
    assert sampled.size == 25
    sampled_groups = strata[sampled]
    assert (sampled_groups == "A").sum() >= 19
    assert (sampled_groups == "B").sum() >= 4


def test_choose_plot_indices_auto_downsamples_with_message():
    filtered = np.arange(200)
    picked, msg = choose_plot_indices(filtered, view_mode="Auto", max_points=50)
    assert picked.size == 50
    assert "downsampled" in msg.lower()


def test_resolve_gene_index_var_names_match():
    idx = resolve_gene_index("COL1A1", "var_names", pd.Index(["ACTB", "COL1A1"]), pd.DataFrame())
    assert idx == 1


def test_validate_core_metadata_reports_missing():
    missing = validate_core_metadata(pd.DataFrame({"cell_type": ["F1", "F2"]}))
    assert missing == ["condition"]
