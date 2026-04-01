from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

CORE_REQUIRED_OBS = ["cell_type", "condition"]


def apply_filters(obs: pd.DataFrame, selected_filters: Dict[str, List[str]]) -> np.ndarray:
    mask = np.ones(obs.shape[0], dtype=bool)
    for col, chosen in selected_filters.items():
        if chosen and col in obs.columns:
            mask &= obs[col].astype(str).isin(chosen).to_numpy()
    return mask


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


def validate_core_metadata(obs: pd.DataFrame) -> List[str]:
    return [col for col in CORE_REQUIRED_OBS if col not in obs.columns]


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
