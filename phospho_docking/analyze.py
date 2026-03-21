"""Selectivity analysis: compute ΔΔG between phospho and unphospho docking results."""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from phospho_docking.dock import DockingResult


@dataclass
class SelectivityResult:
    """Selectivity analysis for a single peptide pair."""

    name: str
    dg_phospho: float  # kcal/mol
    dg_unphospho: float  # kcal/mol
    ddg: float  # ΔΔG = phospho - unphospho (negative = prefers phospho)
    selective: bool  # True if |ΔΔG| > threshold


def compute_selectivity(
    phospho: DockingResult,
    unphospho: DockingResult,
    threshold: float = 0.5,
) -> SelectivityResult:
    """Compute selectivity (ΔΔG) for a phospho/unphospho pair.

    Args:
        phospho: Docking result for phosphorylated peptide.
        unphospho: Docking result for unphosphorylated peptide.
        threshold: Minimum |ΔΔG| (kcal/mol) to consider selective.

    Returns:
        SelectivityResult with ΔΔG and selectivity flag.
    """
    ddg = phospho.best_affinity - unphospho.best_affinity

    # Extract base name (remove _phospho/_unphospho suffix)
    name = phospho.name.replace("_phospho", "").replace("_unphospho", "")

    return SelectivityResult(
        name=name,
        dg_phospho=phospho.best_affinity,
        dg_unphospho=unphospho.best_affinity,
        ddg=ddg,
        selective=abs(ddg) > threshold,
    )


def analyze_all(
    results: list[tuple[DockingResult, DockingResult]],
    threshold: float = 0.5,
) -> pd.DataFrame:
    """Analyze selectivity across all peptide pairs.

    Args:
        results: List of (phospho_result, unphospho_result) tuples.
        threshold: Minimum |ΔΔG| for selectivity.

    Returns:
        DataFrame sorted by ΔΔG (most selective first).
    """
    selectivities = []
    for phospho, unphospho in results:
        sel = compute_selectivity(phospho, unphospho, threshold)
        selectivities.append(sel)

    df = pd.DataFrame([
        {
            "name": s.name,
            "dG_phospho": s.dg_phospho,
            "dG_unphospho": s.dg_unphospho,
            "ddG": s.ddg,
            "selective": s.selective,
        }
        for s in selectivities
    ])

    return df.sort_values("ddG", ascending=True)


def save_scores(df: pd.DataFrame, output_path: Path) -> Path:
    """Save selectivity scores to CSV.

    Args:
        df: DataFrame from analyze_all().
        output_path: Where to write the CSV.

    Returns:
        Path to the written CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    return output_path
