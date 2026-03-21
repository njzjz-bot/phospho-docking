"""Report generation: HTML report with figures and PyMOL scripts."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from jinja2 import Template


REPORT_TEMPLATE = """\
<!DOCTYPE html>
<html>
<head>
    <title>Phosphopeptide Docking Report</title>
    <style>
        body { font-family: Arial, sans-serif; max-width: 1000px; margin: 0 auto; padding: 20px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; border-bottom: 2px solid #ecf0f1; padding-bottom: 5px; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px 12px; text-align: left; }
        th { background-color: #2c3e50; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .selective { color: #27ae60; font-weight: bold; }
        .not-selective { color: #e74c3c; }
        img { max-width: 100%; margin: 10px 0; }
        .summary { background: #eaf2f8; padding: 15px; border-radius: 5px; margin: 20px 0; }
    </style>
</head>
<body>
    <h1>Phosphopeptide Docking Report</h1>

    <div class="summary">
        <h2>Summary</h2>
        <p><strong>Peptide pairs tested:</strong> {{ n_pairs }}</p>
        <p><strong>Selective pairs (|ΔΔG| > {{ threshold }} kcal/mol):</strong> {{ n_selective }}</p>
        <p><strong>Best selectivity:</strong> {{ best_name }} (ΔΔG = {{ best_ddg }} kcal/mol)</p>
    </div>

    <h2>Selectivity Scores</h2>
    <table>
        <tr>
            <th>Peptide</th>
            <th>ΔG phospho (kcal/mol)</th>
            <th>ΔG unphospho (kcal/mol)</th>
            <th>ΔΔG (kcal/mol)</th>
            <th>Selective?</th>
        </tr>
        {% for row in rows %}
        <tr>
            <td>{{ row.name }}</td>
            <td>{{ "%.2f"|format(row.dG_phospho) }}</td>
            <td>{{ "%.2f"|format(row.dG_unphospho) }}</td>
            <td>{{ "%.2f"|format(row.ddG) }}</td>
            <td class="{{ 'selective' if row.selective else 'not-selective' }}">
                {{ "Yes" if row.selective else "No" }}
            </td>
        </tr>
        {% endfor %}
    </table>

    <h2>Binding Affinity Comparison</h2>
    <img src="figures/affinity_comparison.png" alt="Affinity comparison">

    <h2>Selectivity (ΔΔG)</h2>
    <img src="figures/selectivity.png" alt="Selectivity chart">

    <h2>Notes</h2>
    <ul>
        <li>Negative ΔΔG indicates the construct preferentially binds the phosphorylated peptide.</li>
        <li>Binding poses are saved as PDBQT files in <code>poses/</code>.</li>
        <li>PyMOL scripts for visualisation are in <code>pymol/</code>.</li>
    </ul>
</body>
</html>
"""


def plot_affinity_comparison(df: pd.DataFrame, output_path: Path) -> Path:
    """Create a grouped bar chart comparing phospho vs unphospho binding energies.

    Args:
        df: Scores DataFrame with dG_phospho and dG_unphospho columns.
        output_path: Where to save the PNG.

    Returns:
        Path to saved figure.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    x = np.arange(len(df))
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(8, len(df) * 1.5), 5))
    ax.bar(x - width / 2, df["dG_phospho"], width, label="Phospho", color="#2ecc71")
    ax.bar(x + width / 2, df["dG_unphospho"], width, label="Unphospho", color="#e74c3c")

    ax.set_ylabel("Binding affinity (kcal/mol)")
    ax.set_xlabel("Peptide")
    ax.set_title("Phospho vs Unphospho Binding Affinity")
    ax.set_xticks(x)
    ax.set_xticklabels(df["name"], rotation=45, ha="right")
    ax.legend()
    ax.axhline(y=0, color="black", linewidth=0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

    return output_path


def plot_selectivity(df: pd.DataFrame, output_path: Path, threshold: float = 0.5) -> Path:
    """Create a bar chart of ΔΔG values per peptide.

    Args:
        df: Scores DataFrame with ddG column.
        output_path: Where to save the PNG.
        threshold: Selectivity threshold line.

    Returns:
        Path to saved figure.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    colors = ["#2ecc71" if ddg < -threshold else "#e74c3c" if ddg > threshold else "#95a5a6"
              for ddg in df["ddG"]]

    fig, ax = plt.subplots(figsize=(max(8, len(df) * 1.5), 5))
    ax.bar(df["name"], df["ddG"], color=colors)
    ax.axhline(y=-threshold, color="green", linestyle="--", alpha=0.5, label=f"Threshold (±{threshold})")
    ax.axhline(y=threshold, color="red", linestyle="--", alpha=0.5)
    ax.axhline(y=0, color="black", linewidth=0.5)

    ax.set_ylabel("ΔΔG (kcal/mol)")
    ax.set_xlabel("Peptide")
    ax.set_title("Phospho Selectivity (ΔΔG)")
    ax.legend()
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

    return output_path


def generate_pymol_script(
    receptor_pdbqt: Path,
    pose_pdbqt: Path,
    output_path: Path,
    peptide_name: str,
) -> Path:
    """Generate a PyMOL .pml script for visualising a binding pose.

    Args:
        receptor_pdbqt: Path to receptor PDBQT.
        pose_pdbqt: Path to docked pose PDBQT.
        output_path: Where to write the .pml script.
        peptide_name: Label for the peptide.

    Returns:
        Path to the .pml script.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    script = f"""\
# PyMOL visualisation: {peptide_name}
load {receptor_pdbqt}, receptor
load {pose_pdbqt}, {peptide_name}

# Style
hide everything
show cartoon, receptor
show sticks, {peptide_name}
color marine, receptor
color green, {peptide_name}

# Highlight phosphate
select phosphate, {peptide_name} and name P+O1P+O2P+O3P
show spheres, phosphate
color red, phosphate

# Surface around binding site
select binding_site, receptor within 5 of {peptide_name}
show surface, binding_site
set transparency, 0.5, binding_site

zoom {peptide_name}, 5
"""
    output_path.write_text(script)
    return output_path


def generate_report(
    df: pd.DataFrame,
    output_dir: Path,
    receptor_pdbqt: Path | None = None,
    threshold: float = 0.5,
) -> Path:
    """Generate the full HTML report with figures.

    Args:
        df: Scores DataFrame from analyze_all().
        output_dir: Directory for all output files.
        receptor_pdbqt: Path to receptor PDBQT (for PyMOL scripts).
        threshold: Selectivity threshold.

    Returns:
        Path to the HTML report.
    """
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Generate figures
    plot_affinity_comparison(df, figures_dir / "affinity_comparison.png")
    plot_selectivity(df, figures_dir / "selectivity.png", threshold=threshold)

    # Generate PyMOL scripts
    if receptor_pdbqt:
        pymol_dir = output_dir / "pymol"
        pymol_dir.mkdir(parents=True, exist_ok=True)
        poses_dir = output_dir / "poses"
        for _, row in df.iterrows():
            for variant in ["phospho", "unphospho"]:
                pose_path = poses_dir / f"{row['name']}_{variant}_poses.pdbqt"
                if pose_path.exists():
                    generate_pymol_script(
                        receptor_pdbqt=receptor_pdbqt,
                        pose_pdbqt=pose_path,
                        output_path=pymol_dir / f"{row['name']}_{variant}.pml",
                        peptide_name=f"{row['name']}_{variant}",
                    )

    # Render HTML
    template = Template(REPORT_TEMPLATE)
    best_row = df.iloc[0] if len(df) > 0 else None
    html = template.render(
        n_pairs=len(df),
        n_selective=int(df["selective"].sum()),
        best_name=best_row["name"] if best_row is not None else "N/A",
        best_ddg=f"{best_row['ddG']:.2f}" if best_row is not None else "N/A",
        threshold=threshold,
        rows=df.to_dict("records"),
    )

    report_path = output_dir / "report.html"
    report_path.write_text(html)

    return report_path
