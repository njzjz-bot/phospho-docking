"""Integration test: 14-3-3ζ (PDB 1QJB) with Raf-1 pSer259 peptide.

14-3-3 proteins are canonical phosphopeptide readers. The phospho peptide
should dock with a more negative (better) affinity than the unphospho form.

Usage:
    uv run python tests/test_14_3_3.py
"""

from pathlib import Path

from phospho_docking.analyze import compute_selectivity
from phospho_docking.dock import dock_pair
from phospho_docking.generate_peptides import PeptidePair, generate_pair
from phospho_docking.prepare_receptor import DockingBox, prepare
from phospho_docking.report import generate_report

import pandas as pd


DATA_DIR = Path(__file__).parent.parent / "data"
RESULTS_DIR = Path(__file__).parent.parent / "results" / "test_14_3_3"


def main() -> None:
    # Step 1: Receptor prep
    pdb_path = DATA_DIR / "receptors" / "14-3-3_zeta.pdb"
    if not pdb_path.exists():
        print("Downloading 14-3-3ζ structure (PDB 1QJB)...")
        import urllib.request
        pdb_path.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(
            "https://files.rcsb.org/download/1QJB.pdb",
            pdb_path,
        )

    print("Preparing receptor...")
    receptor_pdbqt, _ = prepare(pdb_path, RESULTS_DIR)

    # Focused box around the chain A binding groove
    # (derived from co-crystallized peptide chain Q coordinates)
    box = DockingBox(
        center_x=23.1, center_y=-9.5, center_z=41.6,
        size_x=36.0, size_y=27.0, size_z=26.0,
    )

    # Step 2: Generate peptide pair (Raf-1 around Ser259)
    print("Generating peptide pair (Raf-1 Ser259)...")
    pair = generate_pair(
        name="Raf1_S259",
        sequence="QRSTSTP",
        site=3,
        residue_type="Ser",
        output_dir=RESULTS_DIR / "ligands",
    )

    # Step 3: Dock (exhaustiveness=8 for quick test, use 32 for production)
    print("Docking (exhaustiveness=8, ~2-5 min)...")
    phospho_result, unphospho_result = dock_pair(
        pair=pair,
        receptor_pdbqt=receptor_pdbqt,
        box=box,
        output_dir=RESULTS_DIR / "poses",
        exhaustiveness=8,
        n_poses=5,
    )

    # Step 4: Analyse
    sel = compute_selectivity(phospho_result, unphospho_result)
    print(f"\nResults:")
    print(f"  Phospho:   {sel.dg_phospho:.2f} kcal/mol")
    print(f"  Unphospho: {sel.dg_unphospho:.2f} kcal/mol")
    print(f"  ΔΔG:       {sel.ddg:.2f} kcal/mol")
    print(f"  Selective:  {sel.selective}")

    if sel.ddg < -0.5:
        print("  → PASS: 14-3-3 prefers phosphorylated peptide")
    elif sel.ddg > 0.5:
        print("  → UNEXPECTED: prefers unphosphorylated peptide")
    else:
        print("  → MARGINAL: no strong selectivity detected")

    # Step 5: Generate report
    df = pd.DataFrame([{
        "name": sel.name,
        "dG_phospho": sel.dg_phospho,
        "dG_unphospho": sel.dg_unphospho,
        "ddG": sel.ddg,
        "selective": sel.selective,
    }])
    report = generate_report(df, RESULTS_DIR, receptor_pdbqt)
    print(f"\nReport: {report}")


if __name__ == "__main__":
    main()
