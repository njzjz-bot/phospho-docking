"""CLI entry point: orchestrates the full docking pipeline."""

import argparse
from pathlib import Path

from phospho_docking.analyze import analyze_all, save_scores
from phospho_docking.dock import dock_pair
from phospho_docking.generate_peptides import generate_pair, generate_pairs_from_csv
from phospho_docking.prepare_receptor import prepare
from phospho_docking.report import generate_report


def run_pipeline(
    receptor_pdb: Path,
    peptides_csv: Path,
    output_dir: Path,
    exhaustiveness: int = 32,
    n_poses: int = 5,
    padding: float = 10.0,
    threshold: float = 0.5,
) -> Path:
    """Run the full phosphopeptide docking pipeline.

    Args:
        receptor_pdb: Path to the receptor PDB file.
        peptides_csv: Path to the peptides CSV file.
        output_dir: Directory for all output files.
        exhaustiveness: Vina exhaustiveness parameter.
        n_poses: Number of binding poses per ligand.
        padding: Bounding box padding in Angstroms.
        threshold: ΔΔG threshold for selectivity.

    Returns:
        Path to the HTML report.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Step 1: Prepare receptor
    print("Preparing receptor...")
    receptor_pdbqt, box = prepare(receptor_pdb, output_dir, padding=padding)
    print(f"  Receptor PDBQT: {receptor_pdbqt}")
    print(f"  Box center: ({box.center_x:.1f}, {box.center_y:.1f}, {box.center_z:.1f})")
    print(f"  Box size: ({box.size_x:.1f}, {box.size_y:.1f}, {box.size_z:.1f})")

    # Step 2: Generate peptide pairs
    print("Generating peptide pairs...")
    ligand_dir = output_dir / "ligands"
    pairs = generate_pairs_from_csv(peptides_csv, ligand_dir)
    print(f"  Generated {len(pairs)} peptide pairs")

    # Step 3: Dock all pairs
    print("Docking...")
    poses_dir = output_dir / "poses"
    all_results = []
    for pair in pairs:
        print(f"  Docking {pair.name}...")
        phospho_result, unphospho_result = dock_pair(
            pair=pair,
            receptor_pdbqt=receptor_pdbqt,
            box=box,
            output_dir=poses_dir,
            exhaustiveness=exhaustiveness,
            n_poses=n_poses,
        )
        print(f"    Phospho:   {phospho_result.best_affinity:.2f} kcal/mol")
        print(f"    Unphospho: {unphospho_result.best_affinity:.2f} kcal/mol")
        all_results.append((phospho_result, unphospho_result))

    # Step 4: Analyze selectivity
    print("Analyzing selectivity...")
    df = analyze_all(all_results, threshold=threshold)
    scores_path = save_scores(df, output_dir / "scores.csv")
    print(f"  Scores saved to: {scores_path}")

    # Step 5: Generate report
    print("Generating report...")
    report_path = generate_report(
        df=df,
        output_dir=output_dir,
        receptor_pdbqt=receptor_pdbqt,
        threshold=threshold,
    )
    print(f"  Report: {report_path}")

    return report_path


def run_single_test(
    receptor_pdb: Path,
    sequence: str,
    site: int,
    residue_type: str,
    output_dir: Path,
    exhaustiveness: int = 32,
) -> None:
    """Quick test with a single peptide pair.

    Args:
        receptor_pdb: Path to receptor PDB.
        sequence: Peptide sequence.
        site: 1-indexed phosphorylation site.
        residue_type: Thr, Ser, or Tyr.
        output_dir: Output directory.
        exhaustiveness: Vina exhaustiveness.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Preparing receptor...")
    receptor_pdbqt, box = prepare(receptor_pdb, output_dir)

    print(f"Generating peptide pair: {sequence} (site {site}, {residue_type})...")
    pair = generate_pair(
        name="test_peptide",
        sequence=sequence,
        site=site,
        residue_type=residue_type,
        output_dir=output_dir / "ligands",
    )

    print("Docking...")
    phospho_result, unphospho_result = dock_pair(
        pair=pair,
        receptor_pdbqt=receptor_pdbqt,
        box=box,
        output_dir=output_dir / "poses",
        exhaustiveness=exhaustiveness,
    )

    ddg = phospho_result.best_affinity - unphospho_result.best_affinity
    print(f"\nResults:")
    print(f"  Phospho:   {phospho_result.best_affinity:.2f} kcal/mol")
    print(f"  Unphospho: {unphospho_result.best_affinity:.2f} kcal/mol")
    print(f"  ΔΔG:       {ddg:.2f} kcal/mol")
    if ddg < -0.5:
        print(f"  → Construct PREFERS phosphorylated peptide (selective)")
    elif ddg > 0.5:
        print(f"  → Construct PREFERS unphosphorylated peptide")
    else:
        print(f"  → No significant selectivity detected")


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Phosphopeptide docking pipeline using AutoDock Vina",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Run command
    run_parser = subparsers.add_parser("run", help="Run full docking pipeline")
    run_parser.add_argument("--receptor", type=Path, required=True, help="Receptor PDB file")
    run_parser.add_argument("--peptides", type=Path, required=True, help="Peptides CSV file")
    run_parser.add_argument("--output", type=Path, default=Path("results"), help="Output directory")
    run_parser.add_argument("--exhaustiveness", type=int, default=32)
    run_parser.add_argument("--n-poses", type=int, default=5)
    run_parser.add_argument("--padding", type=float, default=10.0)
    run_parser.add_argument("--threshold", type=float, default=0.5)

    # Test command
    test_parser = subparsers.add_parser("test", help="Quick test with single peptide")
    test_parser.add_argument("--receptor", type=Path, required=True, help="Receptor PDB file")
    test_parser.add_argument("--sequence", type=str, required=True, help="Peptide sequence")
    test_parser.add_argument("--site", type=int, required=True, help="Phosphorylation site (1-indexed)")
    test_parser.add_argument("--residue", type=str, default="Thr", choices=["Thr", "Ser", "Tyr"])
    test_parser.add_argument("--output", type=Path, default=Path("results/test"))
    test_parser.add_argument("--exhaustiveness", type=int, default=32)

    args = parser.parse_args()

    if args.command == "run":
        run_pipeline(
            receptor_pdb=args.receptor,
            peptides_csv=args.peptides,
            output_dir=args.output,
            exhaustiveness=args.exhaustiveness,
            n_poses=args.n_poses,
            padding=args.padding,
            threshold=args.threshold,
        )
    elif args.command == "test":
        run_single_test(
            receptor_pdb=args.receptor,
            sequence=args.sequence,
            site=args.site,
            residue_type=args.residue,
            output_dir=args.output,
            exhaustiveness=args.exhaustiveness,
        )
