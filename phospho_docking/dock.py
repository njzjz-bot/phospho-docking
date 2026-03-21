"""AutoDock Vina docking engine."""

from dataclasses import dataclass
from pathlib import Path

from vina import Vina

from phospho_docking.generate_peptides import PeptidePair
from phospho_docking.prepare_receptor import DockingBox


@dataclass
class DockingResult:
    """Result of a single Vina docking run."""

    name: str
    ligand_type: str  # "phospho" or "unphospho"
    best_affinity: float  # kcal/mol (most negative = best)
    all_affinities: list[float]
    pose_path: Path


def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    box: DockingBox,
    output_path: Path,
    exhaustiveness: int = 32,
    n_poses: int = 5,
) -> DockingResult:
    """Run AutoDock Vina for a single ligand against a receptor.

    Args:
        receptor_pdbqt: Path to receptor PDBQT file.
        ligand_pdbqt: Path to ligand PDBQT file.
        box: Docking search box definition.
        output_path: Where to save the output poses PDBQT.
        exhaustiveness: Search thoroughness (higher = slower, more accurate).
        n_poses: Number of binding poses to generate.

    Returns:
        DockingResult with binding affinities and pose file path.
    """
    v = Vina(sf_name="vina")
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(
        center=[box.center_x, box.center_y, box.center_z],
        box_size=[box.size_x, box.size_y, box.size_z],
    )

    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    v.write_poses(str(output_path), n_poses=n_poses, overwrite=True)

    energies = v.energies(n_poses=n_poses)
    affinities = [float(e[0]) for e in energies]  # First column is affinity

    return DockingResult(
        name=ligand_pdbqt.stem,
        ligand_type="phospho" if "phospho" in ligand_pdbqt.stem else "unphospho",
        best_affinity=affinities[0],
        all_affinities=affinities,
        pose_path=output_path,
    )


def dock_pair(
    pair: PeptidePair,
    receptor_pdbqt: Path,
    box: DockingBox,
    output_dir: Path,
    exhaustiveness: int = 32,
    n_poses: int = 5,
) -> tuple[DockingResult, DockingResult]:
    """Dock both phospho and unphospho peptides from a pair.

    Args:
        pair: PeptidePair with paths to both ligand PDBQT files.
        receptor_pdbqt: Path to receptor PDBQT.
        box: Docking search box.
        output_dir: Directory for output pose files.
        exhaustiveness: Vina exhaustiveness parameter.
        n_poses: Number of poses to generate per ligand.

    Returns:
        Tuple of (phospho_result, unphospho_result).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    phospho_result = run_vina(
        receptor_pdbqt=receptor_pdbqt,
        ligand_pdbqt=pair.phospho_pdbqt,
        box=box,
        output_path=output_dir / f"{pair.name}_phospho_poses.pdbqt",
        exhaustiveness=exhaustiveness,
        n_poses=n_poses,
    )

    unphospho_result = run_vina(
        receptor_pdbqt=receptor_pdbqt,
        ligand_pdbqt=pair.unphospho_pdbqt,
        box=box,
        output_path=output_dir / f"{pair.name}_unphospho_poses.pdbqt",
        exhaustiveness=exhaustiveness,
        n_poses=n_poses,
    )

    return phospho_result, unphospho_result
