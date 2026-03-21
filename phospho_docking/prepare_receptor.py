"""Receptor preparation: PDB → PDBQT with bounding box for Vina."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from meeko import PDBQTWriterLegacy, MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass
class DockingBox:
    """Defines the 3D search space for Vina."""

    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    def to_dict(self) -> dict[str, float]:
        return {
            "center_x": self.center_x,
            "center_y": self.center_y,
            "center_z": self.center_z,
            "size_x": self.size_x,
            "size_y": self.size_y,
            "size_z": self.size_z,
        }


def compute_bounding_box(pdb_path: Path, padding: float = 10.0) -> DockingBox:
    """Compute a bounding box around the protein with padding (Angstroms).

    Args:
        pdb_path: Path to the receptor PDB file.
        padding: Extra space around the protein in each dimension.

    Returns:
        DockingBox centred on the protein with specified padding.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(pdb_path))

    coords = []
    for atom in structure.get_atoms():
        coords.append(atom.get_vector().get_array())
    coords = np.array(coords)

    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    center = (min_coords + max_coords) / 2.0
    size = (max_coords - min_coords) + 2 * padding

    return DockingBox(
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=float(size[0]),
        size_y=float(size[1]),
        size_z=float(size[2]),
    )


def prepare_receptor_pdbqt(pdb_path: Path, output_path: Path) -> Path:
    """Convert a receptor PDB to PDBQT format for AutoDock Vina.

    Uses Open Babel for robust PDB → PDBQT conversion (adds hydrogens,
    assigns Gasteiger charges, writes PDBQT atom types).

    Args:
        pdb_path: Input PDB file path.
        output_path: Where to write the PDBQT file.

    Returns:
        Path to the generated PDBQT file.
    """
    from openbabel import openbabel

    obconv = openbabel.OBConversion()
    obconv.SetInAndOutFormats("pdb", "pdbqt")
    # Receptor-mode flags: no flexible residues, add hydrogens
    obconv.AddOption("r", openbabel.OBConversion.OUTOPTIONS)
    obconv.AddOption("h", openbabel.OBConversion.OUTOPTIONS)

    mol = openbabel.OBMol()
    obconv.ReadFile(mol, str(pdb_path))

    # Remove water molecules
    atoms_to_delete = []
    for atom in openbabel.OBMolAtomIter(mol):
        res = atom.GetResidue()
        if res and res.GetName().strip() in ("HOH", "WAT"):
            atoms_to_delete.append(atom)
    for atom in reversed(atoms_to_delete):
        mol.DeleteAtom(atom)

    # Add hydrogens and assign Gasteiger charges
    mol.AddHydrogens()
    charge_model = openbabel.OBChargeModel.FindType("gasteiger")
    if charge_model:
        charge_model.ComputeCharges(mol)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    obconv.WriteFile(mol, str(output_path))

    return output_path


def prepare(pdb_path: Path, output_dir: Path, padding: float = 10.0) -> tuple[Path, DockingBox]:
    """Full receptor preparation: PDB → PDBQT + bounding box.

    Args:
        pdb_path: Path to receptor PDB file.
        output_dir: Directory to write output files.
        padding: Bounding box padding in Angstroms.

    Returns:
        Tuple of (PDBQT path, DockingBox).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_path = output_dir / f"{pdb_path.stem}_receptor.pdbqt"

    pdbqt = prepare_receptor_pdbqt(pdb_path, pdbqt_path)
    box = compute_bounding_box(pdb_path, padding=padding)

    return pdbqt, box
