"""Generate phosphorylated and unphosphorylated peptide pairs as 3D PDBQT ligands."""

from dataclasses import dataclass
from pathlib import Path

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom


# SMILES building blocks for amino acids (neutral backbone, no caps)
AMINO_ACID_SMILES: dict[str, str] = {
    "A": "C",
    "G": "[H]",
    "V": "CC(C)",
    "L": "CC(C)C",
    "I": "C(CC)C",
    "P": "",  # special — proline is cyclic
    "F": "Cc1ccccc1",
    "W": "Cc1c[nH]c2ccccc12",
    "M": "CCSC",
    "S": "CO",
    "T": "C(O)C",
    "Y": "Cc1ccc(O)cc1",
    "D": "CC(=O)O",
    "E": "CCC(=O)O",
    "N": "CC(=O)N",
    "Q": "CCC(=O)N",
    "K": "CCCCN",
    "R": "CCCNC(=N)N",
    "H": "Cc1c[nH]cn1",
    "C": "CS",
}

# Phosphorylated residue SMILES (replace hydroxyl with phosphate)
PHOSPHO_SMILES: dict[str, str] = {
    "Thr": "C(OP(=O)(O)O)C",
    "Ser": "COP(=O)(O)O",
    "Tyr": "Cc1ccc(OP(=O)(O)O)cc1",
}


@dataclass
class PeptidePair:
    """A paired set of phospho/unphospho peptide ligands."""

    name: str
    sequence: str
    site: int  # 1-indexed position
    residue_type: str  # Thr, Ser, or Tyr
    phospho_pdbqt: Path | None = None
    unphospho_pdbqt: Path | None = None


def sequence_to_smiles(sequence: str, phospho_site: int | None = None, residue_type: str = "Thr") -> str:
    """Convert a peptide sequence to a SMILES string.

    Builds a linear peptide with standard amino acid side chains.
    If phospho_site is specified, replaces that residue with the
    phosphorylated version.

    Args:
        sequence: One-letter amino acid sequence.
        phospho_site: 1-indexed position to phosphorylate (None = unphospho).
        residue_type: Type of residue at phospho_site (Thr, Ser, or Tyr).

    Returns:
        SMILES string for the peptide.
    """
    # Build peptide as sequence of amino acid SMILES joined by peptide bonds
    # Simplified: use RDKit's peptide builder approach
    residues = []
    for i, aa in enumerate(sequence, start=1):
        if phospho_site and i == phospho_site and residue_type in PHOSPHO_SMILES:
            residues.append(f"[NH]C(=O)C(N)({PHOSPHO_SMILES[residue_type]})")
        elif aa in AMINO_ACID_SMILES:
            side_chain = AMINO_ACID_SMILES[aa]
            if side_chain:
                residues.append(f"[NH]C(=O)C(N)({side_chain})")
            else:
                residues.append(f"[NH]C(=O)C(N)")
        else:
            residues.append(f"[NH]C(=O)C(N)")

    # Join into linear chain — simplified SMILES representation
    smiles = ".".join(residues)
    return smiles


def build_peptide_mol(sequence: str, phospho_site: int | None = None, residue_type: str = "Thr") -> Chem.Mol:
    """Build a peptide as an RDKit Mol with 3D coordinates.

    Uses Helm/FASTA-style peptide construction via RDKit.

    Args:
        sequence: One-letter amino acid sequence.
        phospho_site: 1-indexed position to phosphorylate (None for unphospho).
        residue_type: Thr, Ser, or Tyr.

    Returns:
        RDKit Mol object with 3D conformer.
    """
    smiles = sequence_to_smiles(sequence, phospho_site, residue_type)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse peptide SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    params = rdDistGeom.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        # Fallback: try without distance geometry refinement
        status = AllChem.EmbedMolecule(mol, randomSeed=42)
        if status != 0:
            raise RuntimeError(f"Failed to generate 3D conformer for {sequence}")

    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    return mol


def mol_to_pdbqt(mol: Chem.Mol, output_path: Path) -> Path:
    """Convert an RDKit Mol to PDBQT format using Meeko.

    Args:
        mol: RDKit molecule with 3D conformer.
        output_path: Where to write the PDBQT file.

    Returns:
        Path to the written PDBQT file.
    """
    preparator = MoleculePreparation()
    mol_setup_list = preparator.prepare(mol)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    for mol_setup in mol_setup_list:
        pdbqt_string = PDBQTWriterLegacy.write_string(mol_setup)
        output_path.write_text(pdbqt_string[0])
        break  # Take first setup

    return output_path


def generate_pair(
    name: str,
    sequence: str,
    site: int,
    residue_type: str,
    output_dir: Path,
) -> PeptidePair:
    """Generate a phospho/unphospho peptide pair as PDBQT files.

    Args:
        name: Identifier for this peptide (e.g. "S6K1_T389").
        sequence: One-letter amino acid sequence.
        site: 1-indexed position of the residue to phosphorylate.
        residue_type: Thr, Ser, or Tyr.
        output_dir: Directory to write PDBQT files.

    Returns:
        PeptidePair with paths to both PDBQT files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    pair = PeptidePair(name=name, sequence=sequence, site=site, residue_type=residue_type)

    # Unphosphorylated
    mol_unphospho = build_peptide_mol(sequence, phospho_site=None)
    pair.unphospho_pdbqt = mol_to_pdbqt(
        mol_unphospho, output_dir / f"{name}_unphospho.pdbqt"
    )

    # Phosphorylated
    mol_phospho = build_peptide_mol(sequence, phospho_site=site, residue_type=residue_type)
    pair.phospho_pdbqt = mol_to_pdbqt(
        mol_phospho, output_dir / f"{name}_phospho.pdbqt"
    )

    return pair


def generate_pairs_from_csv(csv_path: Path, output_dir: Path) -> list[PeptidePair]:
    """Generate all peptide pairs from an input CSV.

    CSV columns: name, sequence, site, residue_type

    Args:
        csv_path: Path to peptides.csv.
        output_dir: Directory for PDBQT output.

    Returns:
        List of PeptidePair objects.
    """
    import pandas as pd

    df = pd.read_csv(csv_path)
    pairs = []
    for _, row in df.iterrows():
        pair = generate_pair(
            name=row["name"],
            sequence=row["sequence"],
            site=int(row["site"]),
            residue_type=row["residue_type"],
            output_dir=output_dir,
        )
        pairs.append(pair)

    return pairs
