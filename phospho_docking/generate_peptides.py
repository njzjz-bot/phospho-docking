"""Generate phosphorylated and unphosphorylated peptide pairs as 3D PDBQT ligands."""

from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

from meeko import MoleculePreparation, PDBQTWriterLegacy
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdmolops


# Side chain SMARTS for residues that can be phosphorylated
# Used to find and modify the hydroxyl group
PHOSPHO_SMARTS: dict[str, tuple[str, str]] = {
    # (pattern to find OH, replacement SMILES for phospho)
    "Ser": ("[CH2][OH]", "[CH2]OP(=O)(O)O"),
    "Thr": ("[CH]([CH3])[OH]", "[CH]([CH3])OP(=O)(O)O"),
    "Tyr": ("c1cc(O)ccc1", "c1cc(OP(=O)(O)O)ccc1"),
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


def build_peptide_pdb(sequence: str) -> str:
    """Build a linear peptide as PDB text using Open Babel.

    Args:
        sequence: One-letter amino acid sequence.

    Returns:
        PDB file content as a string.
    """
    # Use Open Babel to build a peptide from sequence
    obconv = openbabel.OBConversion()
    obconv.SetInFormat("smi")
    obconv.SetOutFormat("pdb")

    # Build peptide using Open Babel's peptide builder
    mol = openbabel.OBMol()

    # Open Babel doesn't have a direct peptide builder in the Python API,
    # so we use the FASTA → 3D approach via a temporary file
    obconv_fasta = openbabel.OBConversion()
    obconv_fasta.SetInFormat("fasta")
    obconv_fasta.SetOutFormat("pdb")

    fasta_str = f">peptide\n{sequence}\n"
    obconv_fasta.ReadString(mol, fasta_str)

    # Generate 3D coordinates
    builder = openbabel.OBBuilder()
    builder.Build(mol)

    # Add hydrogens
    mol.AddHydrogens()

    # Energy minimise with MMFF94
    ff = openbabel.OBForceField.FindForceField("MMFF94")
    if ff and ff.Setup(mol):
        ff.ConjugateGradients(500)
        ff.GetCoordinates(mol)

    return obconv_fasta.WriteString(mol)


def pdb_to_rdkit_mol(pdb_text: str) -> Chem.Mol:
    """Convert PDB text to an RDKit Mol.

    Args:
        pdb_text: PDB file content as string.

    Returns:
        RDKit Mol with 3D conformer.
    """
    mol = Chem.MolFromPDBBlock(pdb_text, removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError("Failed to parse PDB into RDKit Mol")
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.AtomValenceException:
        # Some PDB structures have valence issues — try partial sanitization
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                         | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                         | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                         | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                         | Chem.SanitizeFlags.SANITIZE_SYMMRINGS)
    return mol


def phosphorylate_pdb(pdb_text: str, site: int, residue_type: str) -> str:
    """Add a phosphate group to a specific residue in a peptide PDB.

    Uses Open Babel to modify the residue at the given site.

    Args:
        pdb_text: PDB text of the unphosphorylated peptide.
        site: 1-indexed residue position to phosphorylate.
        residue_type: Thr, Ser, or Tyr.

    Returns:
        PDB text with phosphorylated residue.
    """
    obconv = openbabel.OBConversion()
    obconv.SetInAndOutFormats("pdb", "pdb")

    mol = openbabel.OBMol()
    obconv.ReadString(mol, pdb_text)

    # Find the target residue's hydroxyl oxygen
    target_oxygen = None
    for atom in openbabel.OBMolAtomIter(mol):
        res = atom.GetResidue()
        if res is None:
            continue
        res_num = res.GetNum()
        if res_num != site:
            continue
        # Look for the side chain oxygen (OG for Ser, OG1 for Thr, OH for Tyr)
        atom_name = res.GetAtomID(atom).strip()
        if residue_type == "Ser" and atom_name == "OG":
            target_oxygen = atom
            break
        elif residue_type == "Thr" and atom_name == "OG1":
            target_oxygen = atom
            break
        elif residue_type == "Tyr" and atom_name == "OH":
            target_oxygen = atom
            break

    if target_oxygen is None:
        raise ValueError(
            f"Could not find hydroxyl oxygen for {residue_type} at position {site}"
        )

    # Find and remove the hydrogen on the target oxygen
    h_to_remove = None
    for neighbor in openbabel.OBAtomAtomIter(target_oxygen):
        if neighbor.GetAtomicNum() == 1:  # Hydrogen
            h_to_remove = neighbor
            break

    if h_to_remove:
        mol.DeleteAtom(h_to_remove)

    # Add phosphorus atom
    p_atom = mol.NewAtom()
    p_atom.SetAtomicNum(15)  # P

    # Position P near the oxygen
    ox = target_oxygen.GetX()
    oy = target_oxygen.GetY()
    oz = target_oxygen.GetZ()
    p_atom.SetVector(ox + 1.6, oy, oz)

    # Bond O-P
    mol.AddBond(target_oxygen.GetIdx(), p_atom.GetIdx(), 1)

    # Add three oxygens to phosphorus (=O, -OH, -OH)
    for i, (dx, dy, dz) in enumerate([(0, 1.5, 0), (1.5, 0, 0), (0, -1.5, 0)]):
        o_atom = mol.NewAtom()
        o_atom.SetAtomicNum(8)
        o_atom.SetVector(ox + 1.6 + dx, oy + dy, oz + dz)

        if i == 0:
            # Double bond (P=O)
            mol.AddBond(p_atom.GetIdx(), o_atom.GetIdx(), 2)
        else:
            # Single bond (P-OH)
            mol.AddBond(p_atom.GetIdx(), o_atom.GetIdx(), 1)
            h_atom = mol.NewAtom()
            h_atom.SetAtomicNum(1)
            h_atom.SetVector(ox + 1.6 + dx + 0.96, oy + dy, oz + dz)
            mol.AddBond(o_atom.GetIdx(), h_atom.GetIdx(), 1)

    # Minimise to relax the phosphate geometry
    ff = openbabel.OBForceField.FindForceField("MMFF94")
    if ff and ff.Setup(mol):
        ff.ConjugateGradients(200)
        ff.GetCoordinates(mol)

    return obconv.WriteString(mol)


def mol_to_pdbqt_via_obabel(pdb_text: str, output_path: Path) -> Path:
    """Convert peptide PDB text to PDBQT format via Open Babel.

    Args:
        pdb_text: PDB file content.
        output_path: Where to write the PDBQT file.

    Returns:
        Path to the written PDBQT file.
    """
    obconv = openbabel.OBConversion()
    obconv.SetInAndOutFormats("pdb", "pdbqt")
    obconv.AddOption("h", openbabel.OBConversion.OUTOPTIONS)

    mol = openbabel.OBMol()
    obconv.ReadString(mol, pdb_text)

    # Assign Gasteiger charges
    charge_model = openbabel.OBChargeModel.FindType("gasteiger")
    if charge_model:
        charge_model.ComputeCharges(mol)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    obconv.WriteFile(mol, str(output_path))

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
        name: Identifier for this peptide.
        sequence: One-letter amino acid sequence.
        site: 1-indexed position of the residue to phosphorylate.
        residue_type: Thr, Ser, or Tyr.
        output_dir: Directory to write PDBQT files.

    Returns:
        PeptidePair with paths to both PDBQT files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    pair = PeptidePair(name=name, sequence=sequence, site=site, residue_type=residue_type)

    # Build unphosphorylated peptide
    pdb_text = build_peptide_pdb(sequence)

    # Unphosphorylated → PDBQT
    pair.unphospho_pdbqt = mol_to_pdbqt_via_obabel(
        pdb_text, output_dir / f"{name}_unphospho.pdbqt"
    )

    # Phosphorylated → modify PDB → PDBQT
    phospho_pdb = phosphorylate_pdb(pdb_text, site, residue_type)
    pair.phospho_pdbqt = mol_to_pdbqt_via_obabel(
        phospho_pdb, output_dir / f"{name}_phospho.pdbqt"
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
