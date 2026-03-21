# Phosphopeptide Docking Pipeline

In silico screening of phosphopeptide–protein interactions using AutoDock Vina.
Compares phosphorylated vs unphosphorylated peptides to measure binding selectivity (on/off discrimination).

## Motivation

Screen whether a protein construct can discriminate phosphorylated from unphosphorylated peptides computationally, before committing to wet lab validation.

## Pipeline

```
Receptor (PDB) ──► Receptor prep (PDBQT)
                                               ├──► AutoDock Vina ──► Scores + Poses
Peptide sequences ──► 3D structure ──► Ligand prep (PDBQT)              │
  (phospho + unphospho pairs)                                           ▼
                                                               Selectivity Report
                                                          (ΔΔG phospho vs unphospho)
```

### Modules

| Module | File | Purpose |
|--------|------|---------|
| Receptor prep | `phospho_docking/prepare_receptor.py` | PDB → PDBQT, define search box |
| Peptide generation | `phospho_docking/generate_peptides.py` | Build phospho/unphospho pairs as 3D conformers |
| Docking | `phospho_docking/dock.py` | Run AutoDock Vina, collect affinities |
| Analysis | `phospho_docking/analyze.py` | Compute ΔΔG selectivity, rank peptides |
| Report | `phospho_docking/report.py` | Generate HTML report with figures |
| CLI | `phospho_docking/cli.py` | Command-line entry point |

### Selectivity metric

For each peptide pair:

- **ΔG_phospho** — best Vina score for phosphorylated peptide (kcal/mol)
- **ΔG_unphospho** — best Vina score for unphosphorylated peptide (kcal/mol)
- **ΔΔG = ΔG_phospho − ΔG_unphospho** — negative means construct prefers phospho (good)

## Example peptide set

1. Define peptides as short windows (±5-7 residues) around a known phosphosite
2. Supports pThr, pSer, and pTyr modifications

## Resources you need to provide

| Resource | Format | Notes |
|----------|--------|-------|
| Receptor structure | `.pdb` | PDB file (experimental or predicted, e.g. AlphaFold) |
| Peptide sequences | text | Residues around the phosphosite (e.g. 7-mer or 11-mer window) |

## Software dependencies (installed automatically via `uv sync`)

| Package | Purpose |
|---------|---------|
| `vina` | AutoDock Vina Python bindings — docking engine |
| `meeko` | Prepares ligands for Vina (PDBQT conversion) |
| `rdkit` | Cheminformatics — build peptide structures, add phospho groups |
| `openbabel-wheel` | 3D coordinate generation, format conversion |
| `biopython` | Parse PDB files, extract sequences |
| `pandas` | Tabulate and analyze docking results |
| `matplotlib` | Binding affinity plots |
| `jinja2` | HTML report templating |

## Optional external tools (for visualization)

| Tool | Purpose | Install |
|------|---------|---------|
| PyMOL | View binding poses | `conda install -c conda-forge pymol-open-source` |
| ADFR Suite | `prepare_receptor` CLI (alternative receptor prep) | https://ccsb.scripps.edu/adfr/ |

## Compute requirements

- Each Vina docking run: ~1-5 min per peptide (single CPU)
- Full pipeline for ~10 peptide pairs: ~30-60 min
- No GPU required

## Usage

```bash
# Install
uv sync

# Run pipeline
phospho-dock run --receptor receptor.pdb --peptides peptides.csv --output results/

# Quick test with single peptide pair
phospho-dock test --receptor receptor.pdb --sequence "FLGFTYVAP" --site 5 --residue Thr
```

## Input format

`peptides.csv`:

```csv
name,sequence,site,residue_type
example_pThr,FLGFTYVAP,5,Thr
```

- `site`: 1-indexed position of the residue to phosphorylate
- `residue_type`: Thr, Ser, or Tyr

## Output

```
results/
├── scores.csv              # All docking scores + ΔΔG
├── poses/                  # Best PDBQT poses per peptide
├── figures/                # Bar charts, selectivity plots
├── pymol/                  # .pml visualization scripts
└── report.html             # Full summary report
```
