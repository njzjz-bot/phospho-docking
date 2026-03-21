# Project Plan: Phosphopeptide Docking Pipeline

## Problem

Given a protein construct with known 3D coordinates, determine whether it selectively
binds phosphorylated peptides over unphosphorylated ones.

Wet lab validation (e.g. reporter assays) is slow and expensive. In silico docking with
AutoDock Vina can screen candidates computationally first, prioritising which peptides
to test experimentally.

## Approach

Dock **paired** peptides (phospho vs unphospho) against the receptor. The difference in
binding energy (ΔΔG) measures selectivity — how well the construct discriminates the
phospho "on" state from the unphospho "off" state.

```
                         ┌─────────────────┐
  Receptor (.pdb) ──────►  Receptor Prep  ├──┐
                         └─────────────────┘  │
                                              ▼
                         ┌─────────────────┐  ┌──────────┐     ┌──────────┐
  Peptide CSV ──────────►│ Peptide Builder  │──► Vina Dock ├────► Analysis │──► Report
                         │ (phospho pair)   │  └──────────┘     └──────────┘
                         └─────────────────┘
```

## Selectivity metric

| Quantity | Definition |
|----------|-----------|
| ΔG_phospho | Best Vina score for phosphorylated peptide (kcal/mol) |
| ΔG_unphospho | Best Vina score for unphosphorylated peptide (kcal/mol) |
| **ΔΔG** | ΔG_phospho − ΔG_unphospho (negative = prefers phospho = good) |

## Implementation steps

### Step 1 — Project scaffold
- [x] `uv init` with Python 3.12
- [x] Git repo
- [x] Directory structure: `phospho_docking/`, `data/`, `results/`
- [ ] Add dependencies to `pyproject.toml`
- [ ] `uv sync` to install

### Step 2 — Receptor preparation (`phospho_docking/prepare_receptor.py`)
- Load receptor PDB file
- Strip waters, add polar hydrogens
- Convert to PDBQT via OpenBabel
- Auto-compute bounding box around protein centre (with padding)
- Save receptor PDBQT + box config

### Step 3 — Peptide generation (`phospho_docking/generate_peptides.py`)
- Read peptide CSV (name, sequence, site, residue_type)
- For each entry, build two SMILES strings:
  - Unphosphorylated peptide
  - Phosphorylated peptide (pThr / pSer / pTyr at specified site)
- Generate 3D conformer (RDKit `EmbedMolecule` + `MMFFOptimizeMolecule`)
- Convert to PDBQT via Meeko `MoleculePreparation`

### Step 4 — Docking engine (`phospho_docking/dock.py`)
- Load receptor PDBQT + box config
- For each ligand PDBQT, run Vina:
  - `exhaustiveness=32` (good balance of speed vs accuracy)
  - `n_poses=5`
- Collect top binding affinity per ligand
- Save pose PDBQT files

### Step 5 — Selectivity analysis (`phospho_docking/analyze.py`)
- Pair up phospho/unphospho results by peptide name
- Compute ΔΔG for each pair
- Rank by selectivity (most negative ΔΔG first)
- Flag pairs with ΔΔG ≈ 0 (no discrimination)
- Output `results/scores.csv`

### Step 6 — Report (`phospho_docking/report.py`)
- Grouped bar chart: phospho vs unphospho binding energy per peptide
- Selectivity bar chart: ΔΔG per peptide
- Summary table in HTML
- PyMOL `.pml` scripts for best poses
- Output `results/report.html`

### Step 7 — CLI (`phospho_docking/cli.py`)
- `phospho-dock run --receptor PDB --peptides CSV --output DIR`
- `phospho-dock test --receptor PDB --sequence SEQ --site N --residue TYPE`
- Orchestrates steps 2–6

## Resources required from user

| What | Where to put it |
|------|----------------|
| Receptor PDB (experimental or predicted) | `data/receptors/receptor.pdb` |
| Peptide definitions (sequence + phosphosite) | `data/peptides/peptides.csv` |

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| vina | ≥1.2.5 | Docking engine (Python bindings) |
| meeko | ≥0.5 | Ligand/receptor → PDBQT |
| rdkit | ≥2024.03 | Build peptide molecules, conformers |
| openbabel-wheel | ≥3.1 | Format conversion fallback |
| biopython | ≥1.84 | PDB parsing |
| pandas | ≥2.0 | Results tabulation |
| matplotlib | ≥3.8 | Figures |
| jinja2 | ≥3.1 | HTML report templating |

## Compute estimate

| Peptide pairs | Exhaustiveness | Approx. time |
|---------------|---------------|-------------|
| 1 | 32 | 2-5 min |
| 10 | 32 | 30-60 min |
| 50 | 32 | 3-5 hours |

No GPU required. Single CPU.

## Verification

1. Run with a known phospho-binding system as positive control
2. Phospho scores should be consistently more negative than unphospho
3. Visual inspection: does the phosphate group contact expected residues?
4. Compare ΔΔG magnitudes with literature values for similar systems

## Open decision

**Peptide conformer strategy** (to decide during implementation):
- **Single conformer** — fastest, may miss better binding poses
- **Ensemble (5-10 conformers)** — more thorough, 5-10x slower
- **Flexible residues in Vina** — compromise, Vina handles flexibility at the phosphosite
