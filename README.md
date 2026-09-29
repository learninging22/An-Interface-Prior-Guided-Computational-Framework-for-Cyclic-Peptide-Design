# An Interface-Prior-Guided Computational Framework for Cyclic Peptide Design

> **Case study: the neuroblastoma-associated target GPC2**

This repository provides an interface-prior-guided computational framework for cyclic peptide design. It covers protein–protein interaction analysis, cyclic peptide backbone construction, sequence designability analysis, sequence design and structural refinement, multidimensional scoring, and molecular dynamics simulation.

---

## Repository Structure

```text
.
├── 1.ppi_analysis/
│   ├── 6WJL_fixed_GHL.pdb
│   ├── ddg.xml
│   └── pack_input_options.txt
│
├── 2.cyclic_peptide_backbone_example/
│   ├── 6WJL_fix_gh_mini_0006.pdb
│   ├── backbone_12.xml
│   └── rosetta.flags
│
├── 3.sequence_designability_analysis_scripts/
│   ├── step1_ligandmpnn_design.py
│   ├── step2_analyze_samples.py
│   └── step3_compute_metrics.py
│
├── 4.mpnn_design_rosetta_relax/
│   ├── ligandmpnn_design.py
│   ├── relax.xml
│   └── rosetta.flags
│
├── 5.multidimensional_scoring/
│   ├── 6WJL.trg
│   ├── adcp_dock.sh
│   └── rosetta_cms_sap_ddg_analysis.xml
│
├── 6.MD_outputs/ 
│   ├── 6WJL_pro/ 
│   ├── cycpep-165/ 
│   ├── cycpep-673/ 
│   ├── cycpep-762/ 
│   ├── cycpep-876/ 
│   ├── cycppep-1760/ 
│   ├── cp-165/ 
│   ├── cp-673/ 
│   ├── cp-1760/ 
│   ├── cp-876/ 
│   └── cp-762/
|
├── .gitattributes
└── README.md
```

---

## 1. `1.ppi_analysis`

This directory is used for protein–protein interaction interface analysis and key-residue evaluation.

| File | Purpose |
|---|---|
| `6WJL_fixed_GHL.pdb` | A repaired and Rosetta-processed three-dimensional structure of the GPC2 complex, used as the structural basis for interface analysis and residue mutation studies. |
| `ddg.xml` | A RosettaScripts protocol for residue mutation and structural relaxation. It mutates a selected residue to alanine and applies coordinate constraints and FastRelax for structural optimization and scoring. |
| `pack_input_options.txt` | A Rosetta InterfaceAnalyzer configuration file that controls side-chain packing, separated-state packing, interface statistics, solvent-accessible surface area calculations, and score reporting. |

---

## 2. `2.cyclic_peptide_backbone_example`

This directory provides a Rosetta example for constructing a 12-residue cyclic peptide backbone from an interface-derived peptide fragment.

| File | Purpose |
|---|---|
| `6WJL_fix_gh_mini_0006.pdb` | A target protein–peptide complex structure used as the starting model for cyclic peptide backbone construction. It contains the initial peptide fragment whose interface conformation is retained. |
| `backbone_12.xml` | A Rosetta protocol for generating a 12-residue cyclic peptide backbone. It selects target residues, extends the peptide chain, assigns backbone torsion angles, declares the terminal covalent bond, performs GeneralizedKIC closure, designs side chains, minimizes energy, relaxes the structure, and filters conformations. |
| `rosetta.flags` | A Rosetta parameter file for cyclic peptide backbone generation. It specifies the protocol file, number of generated structures, full-atom mode, cyclic connectivity output, and logging options. |

---

## 3. `3.sequence_designability_analysis_scripts`

This directory is used to compare sequence designability across cyclic peptide lengths and backbone structures.

| File | Purpose |
|---|---|
| `step1_ligandmpnn_design.py` | A LigandMPNN batch-design script that organizes sequence sampling tasks by cyclic peptide length, backbone, and random seed, while recording task status, progress, and execution logs. |
| `step2_analyze_samples.py` | A LigandMPNN result parsing and statistical analysis script. It extracts designed-chain sequences and confidence information, then calculates sequence duplication, amino-acid frequencies, dominant residues, and positional entropy. |
| `step3_compute_metrics.py` | A sequence designability analysis script that integrates backbone-level and position-level statistics. It calculates sequence recovery, effective amino-acid alphabet size, positional locking, chemical-class coverage, and bootstrap statistics across different cyclic peptide lengths. |

---

## 4. `4.mpnn_design_rosetta_relax`

This directory is used for LigandMPNN sequence design and Rosetta structural refinement of candidate cyclic peptides.

| File | Purpose |
|---|---|
| `ligandmpnn_design.py` | A LigandMPNN batch-task management script. It identifies PDB structures for design, defines design chains and fixed residues, constructs design tasks, supports array-task indexing, and handles skipping or overwriting existing results. |
| `relax.xml` | A RosettaScripts protocol for cyclic peptide conformational processing and optimization. It includes cyclic peptide extension and closure, backbone conformational sampling, side-chain design, energy minimization, FastRelax refinement, and filtering based on hydrogen bonding, steric clashes, and shape complementarity. |
| `rosetta.flags` | A Rosetta configuration file for structural refinement. It links the `relax.xml` protocol, controls the number of refinements per structure, sets full-atom mode, defines output naming and storage locations, and preserves cyclic connectivity information. |

---

## 5. `5.multidimensional_scoring`

This directory is used for Rosetta interface scoring and AutoDock CrankPep redocking of candidate cyclic peptides.

| File | Purpose |
|---|---|
| `6WJL.trg` | An ADCP target file generated from the GPC2 receptor. It stores receptor structural information and precomputed target fields required for cyclic peptide docking. |
| `adcp_dock.sh` | A batch ADCP cyclic peptide redocking script. It identifies candidate peptide structures, extracts amino-acid sequences, schedules parallel docking tasks, and summarizes execution status, runtime, and best predicted affinity. |
| `rosetta_cms_sap_ddg_analysis.xml` | A Rosetta multidimensional interface evaluation protocol. It preserves cyclic peptide closure, optimizes the binding interface, and calculates interface binding energy ΔΔG, contact molecular surface area (CMS), and spatial aggregation propensity (SAP). |

---

## 6. `6.MD_outputs`

This directory contains GROMACS molecular dynamics files for the protein reference system and representative cyclic peptide complexes.

### 6.1 `6WJL_pro`

This subdirectory represents the GPC2 protein reference system without a candidate cyclic peptide.

| File | Purpose |
|---|---|
| `6WJL_nocyc.pdb` | Atomic coordinates of the GPC2 protein without a cyclic peptide, used as the protein-only reference structure. |
| `topol.top` | The main GROMACS topology file for the protein reference system. It defines the force field, protein topology, water model, ions, and molecular composition of the system. |
| `posre.itp` | Position-restraint parameters for protein heavy atoms, used to restrict selected atoms during simulation preparation. |
| `md.mdp` | Production molecular dynamics parameters for the protein reference system, including integration, time step, constraints, neighbor searching, electrostatics, temperature and pressure coupling, and trajectory output settings. |

### 6.2 Cyclic Peptide Complex Subdirectories

The following subdirectories correspond to different candidate cyclic peptide–GPC2 complexes:

| Subdirectory | Purpose |
|---|---|
| `cycpep-165/` | Molecular dynamics system for candidate cyclic peptide 165 bound to GPC2. |
| `cycpep-673/` | Molecular dynamics system for candidate cyclic peptide 673 bound to GPC2. |
| `cycpep-762/` | Molecular dynamics system for candidate cyclic peptide 762 bound to GPC2. |
| `cycpep-876/` | Molecular dynamics system for candidate cyclic peptide 876 bound to GPC2. |
| `cycppep-1760/` | Molecular dynamics system for candidate cyclic peptide 1760 bound to GPC2. |

These cyclic peptide complex subdirectories follow the same file organization:

| File | Purpose |
|---|---|
| `pro.pdb` | Atomic coordinates of the GPC2 target protein. |
| `cp.pdb` | Atomic coordinates of the corresponding candidate cyclic peptide. |
| `cp.itp` | The GROMACS molecular topology file for the cyclic peptide. It defines atom types, charges, masses, bonds, pair interactions, angles, dihedrals, and CMAP correction terms. |
| `topol.top` | The main GROMACS topology file for the cyclic peptide–protein complex. It integrates the protein, cyclic peptide, force field, water model, ions, and position-restraint definitions. |
| `posre.itp` | Position-restraint parameters for GPC2 protein heavy atoms. |
| `posre_cp.itp` | Position-restraint parameters for atoms in the candidate cyclic peptide. |
| `md.mdp` | Production molecular dynamics parameters for the cyclic peptide–protein complex, including integration, constraints, nonbonded interactions, temperature and pressure coupling, and trajectory output settings. |

---

