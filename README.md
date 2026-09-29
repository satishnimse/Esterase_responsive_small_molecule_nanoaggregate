# Molecular dynamics data — Esterase-responsive small-molecule nanoaggregates (BA/BE benzimidazoles)

Analysis data supporting the molecular dynamics (MD) results of the manuscript
*"An Esterase-Responsive Molecular Switch Controls the Biological Consequences of Small-Molecule
Nanoaggregation"* (Yeom, Park, *et al.*).

The archive contains two independent simulation campaigns:

| Folder | Campaign | Manuscript |
|---|---|---|
| `MD_simulation_Protein_ligand/` | 16 protein–ligand complexes: BA-3/7/8/10 and BE-3/7/8/10 at the colchicine (CBS, PDB 4O2B) and nocodazole (NBS, PDB 5CA1) sites of tubulin, 100 ns each, plus MM-GBSA/MM-PBSA | Section 2.2, Figure 2, Figures S83–S87 |
| `MD_simulation_Self-assembly/` | Explicit-solvent self-assembly of BE-7, neutral BA-7 and the BA-7 carboxylate (BA-7C), 27 molecules per box, 3 replicates × 100 ns each | Section 2.3, Figure 3 |

All raw `.xvg`, `.dat` and `.csv` files are unmodified GROMACS, gmx_MMPBSA or analysis-script
output. Trajectories (`.xtc`), run inputs (`.tpr`) and topologies are not included because of their
size; they are available from the corresponding author on request.

---

## 1. `MD_simulation_Protein_ligand/`

### 1.1 Systems and protocol

- **Receptors:** tubulin αβ-heterodimer, 4O2B with colchicine removed (`4O2B_CLC`) and 5CA1 with
  nocodazole removed (`5CA1_NOZ`). Starting ligand poses are the top-ranked docked poses.
- **Ligands:** BA-3, BA-7, BA-8, BA-10 (carboxylic acids) and BE-3, BE-7, BE-8, BE-10 (methyl
  esters). All ligand topologies are net-neutral; the BA carboxylic acids were modelled as the
  **neutral COOH**.
- **Force fields:** CHARMM36 (protein, ions); ligand parameters CGenFF/SwissParam; TIP3P water;
  Na⁺/Cl⁻ at 0.15 M.
- **Equilibration:** steepest-descent minimization; 500 ps NVT (300 K, V-rescale) and 1 ns NPT
  (1 bar, C-rescale), both with protein and ligand position restraints.
- **Production:** 100 ns, unrestrained, 2 fs time step, 300 K (V-rescale), 1 bar (Parrinello–Rahman).
- **Nonbonded settings:** PME electrostatics (1.2 nm real-space cutoff); van der Waals force-switch
  1.0–1.2 nm; LINCS on bonds to hydrogen.
- **Output:** solute-only coordinates every 100 ps (1001 frames per trajectory). Water and ions are
  therefore absent from all analysis output. `energy.xvg` is written every 10 ps.
- **Binding energies:** gmx_MMPBSA 1.6.4, run once per system for GB, PB and per-residue
  decomposition on the same 91 frames (10–100 ns at 1 ns spacing), **no entropy term**.
  - GB: igb = 5, internal/external dielectric 1.0/78.5, 0.15 M salt, surface tension 0.0072.
  - PB: ipb = 2, inp = 1, internal/external dielectric 1.0/80.0, ionic strength 0.15 M,
    grid scale 2.0, fill ratio 4.0, PBRadii = 3, 300 K.
  - Decomposition: idecomp = 2, residues within 6 Å of the ligand.

  The values are effective interaction energies for comparison within the series, not absolute
  binding free energies.

### 1.2 Files

```
MD_simulation/<complex>/<ligand>/        <complex> = 4O2B_CLC | 5CA1_NOZ
    rmsd_backbone.xvg    protein backbone RMSD vs time (nm)
    rmsd_ligand.xvg      ligand RMSD after fitting on the protein backbone (nm)
    rmsf_protein.xvg     per-residue protein RMSF (nm)
    rmsf_ligand.xvg      per-atom ligand RMSF (nm)
    gyrate_protein.xvg   protein radius of gyration vs time (nm)
    gyrate_ligand.xvg    ligand radius of gyration vs time (nm)
    hb.xvg               protein–ligand hydrogen-bond count vs time
    pl_contacts.xvg      protein–ligand atom contacts < 0.35 nm vs time
    pl_min_dist.xvg      protein–ligand minimum distance vs time (nm)
    sasa_ligand.xvg      bound-state ligand SASA (nm²); column 2 = complex, column 3 = ligand
    energy.xvg           potential energy, temperature, pressure, density vs time

MD_simulation/MD_summary_figures_data.xlsx
    Data plotted in Figure 2 and Figures S83–S86, one sheet per ligand pair and site
    (time series and per-residue RMSF). Distances are converted to Å (× 10) and
    areas to Å² (× 100).
MD_simulation/Ligand_SASA_data.xlsx
    Ligand SASA (Å²) vs time for all ligands, one sheet per complex.

MM-PBSA/<complex>/<ligand>/
    FINAL_RESULTS_MMPBSA.dat    GB and PB energy breakdown (complex / receptor / ligand / delta)
    FINAL_RESULTS_MMPBSA.csv    the same, per frame
    FINAL_DECOMP_MMPBSA.dat     per-residue decomposition (idecomp = 2, residues within 6 Å)
```

The ΔTOTAL values in `FINAL_RESULTS_MMPBSA.dat` are the data plotted in Figure S87.

Residues are numbered sequentially across the α- and β-tubulin chains of the simulated topology
(861 residues for 4O2B_CLC, 864 for 5CA1_NOZ). This is not the canonical UniProt numbering.

---

## 2. `MD_simulation_Self-assembly/`

### 2.1 Systems and protocol

| Tag | Species | Net charge | Atoms per molecule |
|---|---|---|---|
| `BE7` | BE-7 (methyl ester) | 0 | 37 |
| `BA7` | BA-7, neutral carboxylic acid | 0 | 34 |
| `BA7C` | BA-7 carboxylate (deprotonated) | −1 | 33 |

- **Starting structures:** DFT-optimized geometries (B3LYP/6-31+G(d,p)).
- **Parameters:** CGenFF 4.0, identical for all three species, with the CHARMM36 force field
  (including the Na⁺–carboxylate NBFIX corrections) and CHARMM-modified TIP3P water.
- **Box:** 27 molecules on a 3 × 3 × 3 lattice with random orientations in a cubic box
  (L ≈ 9.9–10.0 nm), giving 45 mM solute. Solvated, then 0.15 M NaCl plus neutralizing
  counterions were added (BA7C carries 27 extra Na⁺).
- **Equilibration:** steepest-descent minimization; 100 ps NVT and 1 ns NPT (C-rescale, 1 bar),
  both with solute heavy-atom restraints.
- **Production:** 100 ns, unrestrained, 2 fs time step, 310 K (V-rescale, single coupling group),
  1 bar. PME electrostatics (1.2 nm), van der Waals force-switch 1.0–1.2 nm, no dispersion
  correction. Engine: GROMACS 2026.3.
- **Replicates:** three per species. They differ in molecular orientations and initial
  velocities.
- **Output:** coordinates every 10 ps (10 001 frames).
- **Analysis window:** 50–100 ns unless noted. Reported values are the mean ± SD (ddof = 1)
  across the three replicates.

### 2.2 Key analysis definitions

- **Largest cluster, N_max:** `gmx clustsize -cut 0.35` on the solute atoms (single linkage,
  0.35 nm atom–atom cutoff). The atom count of the largest cluster is divided by the atoms per
  molecule to give N_max in molecules.
- **Coordination number:** center-of-mass RDF (`gmx rdf`, bin 0.02 nm), read at r = 0.7 nm.
- **SASA:** `gmx sasa` on the solute, split into carbon (nonpolar), oxygen (polar) and the
  remainder (H/N).
- **Diffusion:** Einstein fits of the solute MSD over 0–20 ns and 50–100 ns.
- **Escape rate and residence time:** based on membership of the largest cluster (0.35 nm,
  single linkage). An escape means a molecule leaves the largest cluster and stays out for more
  than 1 ns.
- **Hydrogen bonds:** donor–acceptor ≤ 0.35 nm and D–H···A ≥ 150°.
- **π-stacking:** centroid–centroid ≤ 0.6 nm and interplanar angle ≤ 30° between the fused
  aromatic cores.
- **Sampling of H-bonds and π-stacking:** computed on 1 frame/ns over 50–100 ns.

### 2.3 Files

```
data/<TAG>/<repN>/                   <TAG> = BE7 | BA7 | BA7C ; <repN> = rep1 | rep2 | rep3
    maxclust_atom.xvg     largest cluster size vs time (atoms; ÷ atoms/molecule = N_max)
    nclust_atom.xvg       number of clusters vs time
    avclust_atom.xvg      average cluster size vs time (atoms)
    histo-clust_atom.xvg  cluster-size distribution (atoms)
    rdf_com.xvg           solute COM–COM radial distribution function
    cn_com.xvg            cumulative coordination number from the RDF
    sasa_total.xvg        solute SASA vs time (nm²): total, carbon, oxygen
    gyration_shape.dat    largest-cluster Rg (nm), asphericity (nm²), relative shape anisotropy κ²
    rmsd_conf_timeseries.dat   per-molecule conformational RMSD vs own t = 0 structure (nm)
    rmsf_atom.dat         per-atom RMSF (nm); atom names in rmsf_atom_names.json
    sec6_dynamics.json    diffusion coefficients, escape events/rate, residence time
    sec7_mechanism.json   hydrogen bonds, π-stacking, group contacts, first-shell water

Simulation_raw_data_avg_figures.xlsx
    Replicate values, means and SDs behind every Figure 3 panel (12 sheets).
build_excel.py
    Regenerates the workbook from data/:  python build_excel.py
    (requires numpy, pandas, openpyxl)

MD_simulation_pdb_files/<TAG>_rep2_t{0,50,100}ns_full.pdb
    Representative full-system snapshots (solute, water, ions) of replicate 2 at 0, 50 and
    100 ns. Coordinates in Å.
```

---

## 3. Units

GROMACS files use GROMACS units (nm, nm², ps or ns as labelled in each file header). The Excel
workbooks convert to Å and Å² where the column header says so. Energies are in kcal mol⁻¹
(gmx_MMPBSA) or kJ mol⁻¹ (`energy.xvg`).

## 4. Contact

Satish Balasaheb Nimse — satish_nimse@hallym.ac.kr

## 5. License

This dataset is released under the Creative Commons Attribution 4.0 International license
(CC BY 4.0); see `LICENSE`.
