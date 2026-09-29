"""
Regenerates Simulation_raw_data_avg_figures.xlsx from the raw .xvg / .dat / .json
files in ./data/<SYSTEM>/<repN>/. Run with: python3 build_excel.py

Each per-metric sheet reports the three replicate values, their mean, and the
sample standard deviation (ddof=1) across replicates.
"""
import json
import numpy as np
import pandas as pd

BASE = "data"
OUT = "Simulation_raw_data_avg_figures.xlsx"

SYS = ["BE7", "BA7", "BA7C"]
REPS = ["rep1", "rep2", "rep3"]
ATOMS_PER_MOL = {"BE7": 37, "BA7": 34, "BA7C": 33}

writer = pd.ExcelWriter(OUT, engine="openpyxl")

# ============================================================
# Sheet 1: Aggregation N_max(t)  [fig 01]
# ============================================================
t_ns = np.loadtxt(f"{BASE}/BE7/rep1/maxclust_atom.xvg", comments=("#", "@"))[:, 0] / 1000.0
minlen = None
data_by_sys = {}
for s in SYS:
    arrs = []
    for r in REPS:
        d = np.loadtxt(f"{BASE}/{s}/{r}/maxclust_atom.xvg", comments=("#", "@"))
        arrs.append(d[:, 1] / ATOMS_PER_MOL[s])
    L = min(len(a) for a in arrs)
    data_by_sys[s] = [a[:L] for a in arrs]
    minlen = L if minlen is None else min(minlen, L)
df = pd.DataFrame({"Time_ns": t_ns[:minlen]})
for s in SYS:
    for i, r in enumerate(REPS):
        df[f"{s}_{r}_Nmax_molecules"] = data_by_sys[s][i][:minlen]
    arr = np.array([a[:minlen] for a in data_by_sys[s]])
    df[f"{s}_mean"] = arr.mean(axis=0)
    df[f"{s}_SD"] = arr.std(axis=0, ddof=1)
df.to_excel(writer, sheet_name="01_Aggregation_Nmax", index=False)

# ============================================================
# Sheet 2: RMSD_conf(t)  [fig 02]
# ============================================================
t_rc = np.loadtxt(f"{BASE}/BE7/rep1/rmsd_conf_timeseries.dat", comments="#")[:, 0]
data_by_sys = {}
minlen = None
for s in SYS:
    arrs = []
    for r in REPS:
        d = np.loadtxt(f"{BASE}/{s}/{r}/rmsd_conf_timeseries.dat", comments="#")
        arrs.append(d[:, 1] * 10.0)  # nm -> Angstrom
    L = min(len(a) for a in arrs)
    data_by_sys[s] = [a[:L] for a in arrs]
    minlen = L if minlen is None else min(minlen, L)
df = pd.DataFrame({"Time_ns": t_rc[:minlen]})
for s in SYS:
    for i, r in enumerate(REPS):
        df[f"{s}_{r}_RMSDconf_Angstrom"] = data_by_sys[s][i][:minlen]
    arr = np.array([a[:minlen] for a in data_by_sys[s]])
    df[f"{s}_mean"] = arr.mean(axis=0)
    df[f"{s}_SD"] = arr.std(axis=0, ddof=1)
df.to_excel(writer, sheet_name="02_RMSD_conf_timeseries", index=False)

# ============================================================
# Sheet 3: RMSF_atom  [fig 03]  -- atom counts differ per system, block layout
# ============================================================
frames = []
for s in SYS:
    arrs = []
    names = json.load(open(f"{BASE}/{s}/rep1/rmsf_atom_names.json"))
    for r in REPS:
        d = np.loadtxt(f"{BASE}/{s}/{r}/rmsf_atom.dat", comments="#")
        arrs.append(d[:, 1])
    arr = np.array(arrs)
    sub = pd.DataFrame({
        "System": s,
        "Atom_index": np.arange(arr.shape[1]),
        "Atom_name": names[:arr.shape[1]],
        f"{s}_rep1_RMSF_nm": arrs[0],
        f"{s}_rep2_RMSF_nm": arrs[1],
        f"{s}_rep3_RMSF_nm": arrs[2],
        "mean_RMSF_nm": arr.mean(axis=0),
        "SD_RMSF_nm": arr.std(axis=0, ddof=1),
    })
    frames.append(sub)
df = pd.concat(frames, ignore_index=True, sort=False)
df.to_excel(writer, sheet_name="03_RMSF_atom_profile", index=False)

# ============================================================
# Sheet 4: RDF COM-COM  [fig 04]
# ============================================================
r_axis = np.loadtxt(f"{BASE}/BE7/rep1/rdf_com.xvg", comments=("#", "@"))[:, 0]
data_by_sys = {}
minlen = None
for s in SYS:
    arrs = []
    for r in REPS:
        d = np.loadtxt(f"{BASE}/{s}/{r}/rdf_com.xvg", comments=("#", "@"))
        arrs.append(d[:, 1])
    L = min(len(a) for a in arrs)
    data_by_sys[s] = [a[:L] for a in arrs]
    minlen = L if minlen is None else min(minlen, L)
df = pd.DataFrame({"r_nm": r_axis[:minlen]})
for s in SYS:
    for i, r in enumerate(REPS):
        df[f"{s}_{r}_g_r"] = data_by_sys[s][i][:minlen]
    arr = np.array([a[:minlen] for a in data_by_sys[s]])
    df[f"{s}_mean"] = arr.mean(axis=0)
    df[f"{s}_SD"] = arr.std(axis=0, ddof=1)
df.to_excel(writer, sheet_name="04_RDF_COMCOM", index=False)

# ============================================================
# Sheet 5: Rg(t)  [fig 05]
# ============================================================
data_by_sys = {}
minlen = None
for s in SYS:
    arrs = []
    for r in REPS:
        d = np.loadtxt(f"{BASE}/{s}/{r}/gyration_shape.dat", comments="#")
        arrs.append(d[:, 0])
    L = min(len(a) for a in arrs)
    data_by_sys[s] = [a[:L] for a in arrs]
    minlen = L if minlen is None else min(minlen, L)
t_rg = np.linspace(0, 100, minlen)
df = pd.DataFrame({"Time_ns": t_rg})
for s in SYS:
    for i, r in enumerate(REPS):
        df[f"{s}_{r}_Rg_nm"] = data_by_sys[s][i][:minlen]
    arr = np.array([a[:minlen] for a in data_by_sys[s]])
    df[f"{s}_mean"] = arr.mean(axis=0)
    df[f"{s}_SD"] = arr.std(axis=0, ddof=1)
df.to_excel(writer, sheet_name="05_Rg_timeseries", index=False)

# ============================================================
# Sheet 6: SASA  [fig 06]  -- long format, per-rep window-mean values
# ============================================================
rows = []
for s in SYS:
    for r in REPS:
        d = np.loadtxt(f"{BASE}/{s}/{r}/sasa_total.xvg", comments=("#", "@"))
        total, nonpolar_c, polar_o = d[:, 1].mean(), d[:, 2].mean(), d[:, 3].mean()
        other = total - nonpolar_c - polar_o
        rows.append({"System": s, "Replicate": r, "Window_ns": "50-100",
                      "Total_SASA_nm2": total, "Nonpolar_C_SASA_nm2": nonpolar_c,
                      "Polar_O_SASA_nm2": polar_o, "Other_HN_SASA_nm2": other})
df = pd.DataFrame(rows)
summary = df.groupby("System")[["Total_SASA_nm2", "Nonpolar_C_SASA_nm2",
                                 "Polar_O_SASA_nm2", "Other_HN_SASA_nm2"]].agg(["mean", "std"])
summary.columns = ["_".join(c) for c in summary.columns]
summary = summary.reset_index()
summary.insert(1, "Replicate", "mean_SD")
df = pd.concat([df, summary], ignore_index=True, sort=False)
df.to_excel(writer, sheet_name="06_SASA", index=False)

# ============================================================
# Sheet 7: Coordination number  [fig 07]
# ============================================================
rows = []
for s in SYS:
    for r in REPS:
        d = np.loadtxt(f"{BASE}/{s}/{r}/cn_com.xvg", comments=("#", "@"))
        idx = np.argmin(np.abs(d[:, 0] - 0.7))
        rows.append({"System": s, "Replicate": r, "r_eval_nm": d[idx, 0],
                      "Coordination_number": d[idx, 1]})
df = pd.DataFrame(rows)
summary = df.groupby("System")["Coordination_number"].agg(["mean", "std"]).reset_index()
summary.insert(1, "Replicate", "mean_SD")
summary["r_eval_nm"] = 0.7
df = pd.concat([df, summary], ignore_index=True, sort=False)
df.to_excel(writer, sheet_name="07_Coordination_number", index=False)

# ============================================================
# Sheets 8-12: sec6/sec7 derived scalars  [figs 08, 09]
# ============================================================
sec6 = {s: [json.load(open(f"{BASE}/{s}/{r}/sec6_dynamics.json")) for r in REPS] for s in SYS}
sec7 = {s: [json.load(open(f"{BASE}/{s}/{r}/sec7_mechanism.json")) for r in REPS] for s in SYS}


def with_summary(df, value_cols):
    summary = df.groupby("System")[value_cols].agg(["mean", "std"])
    summary.columns = ["_".join(c) for c in summary.columns]
    summary = summary.reset_index()
    summary.insert(1, "Replicate", "mean_SD")
    return pd.concat([df, summary], ignore_index=True, sort=False)


rows = []
for s in SYS:
    for i, r in enumerate(REPS):
        d = sec6[s][i]
        rows.append({"System": s, "Replicate": r,
                      "D_early_0-20ns_1e-5cm2_per_s": d["solute_diffusion"]["D_early"]["D_1e-5cm2_per_s"],
                      "D_late_50-100ns_1e-5cm2_per_s": d["solute_diffusion"]["D_late"]["D_1e-5cm2_per_s"]})
df = pd.DataFrame(rows)
with_summary(df, ["D_early_0-20ns_1e-5cm2_per_s", "D_late_50-100ns_1e-5cm2_per_s"]).to_excel(
    writer, sheet_name="08_Diffusion_D_early_late", index=False)

rows = []
for s in SYS:
    for i, r in enumerate(REPS):
        d = sec6[s][i]
        rows.append({"System": s, "Replicate": r,
                      "Escape_rate_per_ns_per_mol": d["cluster_exchange"]["escape_rate_per_ns_per_molecule"],
                      "Escape_events_total": d["cluster_exchange"]["escape_events_total"]})
df = pd.DataFrame(rows)
with_summary(df, ["Escape_rate_per_ns_per_mol", "Escape_events_total"]).to_excel(
    writer, sheet_name="09_Escape_rate", index=False)

rows = []
for s in SYS:
    for i, r in enumerate(REPS):
        d = sec6[s][i]
        rows.append({"System": s, "Replicate": r,
                      "Mean_residence_time_ns": d["cluster_exchange"]["mean_residence_time_ns"],
                      "N_residence_stretches": d["cluster_exchange"]["n_residence_stretches"]})
df = pd.DataFrame(rows)
with_summary(df, ["Mean_residence_time_ns", "N_residence_stretches"]).to_excel(
    writer, sheet_name="10_Residence_time", index=False)

rows = []
for s in SYS:
    for i, r in enumerate(REPS):
        d = sec7[s][i]
        rows.append({"System": s, "Replicate": r,
                      "Solute_solute_Hbonds_per_molecule": d["hbonds"]["solute_solute_per_molecule_mean"],
                      "Solute_water_Hbonds_per_molecule": d["hbonds"]["solute_water_per_molecule_mean"]})
df = pd.DataFrame(rows)
with_summary(df, ["Solute_solute_Hbonds_per_molecule", "Solute_water_Hbonds_per_molecule"]).to_excel(
    writer, sheet_name="11_Hbonds", index=False)

rows = []
for s in SYS:
    for i, r in enumerate(REPS):
        d = sec7[s][i]
        rows.append({"System": s, "Replicate": r,
                      "Pct_molecules_stacked": d["pi_stacking"]["pct_molecules_stacked_mean"],
                      "Centroid_centroid_dist_nm": d["pi_stacking"]["centroid_centroid_dist_nm_mean"],
                      "Interplanar_angle_deg": d["pi_stacking"]["interplanar_angle_deg_mean"]})
df = pd.DataFrame(rows)
with_summary(df, ["Pct_molecules_stacked", "Centroid_centroid_dist_nm", "Interplanar_angle_deg"]).to_excel(
    writer, sheet_name="12_Pi_stacking", index=False)

writer.close()
print("saved", OUT)

import openpyxl
wb = openpyxl.load_workbook(OUT)
print("N sheets:", len(wb.sheetnames))
print(wb.sheetnames)
