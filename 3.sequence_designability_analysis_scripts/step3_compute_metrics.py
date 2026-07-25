#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 3: Metric Calculation Script
Calculate scaffold-level designability metrics from LigandMPNN analysis results.
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

# 20 standard amino acids
AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AROMATIC = list("FYW")
CHARGED = list("KRHDE")
POLAR = list("STNQ")
GP = ["G", "P"]


# ============== Helper Functions ==============
def detect_length_col(df: pd.DataFrame) -> str:
    """Detect the column containing sequence or ring length."""
    for c in ["ring_len", "subdir_folder", "job_folder", "length"]:
        if c in df.columns:
            return c

    raise ValueError(
        "Length column not found. Expected one of: "
        "ring_len / subdir_folder / job_folder / length"
    )


def bootstrap_median_ci(values, n_boot=10000, seed=0):
    """Calculate a bootstrap confidence interval for the median."""
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]

    if len(arr) == 0:
        return np.nan, np.nan, np.nan

    rng = np.random.default_rng(seed)
    point = float(np.median(arr))
    boots = []

    for _ in range(n_boot):
        sample = rng.choice(arr, size=len(arr), replace=True)
        boots.append(np.median(sample))

    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def q10(x):
    """Calculate the 10th percentile."""
    return float(np.quantile(np.asarray(x, dtype=float), 0.10))


def ensure_aa_cols(wide: pd.DataFrame) -> pd.DataFrame:
    """Ensure that columns for all standard amino acids are present."""
    for aa in AA20:
        col = f"p_{aa}"

        if col not in wide.columns:
            wide[col] = 0.0

    return wide


# ============== Main Function ==============
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Calculate designability metrics from LigandMPNN analysis results."
        )
    )

    # Input and output
    ap.add_argument(
        "--analysis_dir",
        required=True,
        help=(
            "Directory containing all_samples.csv and "
            "position_frequencies.csv."
        ),
    )
    ap.add_argument(
        "--outdir",
        default="step3_metrics_out",
        help="Output directory for calculated metrics.",
    )

    # Threshold parameters
    ap.add_argument(
        "--locked_neff_cutoff",
        type=float,
        default=1.5,
        help="Neff cutoff for identifying locked positions. Default: 1.5.",
    )
    ap.add_argument(
        "--locked_top1_cutoff",
        type=float,
        default=0.8,
        help="Top-1 frequency cutoff for locked positions. Default: 0.8.",
    )
    ap.add_argument(
        "--gp_sensitivity_cutoff",
        type=float,
        default=0.6,
        help="G/P enrichment sensitivity cutoff. Default: 0.6.",
    )
    ap.add_argument(
        "--bootstrap_n",
        type=int,
        default=10000,
        help="Number of bootstrap iterations. Default: 10000.",
    )

    args = ap.parse_args()

    analysis_dir = Path(args.analysis_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Validate input files
    all_samples_path = analysis_dir / "all_samples.csv"
    pos_freq_path = analysis_dir / "position_frequencies.csv"

    if not all_samples_path.exists():
        raise SystemExit(f"File does not exist: {all_samples_path}")

    if not pos_freq_path.exists():
        raise SystemExit(f"File does not exist: {pos_freq_path}")

    # Read input data
    all_df = pd.read_csv(all_samples_path)
    pos_df = pd.read_csv(pos_freq_path)

    # Standardize the ring-length column
    length_col_all = detect_length_col(all_df)
    length_col_pos = detect_length_col(pos_df)

    all_df["ring_len"] = pd.to_numeric(
        all_df[length_col_all],
        errors="coerce",
    )
    pos_df["ring_len"] = pd.to_numeric(
        pos_df[length_col_pos],
        errors="coerce",
    )

    print(f"\n{'=' * 60}")
    print("Computing Metrics")
    print(f"{'=' * 60}")
    print(f"Locked Neff cutoff: {args.locked_neff_cutoff}")
    print(f"Locked Top1 cutoff: {args.locked_top1_cutoff}")
    print(f"GP sensitivity cutoff: {args.gp_sensitivity_cutoff}")
    print(f"Bootstrap iterations: {args.bootstrap_n}")
    print(f"{'=' * 60}\n")

    # ============== 1) Scaffold-Level Global Designability ==============
    global_df = (
        all_df.groupby(
            ["ring_len", "scaffold_name"],
            as_index=False,
        )
        .agg(
            n_samples=("designed_seq", "size"),
            median_nll_per_res=("nll_per_res", "median"),
            p10_nll_per_res=("nll_per_res", q10),
            median_avg_log_prob=("avg_log_prob", "median"),
            p10_avg_log_prob=("avg_log_prob", q10),
            best_avg_log_prob=("avg_log_prob", "max"),
            worst_avg_log_prob=("avg_log_prob", "min"),
            mean_seq_rec=("seq_rec", "mean"),
        )
        .sort_values(["ring_len", "scaffold_name"])
        .reset_index(drop=True)
    )

    global_df.to_csv(
        outdir / "per_scaffold_global.csv",
        index=False,
    )
    print(f"Wrote: {outdir / 'per_scaffold_global.csv'}")

    # ============== 2) Position-Level Metrics ==============
    idx_cols = [
        "ring_len",
        "scaffold_name",
        "position_1based",
    ]

    meta_df = (
        pos_df.groupby(
            idx_cols,
            as_index=False,
        )
        .agg(
            top_aa=("top_aa", "first"),
            top_freq=("top_freq", "first"),
            entropy_nats=("entropy_nats", "first"),
            K=("K", "first"),
        )
    )

    wide = pos_df.pivot_table(
        index=idx_cols,
        columns="aa",
        values="freq",
        aggfunc="first",
        fill_value=0.0,
    ).reset_index()

    # Rename amino-acid columns
    wide.columns = [
        f"p_{c}" if c in AA20 else c
        for c in wide.columns
    ]

    wide = ensure_aa_cols(wide)

    per_pos = meta_df.merge(
        wide,
        on=idx_cols,
        how="left",
    )

    # Calculate derived metrics
    per_pos["Neff"] = np.exp(per_pos["entropy_nats"])

    per_pos["locked_position"] = (
        (per_pos["Neff"] < args.locked_neff_cutoff)
        | (per_pos["top_freq"] > args.locked_top1_cutoff)
    )

    per_pos["p_GP"] = per_pos["p_G"] + per_pos["p_P"]

    per_pos["GP_locked"] = (
        per_pos["top_aa"].isin(GP)
        & (per_pos["top_freq"] > args.locked_top1_cutoff)
    )

    per_pos["GP_enriched_sens"] = (
        per_pos["p_GP"] > args.gp_sensitivity_cutoff
    )

    # Chemical-space coverage
    per_pos["aromatic_cov"] = per_pos[
        [f"p_{aa}" for aa in AROMATIC]
    ].sum(axis=1)

    per_pos["charged_cov"] = per_pos[
        [f"p_{aa}" for aa in CHARGED]
    ].sum(axis=1)

    per_pos["polar_cov"] = per_pos[
        [f"p_{aa}" for aa in POLAR]
    ].sum(axis=1)

    per_pos["pharma_cov"] = (
        per_pos["aromatic_cov"]
        + per_pos["charged_cov"]
        + per_pos["polar_cov"]
    )

    per_pos = per_pos.sort_values(
        [
            "ring_len",
            "scaffold_name",
            "position_1based",
        ]
    ).reset_index(drop=True)

    per_pos.to_csv(
        outdir / "per_position_metrics.csv",
        index=False,
    )
    print(f"Wrote: {outdir / 'per_position_metrics.csv'}")

    # ============== 3) Scaffold-Level Position Summary ==============
    scaffold_rows = []

    for (ring_len, scaffold_name), sub in per_pos.groupby(
        ["ring_len", "scaffold_name"]
    ):
        nonlocked = sub.loc[
            ~sub["locked_position"]
        ].copy()

        scaffold_rows.append(
            {
                "ring_len": ring_len,
                "scaffold_name": scaffold_name,
                "n_positions": int(
                    sub["position_1based"].nunique()
                ),
                "locked_positions": int(
                    sub["locked_position"].sum()
                ),
                "locked_fraction": float(
                    sub["locked_position"].mean()
                ),
                "gp_locked_positions": int(
                    sub["GP_locked"].sum()
                ),
                "gp_enriched_positions": int(
                    sub["GP_enriched_sens"].sum()
                ),
                "mean_Neff": float(
                    sub["Neff"].mean()
                ),
                "median_Neff": float(
                    sub["Neff"].median()
                ),
                "min_Neff": float(
                    sub["Neff"].min()
                ),
                "max_Neff": float(
                    sub["Neff"].max()
                ),
                "mean_top_freq": float(
                    sub["top_freq"].mean()
                ),
                "max_top_freq": float(
                    sub["top_freq"].max()
                ),
                "mean_p_G": float(
                    sub["p_G"].mean()
                ),
                "mean_p_P": float(
                    sub["p_P"].mean()
                ),
                "mean_p_GP": float(
                    sub["p_GP"].mean()
                ),
                "mean_aromatic_cov": float(
                    sub["aromatic_cov"].mean()
                ),
                "mean_charged_cov": float(
                    sub["charged_cov"].mean()
                ),
                "mean_polar_cov": float(
                    sub["polar_cov"].mean()
                ),
                "mean_pharma_cov": float(
                    sub["pharma_cov"].mean()
                ),

                # Statistics for non-locked positions
                "n_nonlocked": int(
                    (~sub["locked_position"]).sum()
                ),
                "nonlocked_mean_Neff": (
                    float(nonlocked["Neff"].mean())
                    if len(nonlocked)
                    else np.nan
                ),
                "nonlocked_mean_top_freq": (
                    float(nonlocked["top_freq"].mean())
                    if len(nonlocked)
                    else np.nan
                ),
                "nonlocked_mean_aromatic_cov": (
                    float(nonlocked["aromatic_cov"].mean())
                    if len(nonlocked)
                    else np.nan
                ),
                "nonlocked_mean_charged_cov": (
                    float(nonlocked["charged_cov"].mean())
                    if len(nonlocked)
                    else np.nan
                ),
                "nonlocked_mean_polar_cov": (
                    float(nonlocked["polar_cov"].mean())
                    if len(nonlocked)
                    else np.nan
                ),
                "nonlocked_mean_pharma_cov": (
                    float(nonlocked["pharma_cov"].mean())
                    if len(nonlocked)
                    else np.nan
                ),
            }
        )

    scaffold_pos_df = (
        pd.DataFrame(scaffold_rows)
        .sort_values(
            ["ring_len", "scaffold_name"]
        )
        .reset_index(drop=True)
    )

    scaffold_pos_df.to_csv(
        outdir / "per_scaffold_position_summary.csv",
        index=False,
    )
    print(
        f"Wrote: "
        f"{outdir / 'per_scaffold_position_summary.csv'}"
    )

    # ============== 4) Final Scaffold-Level Metrics Table ==============
    per_scaffold = (
        global_df.merge(
            scaffold_pos_df,
            on=["ring_len", "scaffold_name"],
            how="inner",
        )
        .sort_values(
            ["ring_len", "scaffold_name"]
        )
        .reset_index(drop=True)
    )

    per_scaffold.to_csv(
        outdir / "per_scaffold_metrics.csv",
        index=False,
    )
    print(f"Wrote: {outdir / 'per_scaffold_metrics.csv'}")

    # ============== 5) Ring-Length Summary and Bootstrap CI ==============
    summary_metrics = [
        "median_avg_log_prob",
        "p10_avg_log_prob",
        "median_nll_per_res",
        "p10_nll_per_res",
        "locked_positions",
        "gp_locked_positions",
        "locked_fraction",
        "mean_Neff",
        "median_Neff",
        "nonlocked_mean_Neff",
        "mean_pharma_cov",
        "nonlocked_mean_pharma_cov",
        "mean_seq_rec",
    ]

    length_rows = []

    for ring_len, sub in per_scaffold.groupby("ring_len"):
        base = {
            "ring_len": ring_len,
            "n_scaffolds": len(sub),
        }

        # Proportion-based metrics
        base["frac_scaffolds_locked_le_2"] = float(
            (sub["locked_positions"] <= 2).mean()
        )
        base["frac_scaffolds_locked_le_3"] = float(
            (sub["locked_positions"] <= 3).mean()
        )
        base["frac_scaffolds_gp_locked_le_1"] = float(
            (sub["gp_locked_positions"] <= 1).mean()
        )
        base["frac_scaffolds_gp_locked_le_2"] = float(
            (sub["gp_locked_positions"] <= 2).mean()
        )

        length_rows.append(base)

        for metric in summary_metrics:
            if metric in sub.columns:
                point, lo, hi = bootstrap_median_ci(
                    sub[metric].values,
                    n_boot=args.bootstrap_n,
                    seed=0,
                )

                length_rows.append(
                    {
                        "ring_len": ring_len,
                        "metric": metric,
                        "median": point,
                        "ci_low": lo,
                        "ci_high": hi,
                        "n_scaffolds": len(sub),
                    }
                )

    # Separate base rows and metric rows
    length_summary_long = (
        pd.DataFrame(
            [
                x
                for x in length_rows
                if "metric" in x
            ]
        )
        .sort_values(
            ["ring_len", "metric"]
        )
        .reset_index(drop=True)
    )

    length_summary_long.to_csv(
        outdir / "per_length_summary.csv",
        index=False,
    )
    print(f"Wrote: {outdir / 'per_length_summary.csv'}")

    length_summary_basic = (
        pd.DataFrame(
            [
                x
                for x in length_rows
                if "metric" not in x
            ]
        )
        .sort_values(["ring_len"])
        .reset_index(drop=True)
    )

    length_summary_basic.to_csv(
        outdir / "per_length_basic_flags.csv",
        index=False,
    )
    print(
        f"Wrote: "
        f"{outdir / 'per_length_basic_flags.csv'}"
    )

    # ============== Completion Summary ==============
    print(f"\n{'=' * 60}")
    print("Metrics Computation Complete!")
    print(f"{'=' * 60}")
    print(f"Output directory: {outdir.resolve()}")
    print(f"Scaffolds analyzed: {len(per_scaffold)}")
    print(
        f"Ring lengths: "
        f"{sorted(per_scaffold['ring_len'].unique())}"
    )
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()