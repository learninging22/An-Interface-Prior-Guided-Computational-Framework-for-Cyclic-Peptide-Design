#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1: LigandMPNN Design Script
Perform sequence design on representative PDB structures with different
ring lengths (6-16) using a temperature of T=0.0001.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import time
from pathlib import Path
from typing import List, Optional


# ============== Default Path Configuration ==============
DEFAULT_RUN_PY = "/public/home/jiangjunhao/MD/zml/LigandMPNN-main/run.py"
DEFAULT_CKPT = "/public/home/jiangjunhao/MD/zml/LigandMPNN-main/model_params/proteinmpnn_v_48_020.pt"


# ============== Utility Functions ==============
def format_seconds(sec: float) -> str:
    """Format seconds as a human-readable time string."""
    sec = int(sec)
    h = sec // 3600
    m = (sec % 3600) // 60
    s = sec % 60

    if h > 0:
        return f"{h}h{m:02d}m{s:02d}s"
    if m > 0:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def progress_line(
    done: int,
    total: int,
    start_ts: float,
    prefix: str = "",
) -> str:
    """Return a formatted progress bar string."""
    now = time.time()
    elapsed = now - start_ts
    pct = int(100 * done / total) if total else 100

    if elapsed > 0 and done > 0:
        speed = done / (elapsed / 60.0)
        remain = total - done
        eta_sec = remain / (done / elapsed)
    else:
        speed = 0.0
        eta_sec = 0.0

    bar_len = 36
    filled = int(bar_len * done / total) if total else bar_len
    bar = "#" * filled + " " * (bar_len - filled)

    return (
        f"{prefix}[{bar}] {pct:3d}% ({done}/{total}) | "
        f"{format_seconds(elapsed)} | {speed:.2f} jobs/min | "
        f"ETA {format_seconds(eta_sec)}"
    )


# ============== LigandMPNN Execution Functions ==============
def run_cmd(
    cmd: List[str],
    log_path: Path,
    env: Optional[dict] = None,
) -> int:
    """Execute a command and write its output to a log file."""
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w") as lf:
        lf.write("CMD:\n")
        lf.write(" ".join(shlex.quote(c) for c in cmd) + "\n\n")
        lf.flush()

        process = subprocess.run(
            cmd,
            stdout=lf,
            stderr=lf,
            env=env,
        )
        return process.returncode


def supports_seed_flag(run_py: str) -> bool:
    """Check whether run.py supports the --seed argument."""
    try:
        process = subprocess.run(
            ["python", run_py, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        return "--seed" in process.stdout.lower()
    except Exception:
        return False


def list_pdbs_in_folder(folder: Path, pattern: str) -> List[Path]:
    """List all PDB files in a folder that match the specified pattern."""
    if not folder.is_dir():
        return []

    return sorted(
        [
            pdb_path
            for pdb_path in folder.glob(pattern)
            if pdb_path.is_file()
        ]
    )


def run_ligandmpnn_one_job(
    run_py: str,
    ckpt: str,
    pdb_path: Path,
    chains_to_design: str,
    out_folder: Path,
    temperature: float,
    seed: int,
    num_batches: int,
    use_seed_flag: bool,
    pack_side_chains: int = 0,
    number_of_packs_per_design: int = 0,
    pack_with_ligand_context: int = 0,
) -> int:
    """
    Run a single LigandMPNN design job.

    Args:
        run_py: Path to run.py.
        ckpt: Path to the model checkpoint.
        pdb_path: Path to the input PDB file.
        chains_to_design: Chain or chains to design.
        out_folder: Output directory.
        temperature: Sampling temperature.
        seed: Random seed.
        num_batches: Number of sequences generated for each seed.
        use_seed_flag: Whether to pass the --seed argument.
        pack_side_chains: Side-chain packing setting.
        number_of_packs_per_design: Number of packing runs per design.
        pack_with_ligand_context: Whether to use ligand context during packing.

    Returns:
        Process return code. A return code of 0 indicates success.
    """
    out_folder.mkdir(parents=True, exist_ok=True)
    log_path = out_folder / f"{pdb_path.stem}.log"

    cmd = [
        "python",
        run_py,
        "--pdb_path",
        str(pdb_path),
        "--chains_to_design",
        chains_to_design,
        "--out_folder",
        str(out_folder),
        "--temperature",
        str(temperature),
        "--batch_size",
        "1",
        "--number_of_batches",
        str(num_batches),
        "--pack_side_chains",
        str(pack_side_chains),
        "--number_of_packs_per_design",
        str(number_of_packs_per_design),
        "--pack_with_ligand_context",
        str(pack_with_ligand_context),
        "--checkpoint_protein_mpnn",
        ckpt,
    ]

    env = os.environ.copy()

    if use_seed_flag:
        cmd += ["--seed", str(seed)]
    else:
        env["PYTHONHASHSEED"] = str(seed)

    return run_cmd(
        cmd,
        log_path=log_path,
        env=env,
    )


# ============== Main Function ==============
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Use LigandMPNN at T=0.0001 to perform sequence design on "
            "representative PDB structures with different ring lengths."
        )
    )

    # Path arguments
    parser.add_argument(
        "--run_py",
        default=DEFAULT_RUN_PY,
        help="Path to the LigandMPNN run.py script.",
    )
    parser.add_argument(
        "--ckpt",
        default=DEFAULT_CKPT,
        help="Path to the model checkpoint.",
    )
    parser.add_argument(
        "--pdb_root",
        default=".",
        help="Root directory containing the ring-length folders 6-16.",
    )
    parser.add_argument(
        "--out_root",
        default="mpnn_out_T0p0001",
        help="Root output directory.",
    )

    # Design arguments
    parser.add_argument(
        "--chains_to_design",
        required=True,
        help="Chain or chains to design, such as 'H' or 'A'.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=0.0001,
        help="Sampling temperature. Default: 0.0001.",
    )
    parser.add_argument(
        "--seeds",
        default="0,1,2",
        help="Comma-separated random seed list. Default: '0,1,2'.",
    )
    parser.add_argument(
        "--per_seed",
        type=int,
        default=50,
        help="Number of sequences generated per seed. Default: 50.",
    )

    # Ring-length configuration
    parser.add_argument(
        "--ring_lengths",
        default="6,7,8,9,10,11,12,13,14,15,16",
        help=(
            "Comma-separated ring-length list. "
            "Default: '6,7,8,9,10,11,12,13,14,15,16'."
        ),
    )

    # Optional arguments
    parser.add_argument(
        "--pattern",
        default="*.pdb",
        help="PDB filename pattern. Default: '*.pdb'.",
    )
    parser.add_argument(
        "--pack_side_chains",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--number_of_packs_per_design",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--pack_with_ligand_context",
        type=int,
        default=0,
    )

    args = parser.parse_args()

    # Validate paths
    run_py = Path(args.run_py)
    ckpt = Path(args.ckpt)

    if not run_py.exists():
        raise FileNotFoundError(f"run.py not found: {run_py}")

    if not ckpt.exists():
        raise FileNotFoundError(f"checkpoint not found: {ckpt}")

    pdb_root = Path(args.pdb_root)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    # Parse arguments
    seeds = [
        int(value.strip())
        for value in args.seeds.split(",")
        if value.strip()
    ]
    ring_lengths = [
        int(value.strip())
        for value in args.ring_lengths.split(",")
        if value.strip()
    ]

    if not seeds:
        raise SystemExit("No seeds specified.")

    if not ring_lengths:
        raise SystemExit("No ring lengths specified.")

    # Check seed support
    seed_supported = supports_seed_flag(str(run_py))

    if not seed_supported:
        print("[WARN] LigandMPNN run.py does NOT show '--seed' in --help.")
        print("       Will run without --seed, with limited reproducibility.")
        print(
            "       Recommended fix: add --seed to run.py and set the "
            "random, NumPy, and PyTorch seeds.\n"
        )

    # ============== Collect All PDB Files ==============
    all_pdbs = {}
    total_pdbs = 0

    for ring_len in ring_lengths:
        ring_folder = pdb_root / str(ring_len)
        pdbs = list_pdbs_in_folder(
            ring_folder,
            args.pattern,
        )

        if pdbs:
            all_pdbs[ring_len] = pdbs
            total_pdbs += len(pdbs)

            print(
                f"Ring length {ring_len}: found {len(pdbs)} "
                f"PDBs in {ring_folder}"
            )
        else:
            print(
                f"[WARN] No PDBs found in {ring_folder} "
                f"with pattern {args.pattern}"
            )

    if total_pdbs == 0:
        raise SystemExit(
            "No PDBs found in any ring-length folder."
        )

    # ============== Calculate the Total Number of Jobs ==============
    jobs_total = total_pdbs * len(seeds)
    done = 0
    failed = 0
    start = time.time()

    # ============== Print Configuration ==============
    print("\n" + "=" * 60)
    print("LigandMPNN Design Configuration")
    print("=" * 60)
    print(f"Temperature: {args.temperature}")
    print(f"Seeds: {seeds}")
    print(f"Sequences per seed: {args.per_seed}")
    print(f"Total sequences per PDB: {len(seeds) * args.per_seed}")
    print(f"Total jobs: {jobs_total}")
    print(f"Output root: {out_root.resolve()}")
    print("=" * 60 + "\n")

    # ============== Run Design Jobs ==============
    for ring_len, pdbs in all_pdbs.items():
        print(f"\n>>> Processing ring length: {ring_len}")

        for pdb in pdbs:
            for seed in seeds:
                # Output structure:
                # out_root / ring_length / pdb_name / seed_number / ...
                out_folder = (
                    out_root
                    / str(ring_len)
                    / pdb.stem
                    / f"seed{seed}"
                )

                return_code = run_ligandmpnn_one_job(
                    run_py=str(run_py),
                    ckpt=str(ckpt),
                    pdb_path=pdb,
                    chains_to_design=args.chains_to_design,
                    out_folder=out_folder,
                    temperature=args.temperature,
                    seed=seed,
                    num_batches=args.per_seed,
                    use_seed_flag=seed_supported,
                    pack_side_chains=args.pack_side_chains,
                    number_of_packs_per_design=(
                        args.number_of_packs_per_design
                    ),
                    pack_with_ligand_context=(
                        args.pack_with_ligand_context
                    ),
                )

                done += 1

                if return_code != 0:
                    failed += 1
                    print(
                        f"\n[WARN] Job failed: "
                        f"{pdb.name} seed{seed}"
                    )

                print(
                    "\r"
                    + progress_line(
                        done,
                        jobs_total,
                        start,
                        prefix=f"R{ring_len} ",
                    ),
                    end="",
                    flush=True,
                )

    # ============== Completion Summary ==============
    print("\n\n" + "=" * 60)
    print("Design Complete!")
    print("=" * 60)
    print(f"Total jobs: {jobs_total}")
    print(f"Successful: {jobs_total - failed}")
    print(f"Failed: {failed}")
    print(f"Output directory: {out_root.resolve()}")
    print("=" * 60)


if __name__ == "__main__":
    main()