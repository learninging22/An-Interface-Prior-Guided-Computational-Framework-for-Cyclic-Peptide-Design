#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch sequence design of protein structures using LigandMPNN.

The script processes PDB files from an input directory and stores the
LigandMPNN results in separate output directories. It supports both
sequential execution and Slurm array jobs through the ``--task-id``
argument.

Existing design results are detected automatically and skipped.

Example
-------
Run all PDB files sequentially:

    python run_ligandmpnn_batch.py

Run one PDB file in a Slurm array job:

    python run_ligandmpnn_batch.py --task-id 0

Specify custom input and output directories:

    python run_ligandmpnn_batch.py \
        --input-dir /path/to/input \
        --output-dir /path/to/output
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import Sequence


DEFAULT_INPUT_DIR = Path(
    "/public/home/jiangjunhao/MD/zml/relax_outputs/round_4"
)
DEFAULT_OUTPUT_DIR = Path(
    "/public/home/jiangjunhao/MD/zml/ligandmpnn_outputs/round_4"
)
DEFAULT_LIGANDMPNN_SCRIPT = Path(
    "/public/home/jiangjunhao/MD/zml/LigandMPNN-main/run.py"
)
DEFAULT_MODEL_CHECKPOINT = Path(
    "/public/home/jiangjunhao/MD/zml/"
    "LigandMPNN-main/model_params/proteinmpnn_v_48_020.pt"
)

LOGGER = logging.getLogger(__name__)


def configure_logging() -> None:
    """Configure the logging format used by the script."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def find_pdb_files(directory: Path, recursive: bool = True) -> list[Path]:
    """
    Find PDB files in a directory.

    Parameters
    ----------
    directory
        Directory to search.
    recursive
        Whether to search all subdirectories recursively.

    Returns
    -------
    list of pathlib.Path
        Sorted paths to the detected PDB files.
    """
    if not directory.is_dir():
        return []

    pattern = "**/*.pdb" if recursive else "*.pdb"

    return sorted(
        path
        for path in directory.glob(pattern)
        if path.is_file()
    )


def design_results_exist(output_directory: Path) -> bool:
    """
    Determine whether LigandMPNN output files already exist.

    The function first checks the conventional ``backbones`` directory.
    It then searches the complete output directory recursively for PDB files.

    Parameters
    ----------
    output_directory
        Output directory associated with one input structure.

    Returns
    -------
    bool
        ``True`` when at least one output PDB file is present.
    """
    backbones_directory = output_directory / "backbones"

    if find_pdb_files(backbones_directory, recursive=False):
        return True

    return bool(find_pdb_files(output_directory, recursive=True))


def build_ligandmpnn_command(
    input_pdb: Path,
    output_directory: Path,
    ligandmpnn_script: Path,
    model_checkpoint: Path,
    chains_to_design: str,
    fixed_residues: str,
    num_designs: int,
    temperature: float,
) -> list[str]:
    """
    Construct the LigandMPNN command-line invocation.

    Parameters
    ----------
    input_pdb
        Input PDB structure.
    output_directory
        Directory in which LigandMPNN results will be written.
    ligandmpnn_script
        Path to the LigandMPNN ``run.py`` script.
    model_checkpoint
        Path to the ProteinMPNN checkpoint.
    chains_to_design
        Chain identifiers to redesign.
    fixed_residues
        Residues that must remain fixed during sequence design.
    num_designs
        Number of design batches to generate.
    temperature
        Sampling temperature used during sequence generation.

    Returns
    -------
    list of str
        Command and arguments suitable for ``subprocess.run``.
    """
    return [
        sys.executable,
        str(ligandmpnn_script),
        "--pdb_path",
        str(input_pdb),
        "--chains_to_design",
        chains_to_design,
        "--out_folder",
        str(output_directory),
        "--omit_AA",
        "C",
        "--fixed_residues",
        fixed_residues,
        "--temperature",
        str(temperature),
        "--batch_size",
        "1",
        "--number_of_batches",
        str(num_designs),
        "--pack_side_chains",
        "0",
        "--number_of_packs_per_design",
        "0",
        "--pack_with_ligand_context",
        "0",
        "--checkpoint_protein_mpnn",
        str(model_checkpoint),
    ]


def run_ligandmpnn(
    input_pdb: Path,
    output_root: Path,
    ligandmpnn_script: Path,
    model_checkpoint: Path,
    chains_to_design: str = "H",
    fixed_residues: str = "H:94,97,99",
    num_designs: int = 1,
    temperature: float = 0.0001,
    overwrite: bool = False,
) -> bool:
    """
    Run LigandMPNN for a single PDB structure.

    Parameters
    ----------
    input_pdb
        Path to the input PDB file.
    output_root
        Root directory for all LigandMPNN outputs.
    ligandmpnn_script
        Path to the LigandMPNN execution script.
    model_checkpoint
        Path to the ProteinMPNN model checkpoint.
    chains_to_design
        Chain identifiers to redesign.
    fixed_residues
        Residues excluded from redesign.
    num_designs
        Number of sequence-design batches.
    temperature
        Sampling temperature.
    overwrite
        Whether to run LigandMPNN even when output PDB files already exist.

    Returns
    -------
    bool
        ``True`` if the calculation completed successfully or was skipped
        because results already existed; otherwise ``False``.
    """
    structure_name = input_pdb.stem
    structure_output_directory = output_root / structure_name
    structure_output_directory.mkdir(parents=True, exist_ok=True)

    if not overwrite and design_results_exist(structure_output_directory):
        LOGGER.info(
            "Skipping %s because design results already exist.",
            structure_name,
        )
        return True

    command = build_ligandmpnn_command(
        input_pdb=input_pdb,
        output_directory=structure_output_directory,
        ligandmpnn_script=ligandmpnn_script,
        model_checkpoint=model_checkpoint,
        chains_to_design=chains_to_design,
        fixed_residues=fixed_residues,
        num_designs=num_designs,
        temperature=temperature,
    )

    LOGGER.info("Starting LigandMPNN design for %s.", structure_name)

    try:
        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as error:
        LOGGER.error(
            "LigandMPNN could not be started for %s: %s",
            structure_name,
            error,
        )
        return False

    if result.returncode != 0:
        LOGGER.error(
            "LigandMPNN failed for %s with exit code %d.",
            structure_name,
            result.returncode,
        )

        if result.stdout.strip():
            LOGGER.error("Standard output:\n%s", result.stdout.strip())

        if result.stderr.strip():
            LOGGER.error("Standard error:\n%s", result.stderr.strip())

        return False

    LOGGER.info(
        "LigandMPNN design completed successfully for %s.",
        structure_name,
    )
    return True


def validate_configuration(
    input_directory: Path,
    ligandmpnn_script: Path,
    model_checkpoint: Path,
) -> None:
    """
    Validate all required input paths.

    Raises
    ------
    FileNotFoundError
        If a required directory or file cannot be found.
    """
    if not input_directory.is_dir():
        raise FileNotFoundError(
            f"Input directory does not exist: {input_directory}"
        )

    if not ligandmpnn_script.is_file():
        raise FileNotFoundError(
            f"LigandMPNN script does not exist: {ligandmpnn_script}"
        )

    if not model_checkpoint.is_file():
        raise FileNotFoundError(
            f"Model checkpoint does not exist: {model_checkpoint}"
        )


def parse_arguments(
    arguments: Sequence[str] | None = None,
) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Parameters
    ----------
    arguments
        Optional argument sequence. When omitted, arguments are read from
        ``sys.argv``.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Run LigandMPNN sequence design for a collection of PDB files."
        )
    )

    parser.add_argument(
        "--task-id",
        "--task_id",
        dest="task_id",
        type=int,
        default=-1,
        help=(
            "Zero-based input index for a Slurm array job. "
            "A negative value processes all input files sequentially."
        ),
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing input PDB files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Root directory for LigandMPNN results.",
    )
    parser.add_argument(
        "--ligandmpnn-script",
        type=Path,
        default=DEFAULT_LIGANDMPNN_SCRIPT,
        help="Path to the LigandMPNN run.py script.",
    )
    parser.add_argument(
        "--model-checkpoint",
        type=Path,
        default=DEFAULT_MODEL_CHECKPOINT,
        help="Path to the ProteinMPNN model checkpoint.",
    )
    parser.add_argument(
        "--chains-to-design",
        default="H",
        help="Chain identifiers to redesign.",
    )
    parser.add_argument(
        "--fixed-residues",
        default="H:94,97,99",
        help="Residues to keep fixed during sequence design.",
    )
    parser.add_argument(
        "--num-designs",
        type=int,
        default=1,
        help="Number of design batches generated for each structure.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=0.0001,
        help="Sampling temperature used by LigandMPNN.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Recalculate structures even when output PDB files exist.",
    )

    parsed_arguments = parser.parse_args(arguments)

    if parsed_arguments.num_designs < 1:
        parser.error("--num-designs must be greater than or equal to 1.")

    if parsed_arguments.temperature < 0:
        parser.error("--temperature must be non-negative.")

    return parsed_arguments


def main(arguments: Sequence[str] | None = None) -> int:
    """
    Execute the LigandMPNN batch-design workflow.

    Parameters
    ----------
    arguments
        Optional command-line argument sequence.

    Returns
    -------
    int
        Process exit status. A value of zero indicates successful execution.
    """
    configure_logging()
    args = parse_arguments(arguments)

    try:
        validate_configuration(
            input_directory=args.input_dir,
            ligandmpnn_script=args.ligandmpnn_script,
            model_checkpoint=args.model_checkpoint,
        )
    except FileNotFoundError as error:
        LOGGER.error("%s", error)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    input_pdb_files = find_pdb_files(
        args.input_dir,
        recursive=False,
    )

    LOGGER.info(
        "Detected %d input PDB file(s) in %s.",
        len(input_pdb_files),
        args.input_dir,
    )

    if not input_pdb_files:
        LOGGER.warning("No input PDB files were found.")
        return 0

    if args.task_id >= 0:
        if args.task_id >= len(input_pdb_files):
            LOGGER.warning(
                "Task ID %d is outside the valid range 0-%d; "
                "the task will be skipped.",
                args.task_id,
                len(input_pdb_files) - 1,
            )
            return 0

        selected_pdb_files = [input_pdb_files[args.task_id]]
    else:
        selected_pdb_files = input_pdb_files

    successful_runs = 0

    for input_pdb in selected_pdb_files:
        completed = run_ligandmpnn(
            input_pdb=input_pdb,
            output_root=args.output_dir,
            ligandmpnn_script=args.ligandmpnn_script,
            model_checkpoint=args.model_checkpoint,
            chains_to_design=args.chains_to_design,
            fixed_residues=args.fixed_residues,
            num_designs=args.num_designs,
            temperature=args.temperature,
            overwrite=args.overwrite,
        )

        successful_runs += int(completed)

    failed_runs = len(selected_pdb_files) - successful_runs

    LOGGER.info(
        "Processing finished: %d successful or skipped, %d failed.",
        successful_runs,
        failed_runs,
    )

    return 0 if failed_runs == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())