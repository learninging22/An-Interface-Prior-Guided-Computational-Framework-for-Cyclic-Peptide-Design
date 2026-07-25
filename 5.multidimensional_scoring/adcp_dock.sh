#!/usr/bin/env bash

# ==============================================================================
# Supplementary Script: ADCP Redocking of Cyclic Peptides
#
# Description:
#   This script performs AutoDock CrankPep (ADCP) redocking for cyclic peptide
#   structures located in the round_1, round_2, and round_3 directories.
#
#   Each peptide directory must contain one PDBQT file whose filename includes
#   the string "pep". The peptide sequence is extracted from the PDBQT file and
#   used as input for ADCP.
#
#   For each peptide, the script generates:
#       1. ADCP docking output files
#       2. A complete execution log (slurm.out)
#       3. A summary CSV file:
#          <ligand_name>_result.csv
#
# Requirements:
#   - Bash 4.0 or later
#   - AutoDock CrankPep (adcp) available in PATH
#   - A valid ADCP target file
#
# Usage:
#   bash run_adcp_cyclic_peptide_redocking.sh
# ==============================================================================

set -uo pipefail


# ------------------------------------------------------------------------------
# User-configurable parameters
# ------------------------------------------------------------------------------

TARGET_FILE="6WJL.trg"

ROUND_DIRECTORIES=(
    "round_1"
    "round_2"
    "round_3"
)

MAX_PARALLEL_JOBS=20
NUMBER_OF_MODELS=20
NUMBER_OF_STEPS=500000
NATIVE_CONTACT_WEIGHT=0.8

PDBQT_PATTERN="*pep*.pdbqt"


# ------------------------------------------------------------------------------
# Amino-acid mapping
# ------------------------------------------------------------------------------

declare -A AMINO_ACID_MAP=(
    [ALA]="a"
    [ARG]="r"
    [ASN]="n"
    [ASP]="d"
    [CYS]="c"
    [GLU]="e"
    [GLN]="q"
    [GLY]="g"
    [HIS]="h"
    [ILE]="i"
    [LEU]="l"
    [LYS]="k"
    [MET]="m"
    [PHE]="f"
    [PRO]="p"
    [SER]="s"
    [THR]="t"
    [TRP]="w"
    [TYR]="y"
    [VAL]="v"
)


# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

log_message() {
    local message="$1"
    printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$message"
}


check_dependencies() {
    if ! command -v adcp >/dev/null 2>&1; then
        printf "ERROR: The 'adcp' executable was not found in PATH.\n" >&2
        exit 1
    fi

    if [[ ! -f "$TARGET_FILE" ]]; then
        printf "ERROR: ADCP target file not found: %s\n" "$TARGET_FILE" >&2
        exit 1
    fi
}


extract_peptide_sequence() {
    local pdbqt_file="$1"
    local sequence=""
    local residue_number
    local residue_name
    local amino_acid

    while read -r residue_number residue_name; do
        residue_name="${residue_name^^}"
        amino_acid="${AMINO_ACID_MAP[$residue_name]:-x}"
        sequence+="$amino_acid"
    done < <(
        awk '/^ATOM/ {print $6, $4}' "$pdbqt_file" \
            | sort -n -k1,1 \
            | awk '!seen[$1 FS $2]++'
    )

    printf "%s" "$sequence"
}


extract_runtime() {
    local log_file="$1"

    awk '
        /Docking performed in/ {
            runtime = $4
        }
        END {
            if (runtime != "") {
                print runtime
            }
        }
    ' "$log_file"
}


extract_best_affinity() {
    local log_file="$1"

    awk '
        /^[[:space:]]*[0-9]+[[:space:]]+-[0-9]/ {
            affinity = $2

            if (!found || affinity < best) {
                best = affinity
                found = 1
            }
        }
        END {
            if (found) {
                print best
            }
        }
    ' "$log_file"
}


clean_previous_outputs() {
    local peptide_directory="$1"

    find "$peptide_directory" \
        -maxdepth 1 \
        \( \
            -name "*_redocking*" \
            -o -name "*.log" \
            -o -name "*.pdb" \
            -o -name "*.gz" \
            -o -name "slurm.out" \
            -o -name "*_result.csv" \
        \) \
        -exec rm -rf {} +
}


run_adcp_job() {
    local peptide_directory="$1"
    local pdbqt_file="$2"
    local peptide_sequence="$3"
    local target_file_absolute="$4"

    local ligand_name
    local pdbqt_basename
    local output_prefix
    local log_file
    local csv_file
    local start_time
    local end_time
    local runtime
    local best_affinity
    local return_code

    ligand_name="$(basename "$peptide_directory")"
    pdbqt_basename="$(basename "$pdbqt_file")"
    output_prefix="${ligand_name}_redocking"
    log_file="slurm.out"
    csv_file="${ligand_name}_result.csv"

    (
        cd "$peptide_directory" || exit 1

        start_time="$(date '+%Y-%m-%d %H:%M:%S')"

        {
            printf "ADCP cyclic peptide redocking\n"
            printf "Ligand name: %s\n" "$ligand_name"
            printf "Input PDBQT: %s\n" "$pdbqt_basename"
            printf "Peptide sequence: %s\n" "$peptide_sequence"
            printf "Target file: %s\n" "$target_file_absolute"
            printf "Number of models: %s\n" "$NUMBER_OF_MODELS"
            printf "Number of Monte Carlo steps: %s\n" "$NUMBER_OF_STEPS"
            printf "Native-contact weight: %s\n" "$NATIVE_CONTACT_WEIGHT"
            printf "Start time: %s\n\n" "$start_time"

            printf "Command:\n"
            printf "adcp -t %q -s %q -N %q -n %q -cyc -o %q -ref %q -nc %q\n\n" \
                "$target_file_absolute" \
                "$peptide_sequence" \
                "$NUMBER_OF_MODELS" \
                "$NUMBER_OF_STEPS" \
                "$output_prefix" \
                "$pdbqt_basename" \
                "$NATIVE_CONTACT_WEIGHT"
        } > "$log_file"

        adcp \
            -t "$target_file_absolute" \
            -s "$peptide_sequence" \
            -N "$NUMBER_OF_MODELS" \
            -n "$NUMBER_OF_STEPS" \
            -cyc \
            -o "$output_prefix" \
            -ref "$pdbqt_basename" \
            -nc "$NATIVE_CONTACT_WEIGHT" \
            >> "$log_file" 2>&1

        return_code=$?
        end_time="$(date '+%Y-%m-%d %H:%M:%S')"

        {
            printf "\nEnd time: %s\n" "$end_time"
            printf "Exit status: %s\n" "$return_code"
        } >> "$log_file"

        runtime="$(extract_runtime "$log_file")"
        best_affinity="$(extract_best_affinity "$log_file")"

        {
            printf "Peptide,Sequence,Best_Affinity_kcal_per_mol,Runtime_s,Exit_Status\n"
            printf "%s,%s,%s,%s,%s\n" \
                "$ligand_name" \
                "$peptide_sequence" \
                "${best_affinity:-NA}" \
                "${runtime:-NA}" \
                "$return_code"
        } > "$csv_file"

        if [[ "$return_code" -eq 0 ]]; then
            log_message "Completed ADCP redocking for ${ligand_name}."
        else
            log_message "ADCP redocking failed for ${ligand_name}; see ${peptide_directory}/${log_file}."
        fi

        exit "$return_code"
    )
}


# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main() {
    local target_file_absolute
    local round_directory
    local peptide_directory
    local pdbqt_file
    local peptide_sequence
    local job_count=0
    local failed_jobs=0

    check_dependencies

    target_file_absolute="$(
        cd "$(dirname "$TARGET_FILE")" &&
        printf "%s/%s" "$(pwd)" "$(basename "$TARGET_FILE")"
    )"

    log_message "Starting ADCP cyclic peptide redocking workflow."
    log_message "Target file: ${target_file_absolute}"
    log_message "Maximum parallel jobs: ${MAX_PARALLEL_JOBS}"

    for round_directory in "${ROUND_DIRECTORIES[@]}"; do
        printf "\n"
        log_message "Processing directory: ${round_directory}"

        if [[ ! -d "$round_directory" ]]; then
            log_message "Directory not found; skipping: ${round_directory}"
            continue
        fi

        mapfile -t peptide_directories < <(
            find "$round_directory" \
                -mindepth 1 \
                -maxdepth 1 \
                -type d \
                -print \
                | sort
        )

        if [[ "${#peptide_directories[@]}" -eq 0 ]]; then
            log_message "No peptide directories found in ${round_directory}."
            continue
        fi

        for peptide_directory in "${peptide_directories[@]}"; do
            pdbqt_file="$(
                find "$peptide_directory" \
                    -maxdepth 1 \
                    -type f \
                    -name "$PDBQT_PATTERN" \
                    -print \
                    | sort \
                    | head -n 1
            )"

            if [[ -z "$pdbqt_file" ]]; then
                log_message \
                    "No PDBQT file matching '${PDBQT_PATTERN}' was found in ${peptide_directory}; skipping."
                continue
            fi

            peptide_sequence="$(extract_peptide_sequence "$pdbqt_file")"

            if [[ -z "$peptide_sequence" ]]; then
                log_message \
                    "The peptide sequence could not be extracted from ${pdbqt_file}; skipping."
                continue
            fi

            log_message \
                "Submitting ${peptide_directory} with sequence ${peptide_sequence}."

            clean_previous_outputs "$peptide_directory"

            run_adcp_job \
                "$peptide_directory" \
                "$pdbqt_file" \
                "$peptide_sequence" \
                "$target_file_absolute" &

            ((job_count++))

            if (( job_count % MAX_PARALLEL_JOBS == 0 )); then
                log_message \
                    "Waiting for the current batch of ${MAX_PARALLEL_JOBS} jobs to finish."

                if ! wait; then
                    ((failed_jobs++))
                fi
            fi
        done

        if ! wait; then
            ((failed_jobs++))
        fi

        log_message "Finished processing ${round_directory}."
    done

    printf "\n"

    if (( failed_jobs > 0 )); then
        log_message \
            "The workflow completed with one or more failed docking jobs. Review the individual log files."
        return 1
    fi

    log_message "All ADCP cyclic peptide redocking jobs completed successfully."
}


main "$@"