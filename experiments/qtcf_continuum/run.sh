#!/usr/bin/env bash
set -euo pipefail

usage() {
    printf 'Usage: %s [new-output-directory]\n' "$0"
}
if [[ ${1:-} == --help || ${1:-} == -h ]]; then
    usage
    exit 0
fi
if (( $# > 1 )); then
    usage >&2
    exit 2
fi

source_dir=$(cd "$(dirname "$0")" && pwd -P)
repo=$(git -C "$source_dir" rev-parse --show-toplevel)
revision=e19e600b4d0bca9ed12d6e6fa6369923ee6329b6
command -v julia >/dev/null || { echo 'Install Julia 1.11.5 first.' >&2; exit 1; }
if ! python3 -c 'import pypdf, sympy'; then
    echo 'Install the Python dependencies pypdf and sympy first.' >&2
    exit 1
fi

# Give each run its own output directory so earlier results stay intact.
if (( $# == 1 )); then
    if [[ -z "$1" || -e "$1" || -L "$1" ]]; then
        echo 'Choose a new output directory.' >&2
        exit 2
    fi
    mkdir -p -- "$(dirname "$1")"
    mkdir -- "$1"
    output_dir=$(cd "$1" && pwd -P)
else
    output_dir=$(mktemp -d "$PWD/qtcf-continuum-results.XXXXXX")
fi
scratch=$(mktemp -d "${TMPDIR:-/tmp}/qtcf-continuum.XXXXXX")
finish() {
    local status=$?
    if (( status == 0 )); then
        rm -rf -- "$scratch"
    else
        printf 'Run stopped. Work files: %s\nLogs: %s\n' "$scratch" "$output_dir" >&2
    fi
}
trap finish EXIT

# Work on the pinned package snapshot, using its tracked dependency manifest.
cp "$source_dir/"*.jl "$source_dir/"*.py "$scratch/"
mkdir -p "$scratch/frozen/RFA"
git -C "$repo" archive "$revision" Project.toml LICENSE src ext test |
    tar -xf - -C "$scratch/frozen/RFA"
patch -F 0 -s -p1 -d "$scratch/frozen/RFA" < "$source_dir/rfa_experiment.patch"

cd "$scratch"
export JULIA_NUM_THREADS=1 JULIA_PKG_PRECOMPILE_AUTO=0
run_julia() {
    julia --startup-file=no --compiled-modules=existing --project=frozen/RFA/test "$@"
}
run_logged() {
    local log=$1
    shift
    printf 'Running %s\n' "$log"
    "$@" > "$output_dir/$log.log" 2>&1
}

if [[ ${QTCF_INSTANTIATE:-1} != 0 ]]; then
    run_logged dependencies run_julia -e 'using Pkg; Pkg.instantiate()'
fi
run_julia run.jl 2>&1 | tee "$output_dir/full_run.log"

# Audit the recorded fits before producing tables and pole plots.
run_logged history_degrees run_julia certify_curve_degrees.jl --histories
run_logged returned_degrees run_julia certify_curve_degrees.jl
run_logged interval_certification run_julia certify_intervals.jl
run_logged startup run_julia check_startup.jl
run_logged matched_midpoint run_julia check_matched_midpoint.jl
run_logged pole_quality run_julia check_pole_quality.jl
run_logged curve_spectra run_julia check_curve_spectra.jl
run_logged fixed_nodes run_julia run_audit.jl
run_logged stress_3_3 run_julia stress_3_3.jl
run_logged claim_1 python3 claim_1.py supports.csv
run_logged interval_replay python3 replay_interval_reductions.py
run_logged figures run_julia render_grouped_figures.jl
run_logged tables python3 build_tables.py

cp ./*.csv ./group_*_rows.tex provenance.txt qtcf_normalized_update_figures.pdf "$output_dir/"
cp -R interval_certification "$output_dir/"
printf 'Results: %s\n' "$output_dir"
