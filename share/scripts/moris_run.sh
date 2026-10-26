#!/bin/sh
# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
# Usage: moris_run.sh <dbg|opt> <number of MPI ranks> <input basename or .so> [run prefix]
set -eu
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <dbg|opt> <number of MPI ranks> <input basename or .so> [run prefix]" >&2
    exit 1
fi
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/moris_input_build.sh"
build_dir=$(moris_build_directory_for_mode "$1")
case "$2" in
    ''|*[!0-9]*|0) echo "Error: MPI rank count must be a positive integer" >&2; exit 1 ;;
esac
case "$build_dir" in
    /*) ;;
    *) build_dir=${MORISROOT:?Set MORISROOT}/$build_dir ;;
esac
build_dir=$(CDPATH= cd -- "$build_dir" && pwd -P)
configuration=$(moris_input_configuration "$build_dir")
executable_file=$build_dir/input_build/$configuration/MorisExecutable.txt
if [ -f "$executable_file" ]; then
    IFS= read -r moris_executable < "$executable_file"
else
    # Run-only usage remains compatible with builds configured before the export.
    configurations=$(moris_input_cache_value CMAKE_CONFIGURATION_TYPES "$build_dir/CMakeCache.txt")
    if [ -n "$configurations" ]; then
        moris_executable=$build_dir/projects/mains/$configuration/moris
    else
        moris_executable=$build_dir/projects/mains/moris
    fi
fi
if [ ! -x "$moris_executable" ]; then
    printf 'Error: MORIS executable is missing or not executable: %s\n' "$moris_executable" >&2
    exit 1
fi

input_file=${3%.so}.so
if [ ! -f "$input_file" ]; then
    echo "Error: input library does not exist: $input_file" >&2
    exit 1
fi
input_dir=$(CDPATH= cd -- "$(dirname -- "$input_file")" && pwd)
input_file=$input_dir/$(basename -- "$input_file")

# Preserve the legacy optional multiword run prefix without eval or globbing.
# Shell quoting inside the prefix is not interpreted; use a wrapper script for it.
rank_count=$2
set -f
run_prefix=${4-}
set -- $run_prefix
exec mpirun -np "$rank_count" "$@" "$moris_executable" "$input_file"
