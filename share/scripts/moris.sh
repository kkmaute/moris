#!/bin/sh
# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
# Usage: moris.sh <dbg|opt> <number of MPI ranks> <input basename or .cpp> [run prefix]
set -eu
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <dbg|opt> <number of MPI ranks> <input basename or .cpp> [run prefix]" >&2
    exit 1
fi
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/moris_input_build.sh"
build_dir=$(moris_build_directory_for_mode "$1")

# Configure and build the original input through CMake, then run only on success.
# Compile-only wrappers source the same shell functions; no Python or cmake -P
# orchestration is involved.
moris_build_input . "$build_dir" "$3"
input_name=${3%.cpp}
exec "$script_dir/moris_run.sh" "$1" "$2" "$input_name" "${4-}"
