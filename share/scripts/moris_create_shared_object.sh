#!/bin/sh
# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
# Usage: moris_create_shared_object.sh <dbg|opt> <input basename or .cpp> [dbg]
set -eu
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <dbg|opt> <input basename or .cpp> [dbg]" >&2
    exit 1
fi
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/moris_input_build.sh"
build_dir=$(moris_build_directory_for_mode "$1")
if [ "$#" -eq 3 ]; then
    moris_build_input . "$build_dir" "$2" "$3"
else
    moris_build_input . "$build_dir" "$2"
fi
