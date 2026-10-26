#!/bin/sh
# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
# Usage: create_shared_object.sh <workdir> <MORIS build directory> <input basename or .cpp> [dbg]
# Compile the original source in a problem-local CMake build. The optional legacy
# dbg argument is accepted, but no placeholder, symbolic link, or lock is needed.
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$script_dir/moris_input_build.sh"
moris_build_input "$@"
