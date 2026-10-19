#!/bin/sh

#====================================================================================
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#====================================================================================

# Usage: create_shared_object.sh workdir builddir input_basename [dbg]
# Relative build directories are resolved under MORISROOT. The source and output
# are resolved under workdir. No files in the MORIS source/build tree are changed.
# Stop on failed commands so a failed compilation cannot publish a new .so.
set -eu

fail()
{
    echo "Error: $*" >&2
    exit 1
}

# Validate and resolve paths before allocating scratch space. Quoting supports
# spaces in the working directory, build directory, and input filename.
[ -n "${MORISROOT:-}" ] || fail "MORISROOT environment variable not set"
[ -d "$MORISROOT" ] || fail "$MORISROOT does not exist"
[ "$#" -ge 3 ] && [ "$#" -le 4 ] || fail "Usage: $0 workdir builddir input_basename [dbg]"
[ -d "$1" ] || fail "Working directory $1 does not exist"
workdir=$(cd "$1" && pwd)
case "$2" in
    /*) builddir=$2 ;;
    *) builddir="$MORISROOT/$2" ;;
esac
[ -d "$builddir" ] || fail "Build directory $builddir does not exist"
builddir=$(cd "$builddir" && pwd)
case "$3" in
    /*) source_file="$3.cpp"; output_file="$3.so" ;;
    *) source_file="$workdir/$3.cpp"; output_file="$workdir/$3.so" ;;
esac
[ -f "$source_file" ] || fail "Source file $source_file does not exist"
case "${4:-}" in
    ""|dbg) ;;
    *) fail "Unknown option: $4" ;;
esac

# Read the project exported by the MORIS CMake configure/generate step.
# Runtime compilation never invokes CMake on the original MORIS build tree.
project_dir="$builddir/input_file"
[ -f "$project_dir/default-config" ] || fail "Reconfigure MORIS with CMake to generate the input-file build project"
generator=$(cat "$project_dir/generator")
config=${MORIS_INPUT_BUILD_CONFIG:-$(cat "$project_dir/default-config")}
project_dir="$project_dir/$config"
[ -f "$project_dir/CMakeLists.txt" ] || fail "No input-file build project for configuration $config"

# TMPDIR belongs to the caller; never redirect it into MORISROOT.
# mktemp creates a different directory for every invocation. No shared source,
# object file, or input_file.locked is used, so builds need no common lock.
scratch=$(mktemp -d "${TMPDIR:-/tmp}/moris-input.XXXXXXXX")
cleanup()
{
    if [ "${4:-}" = dbg ]; then
        echo "Input build files retained in $scratch"
    else
        rm -rf "$scratch"
    fi
}
# Pass the original arguments to cleanup so the legacy dbg option keeps scratch
# artifacts available for inspection. Debug symbols still name the real source.
trap 'cleanup "$@"' 0
# Convert termination signals into exits so the cleanup trap also runs on them.
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

# Copy only the generated project, so even CMake's source directory is writable.
cp "$project_dir/CMakeLists.txt" "$scratch/CMakeLists.txt"
cmake -S "$scratch" -B "$scratch/build" -G "$generator" \
    "-DMORIS_INPUT_FILE:FILEPATH=$source_file"
cmake --build "$scratch/build" --target shared_object_file
[ -f "$scratch/build/lib/input_file.so" ] || fail "Compiler did not produce input_file.so"
# Publish only a successful build. Concurrent callers should use distinct output
# paths (different basenames or working directories) to retain both libraries.
mv -f "$scratch/build/lib/input_file.so" "$output_file"
echo "Created $output_file"
