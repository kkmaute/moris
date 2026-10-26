# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
# Shared shell functions for moris.sh and the compile-only entry points.

moris_build_directory_for_mode()
{
    case "$1" in
        dbg) printf '%s\n' "${MORISBUILDDBG:?Set MORISBUILDDBG to the debug build directory}" ;;
        opt) printf '%s\n' "${MORISBUILDOPT:?Set MORISBUILDOPT to the optimized build directory}" ;;
        *) printf '%s\n' 'Error: mode must be dbg or opt' >&2; return 1 ;;
    esac
}

# Read values as data, including spaces and embedded equals signs. Never source
# CMakeCache.txt as shell code. Keys are fixed CMake identifiers from this file.
moris_input_cache_value()
{
    sed -n "s/^$1:[^=]*=//p" "$2"
}

# Use the same configuration selection for compiling an input and launching MORIS.
moris_input_configuration()
(
    set -eu
    moris_input_cache=$1/CMakeCache.txt
    if [ ! -f "$moris_input_cache" ]; then
        printf 'Error: no configured MORIS build at %s\n' "$1" >&2
        exit 1
    fi
    moris_input_build_type=$(moris_input_cache_value CMAKE_BUILD_TYPE "$moris_input_cache")
    moris_input_configuration=${MORIS_INPUT_CONFIGURATION-$moris_input_build_type}
    moris_input_configurations=$(moris_input_cache_value CMAKE_CONFIGURATION_TYPES "$moris_input_cache")
    if [ -n "$moris_input_configurations" ]; then
        if [ -z "$moris_input_configuration" ]; then
            printf 'Error: set MORIS_INPUT_CONFIGURATION to one of: %s\n' "$moris_input_configurations" >&2
            exit 1
        fi
        case ";$moris_input_configurations;" in
            *";$moris_input_configuration;"*) ;;
            *) printf 'Error: configuration is not available in the selected MORIS build: %s\n' "$moris_input_configuration" >&2; exit 1 ;;
        esac
    elif [ "$moris_input_configuration" != "$moris_input_build_type" ]; then
        printf '%s\n' 'Error: MORIS_INPUT_CONFIGURATION must match the selected MORIS build type.' >&2
        exit 1
    fi
    printf '%s\n' "$moris_input_configuration"
)

moris_build_input()
(
    # A subshell keeps build variables and positional arguments out of moris.sh.
    set -eu
    if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
        printf '%s\n' 'Usage: create_shared_object.sh <workdir> <MORIS build directory> <input basename or .cpp> [dbg]' >&2
        exit 1
    fi
    if [ -z "${MORISROOT-}" ]; then
        printf '%s\n' 'Error: MORISROOT must identify the MORIS source directory.' >&2
        exit 1
    fi
    moris_input_root=$(CDPATH= cd -- "$MORISROOT" && pwd -P)
    moris_input_work=$(CDPATH= cd -- "$1" && pwd -P)
    case "$2" in
        /*) moris_input_build=$2 ;;
        *) moris_input_build=$moris_input_root/$2 ;;
    esac
    moris_input_build=$(CDPATH= cd -- "$moris_input_build" && pwd -P)
    moris_input_cache=$moris_input_build/CMakeCache.txt
    if [ ! -f "$moris_input_cache" ]; then
        printf 'Error: no configured MORIS build at %s\n' "$moris_input_build" >&2
        exit 1
    fi

    case "$3" in
        *.cpp) moris_input_source=$3 ;;
        *) moris_input_source=$3.cpp ;;
    esac
    case "$moris_input_source" in
        /*) ;;
        *) moris_input_source=$moris_input_work/$moris_input_source ;;
    esac
    if [ ! -f "$moris_input_source" ]; then
        printf 'Error: input source does not exist: %s\n' "$moris_input_source" >&2
        exit 1
    fi
    moris_input_directory=$(CDPATH= cd -- "$(dirname -- "$moris_input_source")" && pwd -P)
    moris_input_name=$(basename -- "$moris_input_source" .cpp)
    moris_input_source=$moris_input_directory/$moris_input_name.cpp
    if [ "$#" -eq 4 ]; then
        if [ "$4" != dbg ]; then
            printf '%s\n' 'Error: the only supported legacy fourth argument is dbg.' >&2
            exit 1
        fi
        printf '%s\n' 'Note: dbg no longer creates a source link; debug information uses the original source path.'
    fi

    moris_input_configuration=$(moris_input_configuration "$moris_input_build")
    moris_input_config=$moris_input_build/input_build/$moris_input_configuration/MorisInputConfig.cmake
    if [ ! -f "$moris_input_config" ]; then
        moris_input_build_source=$(moris_input_cache_value CMAKE_HOME_DIRECTORY "$moris_input_cache")
        printf 'Error: external-input settings are missing: %s\n' "$moris_input_config" >&2
        printf 'Selected MORIS build: %s\n' "$moris_input_build" >&2
        printf 'Reconfigure the selected build using a checkout with the updated external-input CMake support (BUILD_MAIN must be ON):\n' >&2
        printf '  cmake -S "%s" -B "%s"\n' "${moris_input_build_source:-$moris_input_root}" "$moris_input_build" >&2
        printf 'Check MORISROOT, MORISBUILDDBG/MORISBUILDOPT, and the moris.sh selected by PATH if these paths are unexpected.\n' >&2
        exit 1
    fi

    # Keep the existing cache identity so switching from the Python helper does
    # not discard incremental builds. Distinct MORIS builds/configurations have
    # independent caches, while the .so stays beside the original source.
    moris_input_identity=$(printf '%s\n%s' "$moris_input_build" "$moris_input_configuration" | sha256sum)
    moris_input_identity=$(printf '%.16s' "$moris_input_identity")
    moris_input_local_build=$moris_input_directory/.moris-input-build/$moris_input_name/$moris_input_identity
    moris_input_generator=$(moris_input_cache_value CMAKE_GENERATOR "$moris_input_cache")
    moris_input_platform=$(moris_input_cache_value CMAKE_GENERATOR_PLATFORM "$moris_input_cache")
    moris_input_toolset=$(moris_input_cache_value CMAKE_GENERATOR_TOOLSET "$moris_input_cache")

    set -- -S "$moris_input_root/share/cmake/input" -B "$moris_input_local_build" \
        -G "$moris_input_generator" "-DMORIS_INPUT_CONFIG=$moris_input_config" \
        "-DMORIS_INPUT_SOURCE=$moris_input_source" \
        "-DMORIS_INPUT_OUTPUT_DIRECTORY=$moris_input_directory" \
        "-DCMAKE_BUILD_TYPE=$moris_input_configuration"
    if [ -n "$moris_input_platform" ]; then
        set -- "$@" -A "$moris_input_platform"
    fi
    if [ -n "$moris_input_toolset" ]; then
        set -- "$@" -T "$moris_input_toolset"
    fi
    cmake "$@"

    set -- --build "$moris_input_local_build" --target moris_external_input
    if [ -n "$moris_input_configuration" ]; then
        set -- "$@" --config "$moris_input_configuration"
    fi
    cmake "$@"
    printf 'Input library: %s/%s.so\n' "$moris_input_directory" "$moris_input_name"
)
