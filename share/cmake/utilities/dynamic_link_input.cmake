# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Build an input-file plugin loaded by the MORIS executable with dlopen().
# Keep the unprefixed filename expected by existing example and workflow scripts.
function(dynamic_link_input target_name base_name cpp_name so_includes)
    add_library(${target_name} SHARED "${cpp_name}")
    # Callers pass both quoted lists and expanded lists of include directories.
    # Expanded lists place every directory after the first in ARGN.
    target_include_directories(${target_name} PRIVATE ${so_includes} ${ARGN})
    target_compile_definitions(${target_name} PRIVATE ${MORIS_DEFINITIONS})
    # Public MORIS headers include mpi.h; inherit MPI requirements even when
    # using a regular C++ compiler instead of an MPI compiler wrapper.
    if(TARGET ${MORIS}::mpi)
        target_link_libraries(${target_name} PRIVATE ${MORIS}::mpi)
    endif()
    if(MORIS_HAVE_SHARED)
        target_link_libraries(${target_name} PRIVATE moris_shared)
    endif()
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME "${base_name}" PREFIX "")
endfunction()
