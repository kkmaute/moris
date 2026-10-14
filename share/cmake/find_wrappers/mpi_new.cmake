# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Provide one stable MORIS target for either supported MPI implementation.
# Imported MPI targets carry compiler-specific include paths and link options.
# The logger includes mpi.h even in serial mode. Previously the unconditional
# MPI compiler wrappers supplied these requirements implicitly.
if(NOT TARGET ${MORIS}::mpi)
    add_library(${MORIS}::mpi INTERFACE IMPORTED GLOBAL)
    if(NOT MORIS_USE_MPI IN_LIST MORIS_MPI_LIBS AND MORIS_MPI_LIBS)
        message(FATAL_ERROR "MORIS_USE_MPI supported packages: ${MORIS_MPI_LIBS}")
    endif()
    set(moris_new_mpi_targets)
    foreach(language IN ITEMS C CXX)
        if(NOT TARGET MPI::MPI_${language})
            list(APPEND moris_new_mpi_targets MPI::MPI_${language})
        endif()
    endforeach()
    find_package(MPI REQUIRED COMPONENTS C CXX)
    # FindMPI targets can be created below the top-level directory in static mode.
    # Promote them so all component interfaces can refer to the same targets.
    foreach(mpi_target IN LISTS moris_new_mpi_targets)
        set_property(TARGET ${mpi_target} PROPERTY IMPORTED_GLOBAL TRUE)
    endforeach()
    target_link_libraries(${MORIS}::mpi INTERFACE MPI::MPI_CXX)
    set(MORIS_MPI_INCLUDE_DIRS "${MPI_CXX_INCLUDE_DIRS}")
    set(MORIS_MPI_LIBRARIES "${MPI_CXX_LIBRARIES}")
    if(MORIS_HAVE_PARALLEL)
        set(MORIS_MPI_DEFINITIONS MORIS_HAVE_PARALLEL)
    endif()
    set(MPI_FOUND_ONCE TRUE)
endif()
