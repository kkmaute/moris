# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

message(STATUS "Added structures - contact example sub-directory")

add_subdirectory(Parabolic_Indenter)

# Nonconformal contact mapping requires ArborX.
if(MORIS_HAVE_ARBORX)
    add_subdirectory(Coulomb_Shear)
endif()
