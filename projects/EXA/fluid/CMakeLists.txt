#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

message(STATUS "Added fluid example sub-directory")

add_subdirectory(compressible)
add_subdirectory(laminar)
add_subdirectory(turbulent)
add_subdirectory(BodyFitted)

