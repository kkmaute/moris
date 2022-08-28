#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# moris Doxygen -----------------------------------------------------------
# -------------------------------------------------------------------------

# Find Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
    # Input proper paths #> need to replace CMAKE_SOURCE_DIR in there
    configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    
    # Create a doxygen make target
    add_custom_target(doxy
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating moris documentation with Doxygen."
        VERBATIM )
endif()

