#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Moris Packages ----------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory variables
foreach(MORIS_SOURCE_DIR ${MORIS_SOURCE_DIRS})
    string(REGEX REPLACE "([^\\/]*)\\/" "" DIR_NAME ${MORIS_SOURCE_DIR})
    
    # Build directories
    set(${DIR_NAME}_CONFIG_BUILD_DIR "${CONFIG_BUILD_DIR}/${MORIS_SOURCE_DIR}")
    
    # Install directories
    set(${DIR_NAME}_CONFIG_INSTALL_DIR "${CONFIG_INSTALL_DIR}/${MORIS_SOURCE_DIR}")
    set(${DIR_NAME}_HEADER_INSTALL_DIR "${HEADER_INSTALL_DIR}/${MORIS_SOURCE_DIR}")
    set(${DIR_NAME}_LIB_INSTALL_DIR "${LIB_INSTALL_DIR}/${MORIS_SOURCE_DIR}")
    set(${DIR_NAME}_EXE_INSTALL_DIR "${EXE_INSTALL_DIR}/${MORIS_SOURCE_DIR}")
endforeach()

# Add source directories
# This is done outside of the previous loop because cyclical dependencies may
# use the variables from packages that haven't been called yet
foreach(MORIS_SOURCE_DIR ${MORIS_SOURCE_DIRS})
    add_subdirectory(${MORIS_SOURCE_DIR})
endforeach()

