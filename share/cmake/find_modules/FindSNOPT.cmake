#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# SNOPT Find Module -------------------------------------------------------
# -------------------------------------------------------------------------

set(SNOPT_ENV_VARS
    $ENV{SNOPTDIR}
    $ENV{SNOPT_DIR}
    $ENV{snopt_DIR}
    $ENV{SNOPT_ROOT}
    $ENV{snopt_ROOT}
    $ENV{SNOPT_PATH}
    $ENV{snopt_PATH} )

find_path(SNOPT_LIBRARY_DIRS
    NAMES
    libsnopt.a
    HINTS
    ${SNOPT_ENV_VARS}
    PATH_SUFFIXES
    lib64
    lib )

find_library(SNOPT_LIBRARIES
    NAMES
    snopt
    HINTS
    ${SNOPT_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SNOPT DEFAULT_MSG 
    SNOPT_LIBRARY_DIRS
    SNOPT_LIBRARIES )

mark_as_advanced(SNOPT_LIBRARY_DIRS SNOPT_LIBRARIES)

