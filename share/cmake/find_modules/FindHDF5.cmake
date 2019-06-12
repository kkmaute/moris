# HDF5 Find Module --------------------------------------------------------
# -------------------------------------------------------------------------
# NOTE: There is a FindHDF5.cmake file provided by the people at CMake, but
# the moris cluster does not have cmake-compatible HDF5. If this gets
# remedied, remove this file and use the CMake one.

set(HDF5_ENV_VARS
    $ENV{HDF5DIR}
    $ENV{HDF5_DIR}
    $ENV{hdf5_DIR}
    $ENV{HDF5_ROOT}
    $ENV{hdf5_ROOT}
    $ENV{HDF5_PATH}
    $ENV{hdf5_PATH} )

find_path(HDF5_INCLUDE_DIRS
	NAMES
    hdf5.h
    HINTS
    ${HDF5_ENV_VARS}
    PATH_SUFFIXES
    include )

find_path(HDF5_LIBRARY_DIRS
    NAMES
    libhdf5.a
    HINTS
    ${HDF5_ENV_VARS}
    PATH_SUFFIXES
    lib
    lib64 )

find_library(HDF5_hdf5
	NAMES
	hdf5
	HINTS
	${HDF5_LIBRARY_DIRS} )

find_library(HDF5_fortran
	NAMES
	hdf5_fortran
	HINTS
	${HDF5_LIBRARY_DIRS} )

find_library(HDF5_hl
	NAMES
	hdf5_hl
	HINTS
	${HDF5_LIBRARY_DIRS} )

find_library(HDF5_hl_fortran
	NAMES
	hdf5hl_fortran
	HINTS
	${HDF5_LIBRARY_DIRS} )

set(HDF5_LIBRARIES
    ${HDF5_hdf5}
    ${HDF5_fortran}
    ${HDF5_hl}
    ${HDF5_hl_fortran}
    CACHE FILEPATH "List of library paths." )


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5 DEFAULT_MSG 
	HDF5_INCLUDE_DIRS
	HDF5_LIBRARY_DIRS
	HDF5_LIBRARIES )

mark_as_advanced(HDF5_INCLUDE_DIRS
	HDF5_LIBRARY_DIRS
	HDF5_hdf5 HDF5_fortran
	HDF5_hl HDF5_hl_fortran
	HDF5_LIBRARIES )
