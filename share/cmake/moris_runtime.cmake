# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Common legacy runtime requirements are propagated through MORIS targets, not
# directory-wide link_libraries(), so unrelated targets do not inherit them.
add_library(moris_runtime INTERFACE)
string(CONCAT LIB_SEARCH_PATH "$ENV{LD_LIBRARY_PATH}" ":" "$ENV{SPACK_LINK_DIRS}")
string(REPLACE ":" ";" LIB_SEARCH_PATH "${LIB_SEARCH_PATH}")
list(REMOVE_DUPLICATES LIB_SEARCH_PATH)

# Explicit installation prefixes override stale cache entries in existing builds.
function(moris_find_runtime variable library prefix_variable)
    if(NOT "$ENV{${prefix_variable}}" STREQUAL "")
        unset(${variable} CACHE)
        unset(${variable})
        find_library(
            ${variable}
            NAMES ${library}
            PATHS "$ENV{${prefix_variable}}"
            PATH_SUFFIXES lib lib64
            NO_DEFAULT_PATH
        )
    else()
        find_library(${variable} NAMES ${library} HINTS ${LIB_SEARCH_PATH})
    endif()
    if(NOT ${variable})
        message(
            FATAL_ERROR
                "Cannot find ${library}. Set ${prefix_variable} to its installation prefix or ${variable} to its library file."
        )
    endif()
    message(STATUS "${variable} = ${${variable}}")
    target_link_libraries(moris_runtime INTERFACE "${${variable}}")
endfunction()

if(MORIS_USE_INTEL)
    moris_find_runtime(MORIS_LIB_Ifcore ifcore INTEL_RUNTIME_DIR)
    moris_find_runtime(MORIS_LIB_Ifcoremt ifcoremt INTEL_RUNTIME_DIR)
    moris_find_runtime(MORIS_LIB_Rt rt RT_DIR)
endif()

# Keep explicit Fortran overrides; otherwise ask the active compiler before
# searching the system, avoiding a preference for the oldest installed SONAME.
if(NOT MORIS_LIB_GFortran)
    if(DEFINED ENV{GFORTLIB} AND NOT "$ENV{GFORTLIB}" STREQUAL "")
        find_library(
            MORIS_LIB_GFortran NAMES "$ENV{GFORTLIB}" HINTS "$ENV{GFORTLIB_PATH}"
                                                            ${LIB_SEARCH_PATH}
        )
    else()
        execute_process(
            COMMAND "${CMAKE_CXX_COMPILER}" -print-file-name=libgfortran.so
            OUTPUT_VARIABLE moris_gfortran
            OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
            RESULT_VARIABLE moris_gfortran_status
        )
        if(moris_gfortran_status EQUAL 0 AND IS_ABSOLUTE "${moris_gfortran}"
           AND EXISTS "${moris_gfortran}"
        )
            set(MORIS_LIB_GFortran "${moris_gfortran}" CACHE FILEPATH "Fortran runtime library"
                                                             FORCE
            )
        else()
            find_library(
                MORIS_LIB_GFortran NAMES gfortran libgfortran.so.5 libgfortran.so.4 libgfortran.so.3
                HINTS ${LIB_SEARCH_PATH}
            )
        endif()
    endif()
endif()
if(NOT MORIS_LIB_GFortran)
    message(
        FATAL_ERROR
            "Cannot find the Fortran runtime. Set MORIS_LIB_GFortran or GFORTLIB/GFORTLIB_PATH."
    )
endif()
target_link_libraries(moris_runtime INTERFACE "${MORIS_LIB_GFortran}")

# Normalize the historical mixed-case spelling while preserving old cache overrides.
if(MORIS_Lib_Crypto AND NOT MORIS_LIB_Crypto)
    set(MORIS_LIB_Crypto "${MORIS_Lib_Crypto}" CACHE FILEPATH "Crypto library")
endif()
unset(MORIS_Lib_Crypto CACHE)
moris_find_runtime(MORIS_LIB_Dl dl DL_DIR)
moris_find_runtime(MORIS_LIB_Ssl ssl OPENSSL_DIR)
moris_find_runtime(MORIS_LIB_Crypto crypto OPENSSL_DIR)
moris_find_runtime(MORIS_LIB_Z z ZLIB_DIR)

install(TARGETS moris_runtime EXPORT morisRuntimeTargets)
install(EXPORT morisRuntimeTargets DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS})
