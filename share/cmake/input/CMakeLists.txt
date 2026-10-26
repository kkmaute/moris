# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
cmake_minimum_required(VERSION 3.17...4.4)

if(NOT EXISTS "${MORIS_INPUT_CONFIG}")
    message(FATAL_ERROR "Provide MORIS_INPUT_CONFIG from a configured MORIS build.")
endif()
include("${MORIS_INPUT_CONFIG}")
project(MorisExternalInput LANGUAGES CXX)

if(NOT IS_ABSOLUTE "${MORIS_INPUT_SOURCE}" OR NOT EXISTS "${MORIS_INPUT_SOURCE}")
    message(FATAL_ERROR "MORIS_INPUT_SOURCE must name an existing absolute .cpp path.")
endif()
if(NOT MORIS_INPUT_OUTPUT_DIRECTORY)
    message(FATAL_ERROR "Provide MORIS_INPUT_OUTPUT_DIRECTORY for the resulting .so.")
endif()
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE STREQUAL MORIS_INPUT_CONFIGURATION)
    message(FATAL_ERROR "The input configuration must match the selected MORIS configuration.")
endif()
if(MORIS_INPUT_NEEDS_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()
if(MORIS_INPUT_NEEDS_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()
find_package(Threads REQUIRED)
include("${MORIS_INPUT_TARGETS}")

get_filename_component(input_name "${MORIS_INPUT_SOURCE}" NAME)
string(REGEX REPLACE "\\.cpp$" "" input_name "${input_name}")
add_library(moris_external_input SHARED "${MORIS_INPUT_SOURCE}")
target_link_libraries(moris_external_input PRIVATE moris_input_build_flags)
if(MORIS_INPUT_SHARED_LIBRARY)
    if(NOT EXISTS "${MORIS_INPUT_SHARED_LIBRARY}")
        message(
            FATAL_ERROR "Build the selected MORIS shared library before building external inputs."
        )
    endif()
    target_link_libraries(moris_external_input PRIVATE "${MORIS_INPUT_SHARED_LIBRARY}")
endif()
set_target_properties(
    moris_external_input PROPERTIES PREFIX "" OUTPUT_NAME "${input_name}"
                                    LIBRARY_OUTPUT_DIRECTORY "${MORIS_INPUT_OUTPUT_DIRECTORY}"
)

foreach(configuration IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER "${configuration}" configuration_upper)
    set_target_properties(
        moris_external_input PROPERTIES LIBRARY_OUTPUT_DIRECTORY_${configuration_upper}
                                        "${MORIS_INPUT_OUTPUT_DIRECTORY}"
    )
endforeach()
# Do not implicitly turn all exported MORIS headers into system headers.
set_target_properties(moris_external_input PROPERTIES NO_SYSTEM_FROM_IMPORTED ON)
