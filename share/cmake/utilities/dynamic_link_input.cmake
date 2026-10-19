#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

## This CMake function is used to dynamically link the input file
## Name : dynamic_link_input
## Params: path to input file
## Output: .so

function(dynamic_link_input
    	 target_name
		 base_name 
		 cpp_name 
		 so_includes)

set(SO_LIB_REQS
#    ${HMR}-lib
#    ${COM}-lib
#    ${WRK}-lib
#    ${MORIS_PETSC_LIBS}
#    ${MORIS_BOOST_LIBS}
#    ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
#    ${MORIS_MPI_LIBS}
#    ${MORIS_ARMADILLO_EIGEN_LIBS}
#    ${MORIS_SUPERLU_LIBS}
#    ${MORIS_LDLIBS}
#    ${MORIS_TRILINOS_LIBS}
#    ${MORIS_BASE_LIBS}
    )


add_library(${target_name} SHARED ${cpp_name})

target_include_directories(${target_name} PRIVATE ${SO_INCLUDES})    

target_link_libraries(${target_name} ${SO_LIB_REQS})

target_compile_definitions(${target_name} INTERFACE ${MORIS_DEFINITIONS})                          

set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${base_name}
                                                PREFIX      ""          )

endfunction()

# Capture the input target's build settings while MORIS is configured. Runtime
# input compilation can then use a separate build tree without changing MORIS.
function(generate_input_build_project target_name)
    # Record a default for both single-configuration and multi-configuration
    # generators; "default" also names builds with no CMAKE_BUILD_TYPE.
    set(INPUT_BUILD_CONFIG "${CMAKE_BUILD_TYPE}")
    if(CMAKE_CONFIGURATION_TYPES)
        list(GET CMAKE_CONFIGURATION_TYPES 0 INPUT_BUILD_CONFIG)
    elseif(NOT INPUT_BUILD_CONFIG)
        set(INPUT_BUILD_CONFIG default)
    endif()
    file(WRITE "${CMAKE_BINARY_DIR}/input_file/default-config" "${INPUT_BUILD_CONFIG}\n")
    # Scratch builds use one configuration at a time. Preserve the existing
    # build tool, including Ninja when the parent uses Ninja Multi-Config.
    if(CMAKE_GENERATOR MATCHES "Ninja")
        set(INPUT_BUILD_GENERATOR Ninja)
    else()
        set(INPUT_BUILD_GENERATOR "Unix Makefiles")
    endif()
    file(WRITE "${CMAKE_BINARY_DIR}/input_file/generator" "${INPUT_BUILD_GENERATOR}\n")

    set(INPUT_TARGET "${target_name}")
    set(INPUT_TARGET_PROPERTIES "")
    # Copy optional language/ABI settings only when defined on the parent
    # target, so absent properties retain CMake's normal defaults.
    foreach(property CXX_STANDARD CXX_STANDARD_REQUIRED CXX_EXTENSIONS
            POSITION_INDEPENDENT_CODE INTERPROCEDURAL_OPTIMIZATION)
        get_target_property(value ${target_name} ${property})
        if(NOT "${value}" STREQUAL "value-NOTFOUND")
            string(APPEND INPUT_TARGET_PROPERTIES
                "set_property(TARGET ${target_name} PROPERTY ${property} [==[${value}]==])\n")
        endif()
    endforeach()
    # Expand @variables@ first, then resolve target properties and $<CONFIG>
    # during generation. Each configuration gets a self-contained project.
    configure_file("${MORIS_CMAKE_DIR}/utilities/input_file.CMakeLists.txt.in"
        "${CMAKE_CURRENT_BINARY_DIR}/input_file.CMakeLists.txt.in" @ONLY)
    file(GENERATE
        OUTPUT "${CMAKE_BINARY_DIR}/input_file/$<IF:$<BOOL:$<CONFIG>>,$<CONFIG>,default>/CMakeLists.txt"
        INPUT "${CMAKE_CURRENT_BINARY_DIR}/input_file.CMakeLists.txt.in")
endfunction()





