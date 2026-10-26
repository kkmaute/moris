# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.

# Export the existing input target's requirements for independent, problem-local
# builds. Exported generator expressions are evaluated on the external target,
# preserving C++-only and configuration-specific flags.
function(moris_export_input_build input_target)
    add_library(moris_input_build_flags INTERFACE)
    get_target_property(input_source_dir ${input_target} SOURCE_DIR)
    foreach(property INCLUDE_DIRECTORIES COMPILE_DEFINITIONS COMPILE_OPTIONS COMPILE_FEATURES
                     LINK_OPTIONS LINK_DIRECTORIES
    )
        get_target_property(value ${input_target} ${property})
        if(property STREQUAL "COMPILE_DEFINITIONS")
            # Directory definitions are combined at generation time and are not
            # necessarily present in the target's configure-time property.
            get_directory_property(
                directory_definitions DIRECTORY "${input_source_dir}" COMPILE_DEFINITIONS
            )
            if(NOT value)
                set(value)
            endif()
            list(APPEND value ${directory_definitions})
        endif()
        if(value AND property STREQUAL "INCLUDE_DIRECTORIES")
            # Legacy mains lists contain obsolete include paths. Compilers ignore
            # them, but imported CMake targets reject nonexistent directories.
            set(existing_includes)
            foreach(directory IN LISTS value)
                if(directory MATCHES "^\\$<" OR IS_DIRECTORY "${directory}")
                    list(APPEND existing_includes "${directory}")
                endif()
            endforeach()
            set(value "${existing_includes}")
        endif()
        if(value)
            set_property(TARGET moris_input_build_flags PROPERTY INTERFACE_${property} "${value}")
        endif()
        if(TARGET moris_shared_usage)
            get_target_property(shared_value moris_shared_usage INTERFACE_${property})
            if(shared_value)
                set_property(
                    TARGET moris_input_build_flags APPEND PROPERTY INTERFACE_${property}
                                                                   "${shared_value}"
                )
            endif()
        endif()
    endforeach()

    # MPI is imported by the external project with the selected build's wrapper.
    # Other directory-level links (OpenMP and optional absolute runtime libraries)
    # must remain available to input libraries as well.
    get_target_property(links ${input_target} LINK_LIBRARIES)
    list(REMOVE_ITEM links moris_shared ${MORIS}::mpi)
    if(links)
        target_link_libraries(moris_input_build_flags INTERFACE ${links})
    endif()
    if(TARGET ${MORIS}::mpi)
        target_link_libraries(moris_input_build_flags INTERFACE MPI::MPI_CXX)
        set(MORIS_INPUT_NEEDS_MPI TRUE)
    endif()
    if(TARGET OpenMP::OpenMP_CXX)
        set(MORIS_INPUT_NEEDS_OPENMP TRUE)
    endif()
    if(TARGET moris_shared)
        set(MORIS_INPUT_SHARED_LIBRARY "$<TARGET_FILE:moris_shared>")
    endif()

    # Undefined per-configuration IPO settings inherit the global default.
    # Preserve that distinction when serializing values into the input project.
    foreach(configuration DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
        if(NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION_${configuration})
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_${configuration}
                "${CMAKE_INTERPROCEDURAL_OPTIMIZATION}"
            )
        endif()
    endforeach()

    set(input_config_dir "${CMAKE_BINARY_DIR}/input_build")
    # Export the actual executable path as plain data for the shell runner. This
    # handles configuration subdirectories and custom target output locations.
    if(TARGET moris)
        file(GENERATE OUTPUT "${input_config_dir}/$<CONFIG>/MorisExecutable.txt"
             CONTENT "$<TARGET_FILE:moris>\n"
        )
    endif()
    export(TARGETS moris_input_build_flags FILE "${input_config_dir}/MorisInputTargets.cmake")
    configure_file(
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../input/MorisInputConfig.cmake.in"
        "${input_config_dir}/MorisInputConfig.cmake.in" @ONLY
    )
    file(GENERATE OUTPUT "${input_config_dir}/$<CONFIG>/MorisInputConfig.cmake"
         INPUT "${input_config_dir}/MorisInputConfig.cmake.in"
    )
endfunction()
