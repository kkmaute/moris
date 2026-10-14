# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Build MORIS modules as objects in shared mode, retaining their public target names.
# Keep the static build path unchanged.
function(moris_add_module name result)
    if(MORIS_HAVE_SHARED)
        set(implementation "${name}-objects")
        add_library(${implementation} OBJECT ${ARGN})
        set_target_properties(${implementation} PROPERTIES POSITION_INDEPENDENT_CODE ON)
        add_library(${name} INTERFACE)
        set_property(GLOBAL APPEND PROPERTY MORIS_OBJECT_TARGETS ${implementation})
        set_property(GLOBAL APPEND PROPERTY MORIS_MODULE_TARGETS ${name})
    else()
        set(implementation "${name}")
        add_library(${implementation} STATIC ${ARGN})
    endif()
    if(TARGET moris_runtime)
        target_link_libraries(${implementation} PUBLIC moris_runtime)
    endif()
    # Module headers expose C++17 types to downstream target-based consumers.
    target_compile_features(${implementation} PUBLIC cxx_std_17)
    set(${result} "${implementation}" PARENT_SCOPE)
endfunction()

# Import TPL targets in the top-level directory so they are visible both to
# module compilation and to the final shared-library link.
macro(moris_initialize_shared_library)
    if(MORIS_HAVE_SHARED)
        # Preserve TPL search paths in installed binaries, and locate libmoris
        # relative to installed executables rather than the build tree.
        set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
        if(APPLE)
            set(MORIS_ORIGIN "@loader_path")
        else()
            set(MORIS_ORIGIN "$ORIGIN")
        endif()
        set(CMAKE_INSTALL_RPATH "${MORIS_ORIGIN}/../${LIB_INSTALL_DIR}")
        set(MORIS_SHARED_TPLS)
        foreach(directory IN LISTS MORIS_SOURCE_DIRS)
            get_filename_component(component "${directory}" NAME)
            list(APPEND MORIS_SHARED_TPLS ${${component}_TPL_DEPENDENCIES})
        endforeach()
        list(TRANSFORM MORIS_SHARED_TPLS TOLOWER)
        list(REMOVE_DUPLICATES MORIS_SHARED_TPLS)
        set(MORIS_SHARED_TPL_RPATH)
        foreach(tpl IN LISTS MORIS_SHARED_TPLS)
            include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
            # CMake's automatic install RPATH excludes build-tree dependencies.
            # Preserve paths to locally built shared TPLs as well.
            string(TOUPPER "${tpl}" tpl_upper)
            foreach(library IN LISTS ${tpl_upper}_LIBRARIES)
                if(IS_ABSOLUTE "${library}" AND library MATCHES "\\.(so|dylib)(\\.[0-9.]+)?$")
                    get_filename_component(library_directory "${library}" DIRECTORY)
                    list(APPEND MORIS_SHARED_TPL_RPATH "${library_directory}")
                endif()
            endforeach()
        endforeach()
        list(REMOVE_DUPLICATES MORIS_SHARED_TPL_RPATH)
        list(APPEND CMAKE_INSTALL_RPATH ${MORIS_SHARED_TPL_RPATH})
    endif()
endmacro()

function(moris_finalize_shared_library)
    if(NOT MORIS_HAVE_SHARED)
        return()
    endif()

    get_property(objects GLOBAL PROPERTY MORIS_OBJECT_TARGETS)
    get_property(modules GLOBAL PROPERTY MORIS_MODULE_TARGETS)
    # These header-only modules also participate in the original dependency cycles.
    foreach(component COR ASR EXC CHR LINALG)
        if(TARGET ${component}-lib)
            list(APPEND modules ${component}-lib)
        endif()
    endforeach()
    if(NOT objects)
        return()
    endif()

    add_library(moris_shared_usage INTERFACE)
    add_library(moris_shared SHARED)
    set_target_properties(
        moris_shared
        PROPERTIES OUTPUT_NAME moris
                   LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                   VERSION "${MORIS_VERSION}"
                   SOVERSION "${MORIS_VERSION_MAJOR}"
                   INSTALL_RPATH "${MORIS_ORIGIN};${MORIS_SHARED_TPL_RPATH}"
    )
    if(TARGET SDF-exe)
        set_property(
            TARGET SDF-exe
            PROPERTY INSTALL_RPATH
                     "${MORIS_ORIGIN}/../../../${LIB_INSTALL_DIR};${MORIS_SHARED_TPL_RPATH}"
        )
    endif()

    # Flatten usage requirements before replacing any links. No object target may
    # depend on the public component interfaces, which link back to moris_shared.
    set(external_links)
    foreach(target IN LISTS objects modules)
        foreach(
            property
            INCLUDE_DIRECTORIES
            SYSTEM_INCLUDE_DIRECTORIES
            COMPILE_DEFINITIONS
            COMPILE_OPTIONS
            COMPILE_FEATURES
            LINK_OPTIONS
            LINK_DIRECTORIES
        )
            get_target_property(value ${target} INTERFACE_${property})
            if(value)
                set_property(
                    TARGET moris_shared_usage APPEND PROPERTY INTERFACE_${property} "${value}"
                )
            endif()
        endforeach()
        foreach(property LINK_LIBRARIES INTERFACE_LINK_LIBRARIES)
            get_target_property(links ${target} ${property})
            if(links)
                foreach(link IN LISTS links)
                    if(NOT link IN_LIST modules AND NOT link IN_LIST objects)
                        list(APPEND external_links "${link}")
                    endif()
                endforeach()
            endif()
        endforeach()
    endforeach()
    list(REMOVE_DUPLICATES external_links)
    target_link_libraries(moris_shared_usage INTERFACE ${external_links})
    target_compile_definitions(moris_shared_usage INTERFACE ${MORIS_DEFINITIONS})
    target_compile_features(moris_shared_usage INTERFACE cxx_std_17)

    foreach(target IN LISTS objects)
        set_property(TARGET ${target} PROPERTY LINK_LIBRARIES moris_shared_usage)
        set_property(TARGET ${target} PROPERTY INTERFACE_LINK_LIBRARIES "")
        target_sources(moris_shared PRIVATE $<TARGET_OBJECTS:${target}>)
    endforeach()
    target_link_libraries(moris_shared PUBLIC moris_shared_usage)
    foreach(target IN LISTS modules)
        set_property(TARGET ${target} PROPERTY INTERFACE_LINK_LIBRARIES moris_shared)
    endforeach()

    install(
        TARGETS moris_shared moris_shared_usage
        EXPORT morisSharedTargets
        LIBRARY DESTINATION ${LIB_INSTALL_DIR}
        RUNTIME DESTINATION ${EXE_INSTALL_DIR}
        ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    )
    install(EXPORT morisSharedTargets DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS})

    # All consumers of the monolithic library need its complete TPL dependency set,
    # even when find_package requests only one MORIS component.
    configure_file(
        "${MORIS_CMAKE_DIR}/morisSharedConfig.cmake.in"
        "${CONFIG_BUILD_DIR}/${MORIS}/morisSharedConfig.cmake" @ONLY
    )
    install(FILES "${CONFIG_BUILD_DIR}/${MORIS}/morisSharedConfig.cmake"
            DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS}
    )
endfunction()
