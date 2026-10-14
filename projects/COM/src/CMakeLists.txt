# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Communications Source

# Package version.
set(${COM}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS cl_Communication_Enums.hpp cl_Communication_Manager.hpp cl_Communication_Tools.hpp)

# Communications Library

# Library sources.
set(LIB_SOURCES cl_Communication_Manager.cpp cl_Communication_Tools.cpp)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${LINALG}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${COM}-lib COM_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${COM_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${COM}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${COM_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${COM_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${COM})

# Third-party dependencies.
foreach(TPL IN LISTS COM_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${COM_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${COM_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${COM_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    COM
    TARGET
    ${COM_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${COM_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${HEADER_INSTALL_DIR}/${COM})

# Install and export targets.
install(
    TARGETS ${COM}-lib
    EXPORT ${COM}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${COM}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${COM}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${COM}
)

# Install exported target definitions.
install(EXPORT ${COM}Targets DESTINATION ${CONFIG_INSTALL_DIR}/${COM})
