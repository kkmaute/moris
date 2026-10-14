# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Chronos Source

# Public headers.
set(HEADERS chronos.hpp cl_Stopwatch.hpp cl_Waitbar.hpp)

# Chronos Library

add_library(${CHR}-lib INTERFACE)
target_include_directories(
    ${CHR}-lib INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                         $<INSTALL_INTERFACE:${${CHR}_HEADER_INSTALL_DIR}>
)

foreach(TPL ${CHR_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${CHR}-lib INTERFACE ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${CHR}-lib INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${CHR}-lib INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Chronos Executable

# Generate and install the common component package configuration.
moris_install_component_package(CHR TARGET ${CHR}-lib DEPENDENCIES TPLS ${CHR_TPL_DEPENDENCIES})

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${CHR}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${CHR}-lib DESTINATION ${${CHR}_LIB_INSTALL_DIR} EXPORT ${CHR}Targets)

# Install exported target definitions.
install(EXPORT ${CHR}Targets DESTINATION ${${CHR}_CONFIG_INSTALL_DIR})
