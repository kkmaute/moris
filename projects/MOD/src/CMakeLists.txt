# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Model Source

# Package version.
set(${MOD}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS model.hpp cl_Bspline.hpp cl_EqnManager.hpp cl_Model_Enums.hpp)

# Header directories.
set(MOD_INCLUDES #${ALG}
    #${LINALG}
)

# Include dependency directories.
foreach(MOD_INCLUDE IN LISTS MOD_INCLUDES)
    include_directories(${MORIS_PACKAGE_DIR}/${MOD_INCLUDE}/src)
endforeach()

# Model Library

# Sources.
set(LIB_SOURCES cl_Bspline.cpp cl_EqnManager.cpp)

# MORIS library dependencies.
set(LIB_DEPENDENCIES "")

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE LIB_DEPENDENCIES)
list(REMOVE_DUPLICATES LIB_DEPENDENCIES)
list(REVERSE LIB_DEPENDENCIES)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${MOD}-lib MOD_BUILD_TARGET ${LIB_SOURCES})
target_link_libraries(${MOD_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MOD_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${MOD})

# Third-party dependencies.
foreach(TPL IN LISTS MOD_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${MOD_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_link_libraries(
        ${MOD_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS}
    )
    target_compile_definitions(${MOD_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    MOD
    TARGET
    ${MOD_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${MOD_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${HEADER_INSTALL_DIR}/${MOD})

# Install library
install(TARGETS ${MOD}-lib DESTINATION ${LIB_INSTALL_DIR}/${MOD} EXPORT ${MOD}Targets)

# Install exported target definitions.
install(EXPORT ${MOD}Targets DESTINATION ${CONFIG_INSTALL_DIR}/${MOD})
