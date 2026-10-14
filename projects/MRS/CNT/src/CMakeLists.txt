# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Containers Source

# Public headers.
set(HEADERS
    cl_Array.hpp
    cl_Bi_Map.hpp
    cl_Bitbool.hpp
    cl_Bitset.hpp
    cl_BoostBitset.hpp
    cl_Vector.hpp
    cl_Design_Variable.hpp
    cl_Dist_Map.hpp
    cl_Map.hpp
    cl_Mini_Maps.hpp
    cl_Variant.hpp
    cl_Validator.hpp
    cl_Parameter.hpp
    cl_Parameter_List.hpp
    cl_Tuple.hpp
    containers.hpp
    fn_zip.hpp
    cl_Fortran.hpp
)

# Library sources.
set(LIB_SOURCES cl_Vector.cpp cl_Design_Variable.cpp cl_Variant.cpp cl_Validator.cpp
                cl_Parameter.cpp cl_Parameter_List.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${ENM}-lib ${LINALG}-lib)

# Containers Library

moris_add_module(${CNT}-lib CNT_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${CNT_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${CNT}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${CNT_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(
    ${CNT_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${CNT} ARCHIVE_OUTPUT_DIRECTORY
                                                      "${CMAKE_CURRENT_BINARY_DIR}/${LIB}"
)

foreach(TPL IN LISTS CNT_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${CNT_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${CNT_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${CNT_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    CNT
    TARGET
    ${CNT_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${CNT_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${CNT}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${CNT}-lib DESTINATION ${${CNT}_LIB_INSTALL_DIR} EXPORT ${CNT}Targets)

# Install exported target definitions.
install(EXPORT ${CNT}Targets DESTINATION ${${CNT}_CONFIG_INSTALL_DIR})
