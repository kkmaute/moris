# MAP Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MAP}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    )

# List header dependencies
set(MAP_INCLUDES
    ${LINALG}
    ${COM}
    ${DLA}
    ${NLA}
    ${MTK}
    ${HMR}
    ${MTK}/${MAP}
    ${FEM}/${INT}
    ${FEM}/${MDL}
    ${FEM}/${MSI}
    )
    
# Include depedency directories
foreach(MAP_INCLUDE ${MAP_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${MAP_INCLUDE}/src)
endforeach()
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAP Library

# List library source files
set( LIB_SOURCES
    cl_MTK_Mapper.cpp
    cl_MTK_Mapper_Arguments.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
    ${MTK}-lib
    ${COM}-lib
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_BOOST_LIBS}
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MAP}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${MAP}-lib ${LIB_DEPENDENCIES})
set_target_properties(${MAP}-lib PROPERTIES OUTPUT_NAME ${MAP})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mapper Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# project Executable

# List executable source files
set(EXE_SOURCES
    main.cpp )

set(EXE_DEPENDENCIES
    ${HMR}-lib
    ${GEN}-lib
    ${MAP}-lib
    ${MTK}-lib
    ${MDL}-lib # for L2
    ${INT}-lib # for L2
    ${MSI}-lib # for L2
    ${NLA}-lib # for L2
    ${DLA}-lib # for L2
    ${ALG}-lib #MTK
    ${TOL}-lib #MTK
#    ${XTK}-lib #MTK
#    ${STK}-lib #MTK
    ${TOL}-lib
    ${COM}-lib
#    ${LINALG}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    )
    
# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
#add_executable(${MAP}-exe ${EXE_SOURCES})
#target_link_libraries(${MAP}-exe ${EXE_DEPENDENCIES})
#set_target_properties(${MAP}-exe PROPERTIES OUTPUT_NAME mapper )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# map Config Files

# Build *Config and *ConfigVersion files
#include(CMakePackageConfigHelpers)

#configure_package_config_file(
#    ${MAP}Config.cmake.in
#    ${${MAP}_CONFIG_BUILD_DIR}/${MAP}Config.cmake
#    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${MAP}
#    PATH_VARS LIB_INSTALL_DIR )

#write_basic_package_version_file(
#    ${${MAP}_CONFIG_BUILD_DIR}/${MAP}ConfigVersion.cmake
#    COMPATIBILITY AnyNewerVersion )
    
    