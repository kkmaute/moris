#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# MORIS main -------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party libraries
list(REMOVE_DUPLICATES MAIN_TPL_DEPENDENCIES)
foreach(MAIN_TPL ${MAIN_TPL_DEPENDENCIES})
    include(${MORIS_TPL_DIR}/${MAIN_TPL}_new.cmake)
endforeach()

# Set executable output location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# main moris file
set(MAIN
    main.cpp )

# List include directories
set(MAIN_INCLUDES
    ${COM}
    ${IOS}
    ${WRK}
    ${LINALG}
        )

foreach(MAIN_INCLUDE ${MAIN_INCLUDES})
    include_directories(../${MAIN_INCLUDE}/src)
endforeach()

# Additonal (optional) includes
include_directories(../STK/src/Hierarchical)

# List executable dependencies
set(MAIN_DEPENDENCIES
    ${HMR}-lib
    ${COM}-lib
    ${IOS}-lib
    ${GEN_MAIN}-lib
    ${WRK}-lib
    ${MORIS_PETSC_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_LDLIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
    )

# ------------------------------------------------------------------------

# Create and link main executables
string(REPLACE ".cpp" "-${EXECUTABLE}" EXE_NAME ${MAIN})
add_executable(moris ${MAIN})
target_link_libraries(moris ${MAIN_DEPENDENCIES})
list(APPEND MAIN_EXECUTABLES moris)


# Install executables
install(TARGETS ${MAIN_EXECUTABLES}
    DESTINATION ${EXE_INSTALL_DIR} )

# -------------------------------------------------------------------------
# Making of input file shared object --------------------------------------
# -------------------------------------------------------------------------

## create dynamically linked .so files
set(HARDCODED_INPUT_FILE input_file)

set(SO_TPLS
"trilinos"
${ARMADILLO_EIGEN}
)

# List header dependencies in /FEM/INT/src                                                                                                                                                                                                                                                    
set(INT_INCLUDES
    "CORE"
    "ELEM"
    "IP"
    "MM"
    "CM"
    "SP"
    "IWG"
    "IQI")

set(INT_INTERNAL_INCLUDES)

# Include depedency directories
foreach(INT_INCLUDES ${INT_INCLUDES})
    list(APPEND INT_INTERNAL_INCLUDES "${MORIS_PACKAGE_DIR}/FEM/INT/src/${INT_INCLUDES}")
endforeach()

# List header dependencies in /MTK/src                                                                                                                                                                                                                                                    
set(MTK_INCLUDES
    "IP"
    "IG")

set(MTK_INTERNAL_INCLUDES)

# Include depedency directories
foreach(MTK_INCLUDES ${MTK_INCLUDES})
    list(APPEND MTK_INTERNAL_INCLUDES "${MORIS_PACKAGE_DIR}/MTK/src/${MTK_INCLUDES}")
endforeach()

SET(SO_INCLUDES
    ${Boost_INCLUDE_DIRS}
    ${Trilinos_INCLUDE_DIRS}
    ${MORIS_PACKAGE_DIR}/ALG/src
    ${INT_INTERNAL_INCLUDES}
    ${MTK_INTERNAL_INCLUDES} 
    ${MORIS_PACKAGE_DIR}/PRM/src/ 
    ${MORIS_PACKAGE_DIR}/FEM/MSI/src 
    ${MORIS_PACKAGE_DIR}/FEM/VIS/src 
    ${MORIS_PACKAGE_DIR}/GEN/GEN_CORE/src 
    ${MORIS_PACKAGE_DIR}/HMR/src 
    ${MORIS_PACKAGE_DIR}/SOL/DLA/src 
    ${MORIS_PACKAGE_DIR}/SOL/TSA/src
    ${MORIS_PACKAGE_DIR}/SOL/NLA/src 
    ${MORIS_PACKAGE_DIR}/SOL/SOL_CORE/src 
    ${MORIS_PACKAGE_DIR}/LINALG/src
    ${MORIS_PACKAGE_DIR}/LINALG/src/${LINALG_IMPLEMENTATION_INCLUDES}
    ${MORIS_PACKAGE_DIR}/COM/src 
    ${MORIS_PACKAGE_DIR}/PRM/ENM/src
    ${MORIS_PACKAGE_DIR}/MTK/src
    ${MORIS_PACKAGE_DIR}/HMR/src
    ${MORIS_PACKAGE_DIR}/XTK/src/
    ${MORIS_PACKAGE_DIR}/MRS/COR/src)

# make SO_INCLUDES available as global list
set_property(GLOBAL PROPERTY INT_SRC_LIST "${SO_INCLUDES}")
set_property(GLOBAL PROPERTY MTK_SRC_LIST "${SO_INCLUDES}")
    
foreach(TPL ${SO_TPLS})
    string(TOUPPER ${TPL} TPL)
    list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
endforeach()

dynamic_link_input("shared_object_file" ${HARDCODED_INPUT_FILE} ${HARDCODED_INPUT_FILE}.cpp  ${SO_INCLUDES})


