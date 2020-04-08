# MORIS mains -------------------------------------------------------------
# -------------------------------------------------------------------------

# Note: These execiutables are to be replaced by project-specific mains
#       ...eventaually.

# Third party libraries
list(REMOVE_DUPLICATES MAIN_TPL_DEPENDENCIES)
foreach(MAIN_TPL ${MAIN_TPL_DEPENDENCIES})
    include(${MORIS_TPL_DIR}/${MAIN_TPL}.cmake)
endforeach()

# Set executable output location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# List mains
set(MAINS
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
    ${WRK}-lib
    ${MORIS_PETSC_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
#     ${MORIS_MPI_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_LDLIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
    )

# Create and link main executables
foreach(MAIN ${MAINS})
    string(REPLACE ".cpp" "-${EXECUTABLE}" EXE_NAME ${MAIN})
    add_executable(moris ${MAIN})
    target_link_libraries(moris ${MAIN_DEPENDENCIES})
    list(APPEND MAIN_EXECUTABLES moris)
endforeach()

# Install executables
install(TARGETS ${MAIN_EXECUTABLES}
    DESTINATION ${EXE_INSTALL_DIR} )
# -------------------------------------------------------------------------
# Making of input file shared object -------------------------------------------------------------
# -------------------------------------------------------------------------

## create dynamically linked .so files
set(HARDCODED_INPUT_FILE input_file)

set(SO_TPLS
"trilinos"
${ARMADILLO_EIGEN}
)


SET(SO_INCLUDES 
	${Boost_INCLUDE_DIRS}
	${Trilinos_INCLUDE_DIRS}
	${MORIS_PACKAGE_DIR}/ALG/src
	${MORIS_PACKAGE_DIR}/FEM/INT/src 
	${MORIS_PACKAGE_DIR}/FEM/MSI/src 
	${MORIS_PACKAGE_DIR}/FEM/VIS/src 
	${MORIS_PACKAGE_DIR}/GEN/GEN_CORE/src 
	${MORIS_PACKAGE_DIR}/SOL/DLA/src 
	${MORIS_PACKAGE_DIR}/SOL/TSA/src
	${MORIS_PACKAGE_DIR}/SOL/NLA/src 
	${MORIS_PACKAGE_DIR}/SOL/SOL_CORE/src 
	${MORIS_PACKAGE_DIR}/LINALG/src
	${MORIS_PACKAGE_DIR}/LINALG/src/${LINALG_IMPLEMENTATION_INCLUDES}
	${MORIS_PACKAGE_DIR}/COM/src 
	${MORIS_PACKAGE_DIR}/MTK/src
	${MORIS_PACKAGE_DIR}/HMR/src
	${MORIS_PACKAGE_DIR}/XTK/src/assert                                                                                                                                                                                                                                                                   
    ${MORIS_PACKAGE_DIR}/XTK/src/containers                                                                                                                                                                                                                                                               
    ${MORIS_PACKAGE_DIR}/XTK/src/core                                                                                                                                                                                                                                                                 
    ${MORIS_PACKAGE_DIR}/XTK/src/geomeng                                                                                                                                                                                                                                                                  
    ${MORIS_PACKAGE_DIR}/XTK/src/geometry                                                                                                                                                                                                                                                                 
    ${MORIS_PACKAGE_DIR}/XTK/src/ios                                                                                                                                                                                                                                                                      
    ${MORIS_PACKAGE_DIR}/XTK/src/linalg                                                                                                                                                                                                                                                                  
    ${MORIS_PACKAGE_DIR}/XTK/src/mesh                                                                                                                                                                                                                                                                     
    ${MORIS_PACKAGE_DIR}/XTK/src/tools                                                                                                                                                                                                                                                                    
    ${MORIS_PACKAGE_DIR}/XTK/src/topology                                                                                                                                                                                                                                                                 
    ${MORIS_PACKAGE_DIR}/XTK/src/xtk      
	${MORIS_PACKAGE_DIR}/PRM/src 
    ${MORIS_PACKAGE_DIR}/MRS/COR/src)
    
    message(INFO "\n\n\n\n\nso_includes = ${SO_INCLUDES}")
    
foreach(TPL ${SO_TPLS})
    string(TOUPPER ${TPL} TPL)
    list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
endforeach()

dynamic_link_input(${HARDCODED_INPUT_FILE} ${HARDCODED_INPUT_FILE}.cpp ${HARDCODED_INPUT_FILE}.so  ${SO_INCLUDES})

