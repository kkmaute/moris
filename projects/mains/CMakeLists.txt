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
    add_executable(${EXE_NAME} ${MAIN})
    target_link_libraries(${EXE_NAME} ${MAIN_DEPENDENCIES})
    list(APPEND MAIN_EXECUTABLES ${EXE_NAME})
endforeach()

# Install executables
install(TARGETS ${MAIN_EXECUTABLES}
    DESTINATION ${EXE_INSTALL_DIR} )

