# LINALG Dependencies -------------------------------------------------------
# -------------------------------------------------------------------------

# Check if LINALG has already been included
if(DEFINED LINALG_CONFIGURED_ONCE)
    return()
endif()

set(LINALG_CONFIGURED_ONCE "YES")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Handle Dependencies

# Add LINALG to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${LINALG})

# Include third party libraries directly needed by LINALG
set(LINALG_TPL_DEPENDENCIES
    ${ARMADILLO_EIGEN}
    "superlu"
    "trilinos")
    
# List moris projects directly needed by PROJ
set(LINALG_MORIS_DEPENDENCIES
    ${COM}
    ${ALG} )

foreach(MORIS_DEPENDENCY ${LINALG_MORIS_DEPENDENCIES})
    # Include moris projects directly needed by LINALG
    include(${MORIS_CMAKE_DIR}/depends/${MORIS_DEPENDENCY}_Depends.cmake)

    # Include third party libraries indirectly needed by LINALG
    list(APPEND LINALG_TPL_DEPENDENCIES ${${MORIS_DEPENDENCY}_TPL_DEPENDENCIES} )
endforeach()


# Linear Algebra implementations
set(LINALG_IMPLEMENTATION_INCLUDES ${MORIS_PACKAGE_DIR}/${LINALG}/src)

if(MORIS_USE_EIGEN)
	list(APPEND LINALG_IMPLEMENTATION_INCLUDES ${MORIS_PACKAGE_DIR}/${LINALG}/src/Eigen_Impl )
endif()

if(MORIS_USE_ARMA)
	list(APPEND LINALG_IMPLEMENTATION_INCLUDES ${MORIS_PACKAGE_DIR}/${LINALG}/src/Arma_Impl )
endif()


