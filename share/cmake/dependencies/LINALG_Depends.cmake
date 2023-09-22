#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

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

# Third party libraries directly needed by LINALG
set(LINALG_TPL_DEPENDENCIES
	${ARMADILLO_EIGEN}
	${ACML_LAPACK_MKL_OPENBLAS}
	"superlu"
	"arpack"
	)

if(${ARMADILLO_EIGEN} STREQUAL "armadillo")
    list(APPEND LINALG_TPL_DEPENDENCIES
        #"viennacl"
        )
elseif(${ARMADILLO_EIGEN} STREQUAL "eigen")
    list(APPEND LINALG_TPL_DEPENDENCIES
        "suitesparse"
        )
endif()

if(USE_GPERFTOOLS) #> TEMPORARY SOLUTION
	list(APPEND LINALG_TPL_DEPENDENCIES "gperftools")
endif()

foreach(TPL ${LINALG_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
endforeach()

# List moris projects directly needed by PROJ
set(LINALG_MORIS_DEPENDENCIES
    ${COM}
    ${ALG}
    )

foreach(MORIS_DEPENDENCY ${LINALG_MORIS_DEPENDENCIES})
    # Include moris projects directly needed by LINALG
    include(${MORIS_DEPENDS_DIR}/${MORIS_DEPENDENCY}_Depends.cmake)

    # Include third party libraries indirectly needed by LINALG
    #list(APPEND LINALG_TPL_DEPENDENCIES ${${MORIS_DEPENDENCY}_TPL_DEPENDENCIES} )
endforeach()

list(REMOVE_DUPLICATES LINALG_TPL_DEPENDENCIES)

# Linear Algebra implementations
if(MORIS_USE_EIGEN)
	set(LINALG_IMPLEMENTATION_INCLUDES Eigen_Impl )
endif()

if(MORIS_USE_ARMA)
	set(LINALG_IMPLEMENTATION_INCLUDES Arma_Impl )
endif()



