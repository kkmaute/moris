# Distributed Linear Algebra Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if DLA has already been included
if(DEFINED DLA_CONFIGURED_ONCE)
    return()
endif()

set(DLA_CONFIGURED_ONCE "YES")

# Add DLA to the source directory list
list(APPEND MORIS_SRC_DIRS ${DLA})

# Include libraries needed by DLA
# PETSc and Trilinos; add later
# include(share/cmake/PETSc.cmake)
set(DLA_TPL_DEPENDENCIES
    "PETSc"
    #"trilinos"
    )

include(share/cmake/MatthewCMake/LNA_Depends.cmake)

list(APPEND DLA_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES} )
