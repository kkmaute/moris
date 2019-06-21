# MTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MTK has already been included
if(DEFINED MTK_CONFIGURED_ONCE)
    return()
endif()

set(MTK_CONFIGURED_ONCE "YES")

# Add MTK to the header directory list
list(APPEND MORIS_SOURCE_DIRS ${MTK})

# Third party libraries needed by MTK-lib
set(MTK_LIB_TPL_DEPENDENCIES
    "trilinos"
    #"arpack" # needed for linking order
    )

# Additional third party libraraies needed by MTK-test
set(MTK_TEST_TPL_DEPENDENCIES
    ""
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)

# for test
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

#set(MTK_PROJ ${MORIS_PACKAGE_DIR}/${MTK}/src/stk_impl)


#include_directories(${MTK_PROJ})

# Include third party libraries indirectly needed by MTK
#list(APPEND MTK_TPL_DEPENDENCIES
#    )
    
#list(REMOVE_DUPLICATES MTK_TPL_DEPENDENCIES)