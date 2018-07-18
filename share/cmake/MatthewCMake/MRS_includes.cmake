# MRS includes ------------------------------------------------------------
# -------------------------------------------------------------------------

# Include core MORIS directories needed by nearly everything
include_directories(
    ${PROJECT_SOURCE_DIR}/projects/MRS/ASR/src
    ${PROJECT_SOURCE_DIR}/projects/MRS/CHR/src
    ${PROJECT_SOURCE_DIR}/projects/MRS/CON/src
    ${PROJECT_SOURCE_DIR}/projects/MRS/COR/src
    ${PROJECT_SOURCE_DIR}/projects/MRS/EXC/src )

# Add MRS to the header directory list
list(APPEND MORIS_HEADER_DIRS ${MRS})

#chronos and containers need some tpls for their tests
include(share/cmake/MatthewCMake/LNA_Depends.cmake)

set(ASR_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES} )

set(CHR_TPL_DEPENDENCIES
    "boost"
    ${LNA_TPL_DEPENDENCIES} )

set(CON_TPL_DEPENDENCIES
    "boost"
    ${LNA_TPL_DEPENDENCIES} )

set(EXC_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES} )
