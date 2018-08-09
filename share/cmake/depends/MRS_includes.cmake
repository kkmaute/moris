# MRS includes ------------------------------------------------------------
# -------------------------------------------------------------------------

# Include core MORIS directories needed by nearly everything
include_directories(
    ${MORIS_PACKAGE_DIR}/MRS/ASR/src
    ${MORIS_PACKAGE_DIR}/MRS/CHR/src
    ${MORIS_PACKAGE_DIR}/MRS/CON/src
    ${MORIS_PACKAGE_DIR}/MRS/COR/src
    ${MORIS_PACKAGE_DIR}/MRS/EXC/src )

# Add MRS to the header directory list
list(APPEND MORIS_SOURCE_DIRS ${MRS})

#chronos and containers need some tpls for their tests
include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake)

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
