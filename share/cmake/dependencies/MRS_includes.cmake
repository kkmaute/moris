# MRS includes ------------------------------------------------------------
# -------------------------------------------------------------------------

# Include core MORIS directories needed by nearly everything
include_directories(
    ${MORIS_PACKAGE_DIR}/MRS/ASR/src
    ${MORIS_PACKAGE_DIR}/MRS/CHR/src
    ${MORIS_PACKAGE_DIR}/MRS/CON/src
    ${MORIS_PACKAGE_DIR}/MRS/COR/src
    ${MORIS_PACKAGE_DIR}/MRS/EXC/src
    ${MORIS_PACKAGE_DIR}/MRS/IOS/src )

# Add MRS to the header directory list
list(APPEND MORIS_SOURCE_DIRS ${MRS})

set(ASR_TPL_DEPENDENCIES
    "" )

set(CHR_TPL_DEPENDENCIES
    "boost" )

set(CON_TPL_DEPENDENCIES
    "boost" )

set(EXC_TPL_DEPENDENCIES
    "" )
    
set(IOS_TPL_DEPENDENCIES
    "hdf5"
    "boost" )

 
    