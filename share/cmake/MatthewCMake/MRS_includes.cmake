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
