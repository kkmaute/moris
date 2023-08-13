#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Moris Packages ----------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory variables
foreach(MORIS_SOURCE_DIR ${MORIS_SOURCE_DIRS})
    string(REGEX REPLACE "([^\\/]*)\\/" "" DIR_NAME ${MORIS_SOURCE_DIR})
    
    # Build directories
    set(${DIR_NAME}_CONFIG_BUILD_DIR "${CONFIG_BUILD_DIR}/${MORIS_SOURCE_DIR}")
    
    # Install directories
    set(${DIR_NAME}_CONFIG_INSTALL_DIR "${CONFIG_INSTALL_DIR}/${MORIS_SOURCE_DIR}")
    set(${DIR_NAME}_HEADER_INSTALL_DIR "${HEADER_INSTALL_DIR}/${MORIS_SOURCE_DIR}")
    set(${DIR_NAME}_LIB_INSTALL_DIR "${LIB_INSTALL_DIR}/${MORIS_SOURCE_DIR}")
    set(${DIR_NAME}_EXE_INSTALL_DIR "${EXE_INSTALL_DIR}/${MORIS_SOURCE_DIR}")
endforeach()

#GAH - Add large library "libmorisCombinedLib"
# Add to it in each subdir
add_library(morisCombinedLib ${LIB_LINK_MODE})
target_compile_definitions(morisCombinedLib
      PUBLIC ${MORIS_DEFINITIONS})


add_subdirectory(GEN)
add_subdirectory(FEM)
add_subdirectory(SOL)
add_subdirectory(WRK)
add_subdirectory(XTK)
add_subdirectory(TOL)
add_subdirectory(MTK)
add_subdirectory(MIG)
add_subdirectory(ALG)
add_subdirectory(COM)
#add_subdirectory(EXA)
add_subdirectory(HMR)
add_subdirectory(LINALG)
add_subdirectory(MOD)
add_subdirectory(MRS)
add_subdirectory(OPT)
add_subdirectory(PRM)

# Install big library and executable
install(TARGETS morisCombinedLib
#    EXPORT ${COM}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR} )


