#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

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
list(APPEND MORIS_SOURCE_DIRS 
	${MRS}/${ASR} 
	${MRS}/${CHR} 
	${MRS}/${CON} 
	${MRS}/${COR} 
	${MRS}/${EXC} 
	${MRS}/${IOS} )

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
    "boost"
    "mpi" )

 if(USE_GPERFTOOLS) #> TEMPORARY SOLUTION
	list(APPEND ASR_TPL_DEPENDENCIES "gperftools")
	list(APPEND CHR_TPL_DEPENDENCIES "gperftools")
	list(APPEND CON_TPL_DEPENDENCIES "gperftools")
	list(APPEND COR_TPL_DEPENDENCIES "gperftools")
	list(APPEND EXC_TPL_DEPENDENCIES "gperftools")
	list(APPEND IOS_TPL_DEPENDENCIES "gperftools")
endif()

