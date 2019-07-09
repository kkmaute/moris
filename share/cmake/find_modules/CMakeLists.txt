# Install moris custom Find modules

set(MORIS_FIND_MODULES
	FindACML.cmake
	FindARPACK.cmake
	FindEigen.cmake
	FindEigen3.cmake
	FindGCMMA.cmake
	FindHDF5.cmake
	FindMKL.cmake
	FindPETSc.cmake
	FindSNOPT.cmake
	FindSuiteSparse.cmake
	FindSuperLU_DIST.cmake
	FindSuperLU.cmake
	)

install(FILES ${MORIS_FIND_MODULES}
	DESTINATION ${CONFIG_INSTALL_DIR}/find_modules )
