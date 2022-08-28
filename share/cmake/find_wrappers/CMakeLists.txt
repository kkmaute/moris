#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

set(MORIS_FIND_WRAPPERS
	acml_new.cmake
	#acml.cmake
	armadillo_new.cmake
	#armadillo.cmake
	arpack_new.cmake
	#arpack.cmake
	boost_new.cmake
	#boost.cmake
	eigen_new.cmake
	#eigen.cmake
	gcmma_new.cmake
	#gcmma.cmake
	hdf5_new.cmake
	lapack_new.cmake
	#lapack.cmake
	mkl_new.cmake
	#mkl.cmake
	mpi_new.cmake
	#mpi.cmake
    openblas_new.cmake
	petsc_new.cmake
	#PETSc.cmake
	snopt_new.cmake
	#snopt.cmake
	suitesparse_new.cmake
	#suitesparse.cmake
	superlu_new.cmake
	#superlu.cmake
	trilinos_new.cmake
	#trilinos.cmake
	viennacl_new.cmake
	#viennacl.cmake
	)

install(FILES ${MORIS_FIND_WRAPPERS}
	DESTINATION ${CONFIG_INSTALL_DIR}/find_wrappers )
