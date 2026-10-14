# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

set(MORIS_FIND_WRAPPERS
    acml_new.cmake
    armadillo_new.cmake
    arborx_new.cmake
    arpack_new.cmake
    boost_new.cmake
    eigen_new.cmake
    gperftools_new.cmake
    gcmma_new.cmake
    hdf5_new.cmake
    lbfgsb_new.cmake
    lapack_new.cmake
    mkl_new.cmake
    mpi_new.cmake
    openblas_new.cmake
    petsc_new.cmake
    slepc_new.cmake
    snopt_new.cmake
    suitesparse_new.cmake
    superlu_new.cmake
    viennacl_new.cmake
)

install(FILES ${MORIS_FIND_WRAPPERS} DESTINATION ${CONFIG_INSTALL_DIR}/find_wrappers)

set(MORIS_TRILINOS_WRAPPER trilinos_new_install.cmake)

install(FILES ${MORIS_TRILINOS_WRAPPER} DESTINATION ${CONFIG_INSTALL_DIR}/find_wrappers
        RENAME trilinos_new.cmake
)
