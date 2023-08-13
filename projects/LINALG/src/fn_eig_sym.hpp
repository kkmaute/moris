/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_eig_sym.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_EIG_SYM_HPP_
#define PROJECTS_LINALG_SRC_FN_EIG_SYM_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_eig_sym_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_eig_sym_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief Eigenvalue decomposition of dense symmetric/Hermitian Matrix A.
 *
 * @param[in] aA Matrix.
 *
 * @param[out] eigenvalues The eigenvalues of Matrix aA.
 *
 * @param[out] eigenvectors The eigenvectors of Matrix aA.
 */
    template< typename Matrix_Type >
    void
    eig_sym(      Matrix< Matrix_Type >  & aEigenvalues,
                  Matrix< Matrix_Type >  & aEigenvectors,
            const Matrix< Matrix_Type >  & aA )
    {
        eig_sym( aEigenvalues.matrix_data(), aEigenvectors.matrix_data(), aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_EIG_SYM_HPP_ */

