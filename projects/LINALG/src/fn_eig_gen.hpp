/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_eig_gen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_EIG_GEN_HPP_
#define PROJECTS_LINALG_SRC_FN_EIG_GEN_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_eig_gen_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_eig_gen_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief Eigenvalue decomposition of dense general Matrix A.
 *
 * @param[in] aA Matrix.
 *
 * @param[out] eigenvalues The eigenvalues of Matrix aA.
 *
 * @param[out] eigenvectors The eigenvectors of Matrix aA.
 */
    template<typename Matrix_Type1, typename Matrix_Type2 >
    void
    eig_gen(      Matrix< Matrix_Type1  >  & aEigenvalues,
                  Matrix< Matrix_Type1  >  & aEigenvectors,
            const Matrix< Matrix_Type2  >  & aA )
    {
        eig_gen( aEigenvalues.matrix_data(), aEigenvectors.matrix_data(), aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_EIG_GEN_HPP_ */

