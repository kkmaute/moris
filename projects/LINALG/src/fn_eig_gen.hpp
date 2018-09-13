/*
 * fn_eig_gen.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_FN_EIG_GEN_HPP_
#define PROJECTS_LINALG_SRC_FN_EIG_GEN_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_eig_gen_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_eig_gen_Arma.hpp"
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
    template< typename Type1, typename Matrix_Type1, typename Type2, typename Matrix_Type2 >
    void
    eig_gen(      Matrix< Type1, Matrix_Type1  >  & aEigenvalues,
                  Matrix< Type1, Matrix_Type1  >  & aEigenvectors,
            const Matrix< Type2, Matrix_Type2  >  & aA )
    {
        eig_gen( aEigenvalues.matrix_data(), aEigenvectors.matrix_data(), aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_EIG_GEN_HPP_ */
