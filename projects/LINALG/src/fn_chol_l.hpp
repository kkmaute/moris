/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_chol_l.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_CHOL_L_HPP_
#define PROJECTS_LINALG_SRC_FN_CHOL_L_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_chol_l_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_chol_l_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Cholesky decomposition of a Hermitian positive-definite matrix
     *
     *
     * @param[in]  aA Input matrix
     *
     * @return aL Lower triangle matrix
     *
     * This Cholesky decomposition uses only the diagonal and the lower triangle
     * of A to produce a lower triangular L so that
     * @f$ \mathbf{A}_{ij} = \mathbf{L}_{ik} \mathbf{L}_{jk} @f$
     * If A is not positive definite, an error message is printed.
     */
    template< typename Matrix_Type >
    auto
    chol_l(  const Matrix< Matrix_Type >& aA )
        -> decltype( chol_l( aA.matrix_data() ) )
    {
        return chol_l( aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_CHOL_L_HPP_ */

