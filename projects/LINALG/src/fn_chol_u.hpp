/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_chol_u.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_CHOL_U_HPP_
#define PROJECTS_LINALG_SRC_FN_CHOL_U_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_chol_u_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_chol_u_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Cholesky decomposition of a Hermitian positive-definite matrix
     *
     *
     * @param[in]  aA Input Matrix
     *
     * @return aU Upper triangle matrix
     * This Cholesky decomposition uses only the diagonal and the upper triangle
     * of A to produce an upper triangular U so that
     * @f$ \mathbf{A}_{ij} = \mathbf{U}_{ki} \mathbf{U}_{kj} @f$
     * If A is not positive definite, an error message is printed.
     *
     */
    template< typename Matrix_Type >
    auto
    chol_u(  const Matrix< Matrix_Type >& aA )
        -> decltype( chol_u( aA.matrix_data() ) )
    {
        return chol_u( aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_CHOL_U_HPP_ */

