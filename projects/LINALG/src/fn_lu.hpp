/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_lu.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_LU_HPP_
#define PROJECTS_LINALG_SRC_FN_LU_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_lu_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_lu_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief LU decomposition with partial pivoting.
 *
 *@param[out] aL Lower-triangular output matrix
 *@param[out] aU Upper-triangular output matrix
 *@param[out] aP Permutation matrix
 *@param[in]  aA Input matrix
 *
 * The LU decomposition of a matrix aA is such that \f$ P*A = L*U \f$ where L is a lower
 * triangular matrix, U is an upper triangular matrix, and P is a permutation matrix.
 *
 */
    template< typename Matrix_Type >
    void
    lu(
            Matrix< Matrix_Type >       & aL,
            Matrix< Matrix_Type >       & aU,
            Matrix< Matrix_Type >       & aP,
            Matrix< Matrix_Type > const & aA )
    {
        lu( aL.matrix_data(),
            aU.matrix_data(),
            aP.matrix_data(),
            aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_LU_HPP_ */

