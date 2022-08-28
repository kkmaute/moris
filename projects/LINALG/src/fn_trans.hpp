/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_trans.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_TRANS_HPP_
#define PROJECTS_LINALG_SRC_FN_TRANS_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_trans_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_trans_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief Computes transpose of Matrix A.
 *
 * @param[in] A Matrix.
 *
 * @return The transpose of the Matrix A.
 *
 * @note In case of complex numbers, no conjugate is taken.
 *
 * Example:
 * @include LNA/src/fn_trans/trans_real.inc
 * @include LNA/src/fn_trans/trans_complex.inc
 */
template< typename Matrix_Type >
auto
trans( const Matrix< Matrix_Type > & A )
-> decltype( linalg_internal::trans(A.matrix_data()) )
{
    return linalg_internal::trans(A.matrix_data());
}

template< typename Matrix_Type >
auto
trans( Matrix< Matrix_Type > & A )
-> decltype( linalg_internal::trans(A.matrix_data()) )
{
    return linalg_internal::trans(A.matrix_data());
}

}

#endif /* PROJECTS_LINALG_SRC_FN_TRANS_HPP_ */

