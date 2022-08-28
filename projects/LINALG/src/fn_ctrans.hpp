/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_ctrans.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_CTRANS_HPP_
#define PROJECTS_LINALG_SRC_FN_CTRANS_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_ctrans_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_ctrans_Arma.hpp"
#endif

namespace moris
{
/*
 * @brief Computes the conjugate transpose of Matrix A. This functionality
 * is more costly than the 'trans' function, but may be necessary for complex arithmetic
 *
 * @param[in] A Matrix.
 *
 * @return The conjugate transpose of the Matrix A.
 *
 * Examples:
 * @include LNA/src/fn_ctrans/ctrans_real.inc
 * @include LNA/src/fn_ctrans/ctrans_complex.inc
 */
template< typename Matrix_Type >
auto
ctrans( const Matrix< Matrix_Type > & aA )
-> decltype( ctrans( aA.matrix_data() ) )
{
    return ctrans( aA.matrix_data() );
}

}

#endif /* PROJECTS_LINALG_SRC_FN_CTRANS_HPP_ */

