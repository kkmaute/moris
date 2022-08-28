/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_inv.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_INV_HPP_
#define PROJECTS_LINALG_SRC_FN_INV_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_inv_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_inv_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief Compute the inverse of Matrix A.
 *
 * @param[in] aA Matrix.
 *
 * @return The inverse of the Matrix aA.
 */
template< typename Matrix_Type >
auto
inv( const Matrix< Matrix_Type > & aA )
-> decltype( inv( aA.matrix_data() ) )
{
    MORIS_ASSERT( aA.n_rows() < 15, "For matrices larger than 15x15 use solve() instead of inv().\n");

    return inv( aA.matrix_data() );
}

}

#endif /* PROJECTS_LINALG_SRC_FN_INV_HPP_ */

