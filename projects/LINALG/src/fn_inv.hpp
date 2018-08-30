/*
 * fn_inv.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_FN_INV_HPP_
#define PROJECTS_LINALG_SRC_FN_INV_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_inv_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_inv_Arma.hpp"
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
template< typename Type, typename Matrix_Type>
auto
inv( const Matrix< Type, Matrix_Type > & aA )
-> decltype( inv( aA.matrix_data() ) )
{
    return inv( aA.matrix_data() );
}

}



#endif /* PROJECTS_LINALG_SRC_FN_INV_HPP_ */
