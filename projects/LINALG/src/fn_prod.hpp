/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_prod.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_PROD_HPP_
#define PROJECTS_LINALG_SRC_FN_PROD_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_prod_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_prod_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief For vector V, return the product of all elements
     * For matrix M, return the product of elements in each column
     *
     *@param[in] aA A matrix or vector
     *
     *
     */
    template< typename Matrix_Type >
    auto
    prod( Matrix< Matrix_Type > const & aA )
    -> decltype( prod(aA.matrix_data()) ) const
    {
        return prod(aA.matrix_data());
    }
}
#endif /* PROJECTS_LINALG_SRC_FN_PROD_HPP_ */

