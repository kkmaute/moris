/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_find_unique.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_FIND_UNIQUE_HPP_
#define PROJECTS_LINALG_SRC_FN_FIND_UNIQUE_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_find_unique_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_find_unique_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Return a column vector containing the indices of unique elements of the input vector ( first position of the specific value )
     *
     * @param[in] aMat     Vector.
     *
     * @return  Column vector containing the indices of unique elements
     */
    template< typename Matrix_Type >
    auto
    find_unique( const Matrix< Matrix_Type > & aA)
    -> decltype( find_unique( aA.matrix_data()) )
    {
    	return find_unique( aA.matrix_data() );
    }
}
#endif /* PROJECTS_LINALG_SRC_FN_FIND_UNIQUE_HPP_ */

