/*
 * fn_find_unique.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
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
    template< typename Type, typename Matrix_Type >
    auto
    find_unique( const Matrix< Type, Matrix_Type > & aA)
    -> decltype( find_unique( aA.matrix_data()) )
    {
    	return find_unique( aA.matrix_data() );
    }
}
#endif /* PROJECTS_LINALG_SRC_FN_FIND_UNIQUE_HPP_ */
