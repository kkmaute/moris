/*
 * fn_norm.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_FN_NORM_HPP_
#define PROJECTS_LINALG_SRC_FN_NORM_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_norm_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_norm_Arma.hpp"
#endif


namespace moris
{
    /**
     * @brief Calculate the l2 norm of a matrix.
     *
     *@param[in] aA A given matrix
     *
     *
     */
    template< typename Type, typename Matrix_Type >
    auto
    norm( Matrix< Type, Matrix_Type  > const & aA )
    -> decltype( norm( aA.matrix_data() ) )
    {
        return norm( aA.matrix_data() );
    }
}



#endif /* PROJECTS_LINALG_SRC_FN_NORM_HPP_ */
