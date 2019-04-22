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
#include "fn_norm_Arma.hpp"
#endif


namespace moris
{
    /**
     * @brief Calculate the L2 norm of a matrix.
     *
     *@param[in] aA A given matrix
     *
     *
     */
    template< typename Matrix_Type >
    auto
    norm( const Matrix< Matrix_Type > & aA )
        -> decltype( norm( aA.matrix_data() ) )
    {
        return norm( aA.matrix_data() );
    }

    template< typename Matrix_Type >
    auto
    norm( Matrix< Matrix_Type > & aA )
        -> decltype( norm( aA.matrix_data() ) )
    {
        return norm( aA.matrix_data() );
    }

}



#endif /* PROJECTS_LINALG_SRC_FN_NORM_HPP_ */
