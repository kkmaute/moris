/*
 * fn_unique.hpp
 *
 *  Created on: Sep 5, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_FN_UNIQUE_HPP_
#define PROJECTS_LINALG_SRC_FN_UNIQUE_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_unique_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_unique_Arma.hpp"
#endif

namespace moris
{
    template< typename Type, typename Matrix_Type >
    auto
    unique( const Matrix< Type, Matrix_Type > & aMatrix )
    -> decltype( unique( aMatrix.matrix_data() ) )
    {
        return unique( aMatrix.matrix_data() );
    }
}



#endif /* PROJECTS_LINALG_SRC_FN_UNIQUE_HPP_ */
