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
#include "fn_unique_Arma.hpp"
#endif

namespace moris
{
    template< typename Matrix_Type >
    void
    unique( const Matrix< Matrix_Type > & aMatrix,
                  Matrix< Matrix_Type > & aUniqueMatrix )
    {
        unique( aMatrix.matrix_data(), aUniqueMatrix );
    }
}



#endif /* PROJECTS_LINALG_SRC_FN_UNIQUE_HPP_ */
