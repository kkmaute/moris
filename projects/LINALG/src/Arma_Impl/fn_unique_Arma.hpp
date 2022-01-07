/*
 * fn_unique_Arma.hpp
 *
 *  Created on: Sep 5, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_UNIQUE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_UNIQUE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template<typename ET, typename Matrix_Type>
    void
    unique( const                           ET &  aMatrix,
            moris::Matrix< Matrix_Type > & aUniqueMatrix )
    {
        aUniqueMatrix = unique( aMatrix );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_UNIQUE_ARMA_HPP_ */
