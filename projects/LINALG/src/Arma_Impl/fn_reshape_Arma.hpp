/*
 * fn_reshape_Arma.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename Matrix_Type >
    arma::Mat< Matrix_Type >
    reshape(
            const arma::Mat< Matrix_Type > & aA,
            const size_t                   & aB,
            const size_t                   & aC )
    {
        return arma::reshape(aA, aB, aC);;
    }
} /* moris namespace */

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_ */
