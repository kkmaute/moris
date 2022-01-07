/*
 * fn_norm_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_NORM_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_NORM_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template<typename ET >
    auto
    norm( ET &  A)
        ->decltype( arma::norm( A, 2 ) )
    {
        return arma::norm( A, 2 );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_NORM_ARMA_HPP_ */
