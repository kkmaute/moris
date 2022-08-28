/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_norm_Arma.hpp
 *
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

