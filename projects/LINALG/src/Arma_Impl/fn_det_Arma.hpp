/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_det_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DET_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DET_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template<typename ET >
    auto
    det( const ET &  A)
    ->decltype( arma::det( A ))
    {
        return arma::det( A );
    }

    template<typename ET >
    auto
    det( ET &  A)
    ->decltype( arma::det( A ))
    {
        return arma::det( A );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DET_ARMA_HPP_ */

