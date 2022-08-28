/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_inv_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_INV_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_INV_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    inv( const ET & aA )
    ->decltype( arma::inv(aA) )
    {
        return arma::inv(aA);
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_INV_ARMA_HPP_ */

