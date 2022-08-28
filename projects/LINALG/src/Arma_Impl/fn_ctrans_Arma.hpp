/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_ctrans_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CTRANS_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CTRANS_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    ctrans( const ET & aA )
    ->decltype( arma::trans(aA) )
    {
        return arma::trans(aA);
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CTRANS_ARMA_HPP_ */

