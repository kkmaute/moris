/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_cond_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_COND_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_COND_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET>
    auto
    cond(ET const & aA)
    -> decltype( arma::cond( aA ) )
    {
        return arma::cond( aA );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_COND_ARMA_HPP_ */

