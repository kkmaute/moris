/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_rank_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RANK_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RANK_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET>
    auto
    rank(ET const & aA)
    -> decltype( arma::rank( aA ) )
    {
        return arma::rank( aA );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RANK_ARMA_HPP_ */

