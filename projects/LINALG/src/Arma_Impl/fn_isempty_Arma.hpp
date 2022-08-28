/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isempty_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISEMPTY_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISEMPTY_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    isempty( ET const & aA )
    ->decltype( aA.is_empty() )
    {
        return aA.is_empty();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISEMPTY_ARMA_HPP_ */

