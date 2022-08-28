/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_find_unique_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_FIND_UNIQUE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_FIND_UNIQUE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    //Find all indices
    template< typename ET >
    auto
    find_unique( ET & aA )
    ->decltype( arma::find_unique( aA ) )
    {
        return arma::find_unique( aA );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_FIND_UNIQUE_ARMA_HPP_ */

