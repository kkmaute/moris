/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_find_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_FIND_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_FIND_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    //Find all indices
    template< typename ET >
    auto
    find( ET & aA )
    ->decltype( arma::find( aA ) )
    {
        return arma::find( aA );
    }

    // Find a specific number (ab) of indices
    template< typename ET >
    auto
    find( ET          const & aA,
          moris::uint const & ab )
    ->decltype( arma::find( aA, ab ) )
    {
        return ( arma::find( aA, ab ) );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_FIND_ARMA_HPP_ */

