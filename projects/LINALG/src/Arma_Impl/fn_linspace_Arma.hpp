/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_linspace_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LINSPACE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LINSPACE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename T >
    auto
    linspace_base( T             const & aStart,
                   T             const & aEnd,
                   moris::size_t const & aN )
    ->decltype( arma::linspace< arma::Mat< T > >( aStart, aEnd, aN ) )
    {
        return arma::linspace< arma::Mat< T > >( aStart, aEnd, aN );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LINSPACE_ARMA_HPP_ */

