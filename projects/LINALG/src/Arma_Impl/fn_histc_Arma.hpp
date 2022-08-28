/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_histc_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_HISTC_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_HISTC_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET1, typename ET2>
    auto
    histc( ET1 const & aA,
           ET2 const & aB )
    ->decltype( arma::histc( aA, aB ) )
    {
        return ( arma::histc( aA, aB ) ) ;
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_HISTC_ARMA_HPP_ */

