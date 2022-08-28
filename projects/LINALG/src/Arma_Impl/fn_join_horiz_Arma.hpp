/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_join_horiz_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_JOIN_HORIZ_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_JOIN_HORIZ_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    join_horiz( const ET          & aA,
                const ET          & aB)
    -> decltype( arma::join_horiz( aA, aB ) )
    {
        return arma::join_horiz( aA, aB );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_JOIN_HORIZ_ARMA_HPP_ */

