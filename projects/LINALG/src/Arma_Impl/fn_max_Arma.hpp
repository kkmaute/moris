/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_max_Arma.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ETA, typename ETB >
    auto
    max( ETA &    aA,
              ETB & aB )
            -> decltype( arma::max( aA, aB ) )
    {
        return arma::max( aA, aB );
    }
}    // namespace moris
