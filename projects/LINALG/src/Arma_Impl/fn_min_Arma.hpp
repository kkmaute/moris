/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_min_Arma.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"

namespace moris
{

    template< typename ETA, typename ETB >
    auto
    min( ETA&    aA,
            ETB& aB )
            -> decltype( arma::min( aA, aB ) )
    {
        return arma::min( aA, aB );
    }
}    // namespace moris
