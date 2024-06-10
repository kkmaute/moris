/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for traceails.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_trace_Arma.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"

namespace moris
{
    template<typename ET >
    auto
    trace( const ET &  A)
    ->decltype( arma::trace( A ))
    {
        return arma::trace( A );
    }

    template<typename ET >
    auto
    trace( ET &  A)
    ->decltype( arma::trace( A ))
    {
        return arma::trace( A );
    }
}



