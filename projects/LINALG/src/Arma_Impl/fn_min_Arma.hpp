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
    // Reuse native overloads for both reductions and coefficient-wise extrema.
    using arma::min;
}    // namespace moris
