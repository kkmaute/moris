/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include "fn_reshape_Eigen.hpp"

namespace moris
{
    template< typename E >
    auto vectorize( const Eigen::MatrixBase< E >& aA, bool aRow = false )
    {
        const Eigen::Index tRows = aRow ? 1 : aA.size();
        const Eigen::Index tCols = aRow ? aA.size() : 1;
        return linalg_internal::eigen_plain< E >::NullaryExpr( tRows, tCols, linalg_internal::eigen_reshape< E >{ aA.derived(), tRows, aRow } );
    }
}    // namespace moris
