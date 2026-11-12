/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once
#include <catch.hpp>
#include "cl_Matrix.hpp"
#include <complex>

namespace moris::linalg_test
{
    template< typename Actual, typename Expected >
    void check_matrix( const Matrix< Actual >& aActual, const Matrix< Expected >& aExpected, real aTolerance = 1e-12 )
    {
        REQUIRE( aActual.n_rows() == aExpected.n_rows() );
        REQUIRE( aActual.n_cols() == aExpected.n_cols() );
        for ( uint j = 0; j < aActual.n_cols(); ++j )
        {
            for ( uint i = 0; i < aActual.n_rows(); ++i )
            {
                INFO( "row=" << i << " column=" << j << " actual=" << aActual( i, j ) << " expected=" << aExpected( i, j ) );
                CHECK( std::abs( cplx( aActual( i, j ) ) - cplx( aExpected( i, j ) ) ) <= aTolerance );
            }
        }
    }
}    // namespace moris::linalg_test
