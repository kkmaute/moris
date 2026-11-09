/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Lumped_Matrix.cpp
 *
 */

#include "catch.hpp"
#include "cl_FEM_Lumped_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_vectorize.hpp"

namespace moris::fem
{
    TEST_CASE( "IWG Lumped Matrix Operations", "[fem],[IWG_Lumped_Matrix]" )
    {
        // Two spatial nodes and two time bases give four independent diagonals.
        // Asymmetric time blocks detect swapped time indices. A second field
        // exercises the field-major coefficient ordering without cross coupling.
        for ( uint tNumFields : { 1u, 2u } )
        {
            Lumped_Matrix tLump;
            tLump.set_size( 2, 2, tNumFields );

            Matrix< DDRMat >       tReference( 4 * tNumFields, 4 * tNumFields, 0.0 );
            const Matrix< DDRMat > tFieldReference = {
                { 1.0, 0.0, 3.0, 0.0 },
                { 0.0, 2.0, 0.0, 4.0 },
                { 5.0, 0.0, 7.0, 0.0 },
                { 0.0, 6.0, 0.0, 8.0 }
            };

            for ( uint iField = 0; iField < tNumFields; ++iField )
            {
                const real tFieldScale = iField + 1.0;
                const uint tStart      = 4 * iField;
                tReference( { tStart, tStart + 3 }, { tStart, tStart + 3 } ) =
                        tFieldScale * tFieldReference;

                for ( uint iTime = 0; iTime < 2; ++iTime )
                {
                    for ( uint jTime = 0; jTime < 2; ++jTime )
                    {
                        tLump.diagonal( iField, iTime, jTime, 0 ) = tFieldScale * ( 1.0 + 4 * iTime + 2 * jTime );
                        tLump.diagonal( iField, iTime, jTime, 1 ) = tFieldScale * ( 2.0 + 4 * iTime + 2 * jTime );
                    }
                }
            }

            CHECK( tLump.num_stored_entries() == 8 * tNumFields );

            for ( uint tNumVectors : { 1u, 3u } )
            {
                Matrix< DDRMat > tValues( tReference.n_cols(), tNumVectors );

                for ( uint i = 0; i < tValues.n_rows(); ++i )
                {
                    for ( uint j = 0; j < tNumVectors; ++j )
                    {
                        tValues( i, j ) = 0.7 * i - 1.3 * j - 0.4;
                    }
                }

                for ( real tScale : { 0.0, -0.3, 2.5 } )
                {
                    INFO( "Fields=" << tNumFields << " RHS=" << tNumVectors << " scale=" << tScale );

                    const Matrix< DDRMat > tExpected = tScale * tReference * tValues;
                    const Matrix< DDRMat > tActual   = tScale * tLump * tValues;

                    CHECK( norm( tActual - tExpected ) < 1e-11 );

                    const Matrix< DDRMat > tChained = 0.5 * ( tLump * tScale ) * 2.0 * tValues;

                    CHECK( norm( tChained - tExpected ) < 1e-11 );
                }

                const Matrix< DDRMat > tUnscaled         = tLump * tValues;
                const Matrix< DDRMat > tUnscaledExpected = tReference * tValues;

                CHECK( norm( tUnscaled - tUnscaledExpected ) < 1e-11 );
            }

            // Products also accept backend expressions without a dense lumped
            // intermediate. Native time-continuity tests exercise vectorize().
            Matrix< DDRMat >       tOnes( tReference.n_cols(), 1, 1.0 );
            const Matrix< DDRMat > tExpressionProduct  = tLump * ( 1.7 * tOnes );
            const Matrix< DDRMat > tExpressionExpected = 1.7 * tReference * tOnes;

            CHECK( norm( tExpressionProduct - tExpressionExpected ) < 1e-11 );
#ifdef MORIS_USE_ARMA
            Matrix< DDRMat >       tCoefficients( 4, tNumFields, 1.0 );
            const Matrix< DDRMat > tVectorProduct = 1.7 * tLump * vectorize( tCoefficients );

            CHECK( norm( tVectorProduct - tExpressionExpected ) < 1e-11 );
#endif

            // Assembly must accumulate in a nonzero target, preserve entries
            // outside an offset block, and support both named and temporary views.
            const uint tSize = tLump.n_rows();

            for ( real tScale : { 0.0, -0.3, 2.5 } )
            {
                Matrix< DDRMat > tTarget( tSize + 4, tSize + 5, 0.125 );
                Matrix< DDRMat > tExpected = tTarget;
                tExpected( { 2, tSize + 1 }, { 3, tSize + 2 } ) += tScale * tReference;
                tTarget( { 2, tSize + 1 }, { 3, tSize + 2 } ) += tScale * tLump;

                CHECK( norm( tTarget - tExpected ) < 1e-11 );

                auto tBlock = tTarget( { 2, tSize + 1 }, { 3, tSize + 2 } );
                tBlock -= tScale * tLump;

                CHECK( norm( tTarget - Matrix< DDRMat >( tSize + 4, tSize + 5, 0.125 ) ) < 1e-11 );

                Matrix< DDRMat >       tWhole( tSize, tSize, 0.25 );
                const Matrix< DDRMat > tScaledReference = tScale * tReference;
                const Matrix< DDRMat > tWholeExpected   = tWhole + tScaledReference;
                tWhole += tScale * tLump;

                CHECK( norm( tWhole - tWholeExpected ) < 1e-11 );
                tWhole -= tScale * tLump;

                CHECK( norm( tWhole - Matrix< DDRMat >( tSize, tSize, 0.25 ) ) < 1e-11 );
            }

            Matrix< DDRMat > tUnscaledTarget( tSize, tSize, 0.0 );
            tUnscaledTarget += tLump;

            CHECK( norm( tUnscaledTarget - tReference ) < 1e-11 );
            tUnscaledTarget( { 0, tSize - 1 }, { 0, tSize - 1 } ) -= tLump;

            CHECK( norm( tUnscaledTarget ) < 1e-11 );
        }
    }

    TEST_CASE( "IWG Lumped Matrix Storage And Validation", "[fem],[IWG_Lumped_Matrix]" )
    {
        // Three-component HEX64 with two time bases stores 768 values instead
        // of a 384x384 matrix. Reusing the object for constant time resets storage.
        Lumped_Matrix tLump;
        tLump.set_size( 64, 2, 3 );

        CHECK( tLump.n_rows() == 384 );
        CHECK( tLump.n_cols() == 384 );
        CHECK( tLump.num_stored_entries() == 768 );
        tLump.set_size( 2, 1, 1 );
        tLump.diagonal( 0, 0, 0, 0 ) = 2.0;
        tLump.diagonal( 0, 0, 0, 1 ) = -3.0;

        Matrix< DDRMat >       tValues  = { { 4.0 }, { 5.0 } };
        const Matrix< DDRMat > tProduct = tLump * tValues;

        CHECK( tProduct( 0 ) == Approx( 8.0 ) );
        CHECK( tProduct( 1 ) == Approx( -15.0 ) );
        CHECK( tLump.num_stored_entries() == 2 );

        // Reject invalid sizes before touching the destination or coefficients.
        Matrix< DDRMat > tWrongSize( 3, 3, 1.0 );

        REQUIRE_THROWS( tLump.multiply( tWrongSize ) );
        REQUIRE_THROWS( tWrongSize += tLump );
        REQUIRE_THROWS( tWrongSize( { 0, 2 }, { 0, 2 } ) += 2.0 * tLump );
        REQUIRE_THROWS( tLump.set_size( 0, 1, 1 ) );
    }

}    // namespace moris::fem
