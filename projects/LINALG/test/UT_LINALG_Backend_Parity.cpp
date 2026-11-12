/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#include <catch.hpp>
#include "cl_Matrix.hpp"
#include "fn_chol_l.hpp"
#include "fn_diag_mat.hpp"
#include "fn_diag_vec.hpp"
#include "fn_eig_gen.hpp"
#include "fn_find.hpp"
#include "fn_histc.hpp"
#include "fn_interp1.hpp"
#include "fn_join_horiz.hpp"
#include "fn_lu.hpp"
#include "fn_max.hpp"
#include "fn_min.hpp"
#include "fn_norm.hpp"
#include "fn_prod.hpp"
#include "fn_qr.hpp"
#include "fn_rank.hpp"
#include "fn_sort.hpp"
#include "fn_sort_index.hpp"
#include "fn_sum.hpp"
#include "fn_svd.hpp"
#include "fn_trace.hpp"
#include "fn_trans.hpp"
#include "fn_vectorize.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
    TEST_CASE( "LINALG general eigenvalues include imaginary pairs", "[linalgebra],[backend_parity]" )
    {
        // A rotation has eigenvalues +i and -i; a symmetric solver is invalid.
        Matrix< DDCMat > tA = { { 0.0, -1.0 }, { 1.0, 0.0 } };
        Matrix< DDCMat > tValues, tVectors;
        eig_gen( tValues, tVectors, tA );
        REQUIRE( tValues.numel() == 2 );
        for ( uint j = 0; j < 2; ++j )
        {
            CHECK( std::abs( tValues( j ).real() ) < 1e-12 );
            CHECK( std::abs( tValues( j ).imag() ) == Approx( 1.0 ) );
            for ( uint i = 0; i < 2; ++i )
            {
                cplx tResidual = -tValues( j ) * tVectors( i, j );
                for ( uint k = 0; k < 2; ++k )
                {
                    tResidual += tA( i, k ) * tVectors( k, j );
                }
                CHECK( std::abs( tResidual ) < 1e-12 );
            }
        }
        CHECK( std::abs( tValues( 0 ) + tValues( 1 ) ) < 1e-12 );
    }

    TEST_CASE( "LINALG diagonals and histogram bins match documented shapes", "[linalgebra],[backend_parity]" )
    {
        Matrix< DDRMat > tVector   = { { 1.0 }, { 2.0 }, { 3.0 } };
        Matrix< DDRMat > tDiagonal = diag_mat( tVector, 1 );
        REQUIRE( tDiagonal.n_rows() == 4 );
        REQUIRE( tDiagonal.n_cols() == 4 );
        CHECK( tDiagonal( 0, 1 ) == 1.0 );
        CHECK( tDiagonal( 2, 3 ) == 3.0 );
        CHECK( sum( tDiagonal ) == 6.0 );
        Matrix< DDRMat > tA      = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }, { 7.0, 8.0, 9.0 } };
        Matrix< DDRMat > tOffset = diag_vec( tA, 1 );
        REQUIRE( tOffset.numel() == 2 );
        CHECK( tOffset( 0 ) == 2.0 );
        CHECK( tOffset( 1 ) == 6.0 );

        Matrix< DDUMat >  tSamples = { { 0 }, { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } };
        Matrix< DDUMat >  tEdges   = { { 0 }, { 3 }, { 5 } };
        Matrix< DDNIMat > tCounts  = histc( tSamples, tEdges );
        REQUIRE( tCounts.numel() == 3 );
        CHECK( tCounts( 0 ) == 3 );
        CHECK( tCounts( 1 ) == 2 );
        CHECK( tCounts( 2 ) == 1 );
    }

    TEST_CASE( "LINALG matrix reductions and spectral norm", "[linalgebra],[backend_parity]" )
    {
        const Matrix< DDRMat > tA       = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        Matrix< DDRMat >       tProduct = prod( tA );
        REQUIRE( tProduct.n_rows() == 1 );
        REQUIRE( tProduct.n_cols() == 2 );
        CHECK( tProduct( 0 ) == 3.0 );
        CHECK( tProduct( 1 ) == 8.0 );
        Matrix< DDRMat > tMinimum = min( tA );
        Matrix< DDRMat > tMaximum = max( tA );
        CHECK( tMinimum( 0 ) == 1.0 );
        CHECK( tMinimum( 1 ) == 2.0 );
        CHECK( tMaximum( 0 ) == 3.0 );
        CHECK( tMaximum( 1 ) == 4.0 );
        CHECK( trace( tA ) == 5.0 );
        const Matrix< DDRMat > tD = { { 3.0, 0.0 }, { 0.0, 4.0 } };
        CHECK( norm( tD ) == Approx( 4.0 ) );
        Matrix< DDNIMat > tIndices = find( tA, 0 );
        CHECK( tIndices.numel() == 4 );
    }

    TEST_CASE( "LINALG rectangular factorizations reconstruct the input", "[linalgebra],[backend_parity]" )
    {
        const Matrix< DDRMat > tA = { { 1.0, 2.0 }, { 3.0, 4.0 }, { 5.0, 7.0 } };
        Matrix< DDRMat >       tU, tS, tV;
        svd( tU, tS, tV, tA );
        REQUIRE( tU.n_rows() == 3 );
        REQUIRE( tU.n_cols() == 3 );
        REQUIRE( tV.n_rows() == 2 );
        REQUIRE( tV.n_cols() == 2 );
        Matrix< DDRMat > tSigma( 3, 2, 0.0 );
        tSigma( 0, 0 ) = tS( 0 );
        tSigma( 1, 1 ) = tS( 1 );
        CHECK( norm( tU * tSigma * trans( tV ) - tA ) < 1e-12 );
        Matrix< DDRMat > tQ, tR, tL, tP;
        qr( tQ, tR, tA );
        CHECK( norm( tQ * tR - tA ) < 1e-12 );
        lu( tL, tU, tP, tA );
        CHECK( norm( tP * tA - tL * tU ) < 1e-12 );
        CHECK( rank( tA ) == 2 );
        Matrix< DDRMat > tInvalid = { { 1.0, 0.0 }, { 0.0, -1.0 } };
        CHECK_THROWS( Matrix< DDRMat >( chol_l( tInvalid ) ) );
    }

    TEST_CASE( "LINALG sorting interpolation and concatenation", "[linalgebra],[backend_parity]" )
    {
        const Matrix< DDRMat > tA = { { 3.0 }, { 1.0 }, { 2.0 } };
        Matrix< DDRMat >       tSorted;
        sort( tA, tSorted );
        REQUIRE( tSorted.numel() == 3 );
        CHECK( tSorted( 0 ) == 1.0 );
        CHECK( tSorted( 2 ) == 3.0 );
        Matrix< DDUMat > tOrder = sort_index( tA, "descend" );
        CHECK( tOrder( 0 ) == 0 );
        CHECK( tOrder( 1 ) == 2 );
        CHECK( tOrder( 2 ) == 1 );
        Matrix< DDRMat > tJoined = join_horiz( tA, tSorted );
        REQUIRE( tJoined.n_rows() == 3 );
        REQUIRE( tJoined.n_cols() == 2 );
        CHECK( tJoined( 0, 1 ) == 1.0 );
        const Matrix< DDRMat > tX     = { { 0.0 }, { 1.0 }, { 2.0 } };
        const Matrix< DDRMat > tY     = { { 0.0 }, { 2.0 }, { 4.0 } };
        const Matrix< DDRMat > tQuery = { { -1.0 }, { 0.5 }, { 1.5 }, { 3.0 } };
        Matrix< DDRMat >       tInterpolated;
        interp1( tX, tY, tQuery, tInterpolated );
        REQUIRE( tInterpolated.numel() == 4 );
        CHECK( std::isnan( tInterpolated( 0 ) ) );
        CHECK( tInterpolated( 1 ) == 1.0 );
        CHECK( tInterpolated( 2 ) == 3.0 );
        CHECK( std::isnan( tInterpolated( 3 ) ) );
    }

    TEST_CASE( "LINALG borrowed storage and column views retain bindings", "[linalgebra],[backend_parity]" )
    {
        double           tStorage[] = { 1.0, 2.0, 3.0, 4.0 };
        Matrix< DDRMat > tBorrowed( tStorage, 2, 2, false, true );
        auto             tColumn = tBorrowed.get_column( 1 );
        tColumn( 0 )             = 9.0;
        CHECK( tStorage[ 2 ] == 9.0 );
        const Matrix< DDRMat > tOther( 2, 2, 1.0 );
        tBorrowed = tBorrowed + tOther;
        CHECK( tBorrowed.data() == tStorage );
        CHECK( tStorage[ 2 ] == 10.0 );
        CHECK_THROWS( tBorrowed.set_size( 3, 2 ) );
        Matrix< DDRMat > tCopy = tBorrowed;
        tCopy( 0, 0 )          = -5.0;
        CHECK( tStorage[ 0 ] == 2.0 );
        Matrix< DDRMat > tLoose( tStorage, 2, 2, false, false );
        tLoose.set_size( 3, 2 );
        CHECK( tLoose.data() != tStorage );
        Matrix< DDRMat > tRow = { { 1.0 }, { 2.0 } };
        tBorrowed.set_row( 0, tRow );
        CHECK( tStorage[ 0 ] == 1.0 );
        CHECK( tStorage[ 2 ] == 2.0 );
    }
    TEST_CASE( "LINALG fixed matrices share the arithmetic and view interface", "[linalgebra],[backend_parity]" )
    {
        Matrix< F33RMat > tA( 3, 3, 0.0 );
        tA( 0, 0 )                 = 1.0;
        tA( 1, 1 )                 = 2.0;
        tA( 2, 2 )                 = 3.0;
        const Matrix< F33RMat > tB = tA + tA;
        CHECK( tB( 2, 2 ) == 6.0 );
        Matrix< F31RMat > tColumn = tA.get_column( 1 );
        CHECK( tColumn( 1 ) == 2.0 );
        auto tView = tA.get_column( 1 );
        tView( 1 ) = 7.0;
        CHECK( tA( 1, 1 ) == 7.0 );
        Matrix< F33RMat > tProduct = tA * tA;
        CHECK( tProduct( 1, 1 ) == 49.0 );
        CHECK( tA.capacity() == 9 * sizeof( real ) );
    }
}    // namespace moris
