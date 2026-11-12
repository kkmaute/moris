/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#include "UT_LINALG_Test_Helpers.hpp"
#include "fn_all_true.hpp"
#include "fn_cond.hpp"
#include "fn_diag_mat.hpp"
#include "fn_diag_vec.hpp"
#include "fn_eye.hpp"
#include "fn_find.hpp"
#include "fn_find_unique.hpp"
#include "fn_interp1.hpp"
#include "fn_intersect.hpp"
#include "fn_join_horiz.hpp"
#include "fn_linsolve.hpp"
#include "fn_linspace.hpp"
#include "fn_lu.hpp"
#include "fn_max.hpp"
#include "fn_min.hpp"
#include "fn_norm.hpp"
#include "fn_prod.hpp"
#include "fn_qr.hpp"
#include "fn_rank.hpp"
#include "fn_reindex_mat.hpp"
#include "fn_sort_index.hpp"
#include "fn_sum.hpp"
#include "fn_trace.hpp"
#include "fn_trans.hpp"
#include <limits>
#include <vector>

namespace moris
{
    using linalg_test::check_matrix;

    TEST_CASE( "LINALG native vector extrema return scalar values", "[linalgebra],[function_edges]" )
    {
        Matrix< DDRMat > tA       = { { 1.0, -2.0 }, { 3.0, 4.0 }, { -5.0, 6.0 } };
        const real       tMinimum = min( tA.get_column( 0 ) );
        const real       tMaximum = max( tA.get_row( 1 ) );
        CHECK( tMinimum == -5.0 );
        CHECK( tMaximum == 4.0 );
        const real tExpressionMinimum = min( tA.get_column( 0 ) + tA.get_column( 0 ) );
        CHECK( tExpressionMinimum == -10.0 );

        Matrix< DDCMat > tComplex        = { { cplx( 1.0, 1.0 ), cplx( 4.0, -2.0 ) },
                   { cplx( 2.0, 1.0 ), cplx( 0.5, 0.0 ) } };
        const cplx       tComplexMinimum = min( tComplex.get_column( 1 ) );
        const cplx       tComplexMaximum = max( tComplex.get_row( 0 ) );
        CHECK( tComplexMinimum == cplx( 0.5, 0.0 ) );
        CHECK( tComplexMaximum == cplx( 4.0, -2.0 ) );
    }

    TEST_CASE( "LINALG reindexing selects a row and preserves repeated indices", "[linalgebra],[function_edges]" )
    {
        const Matrix< DDRMat > tSource  = { { 11.0, 12.0, 13.0, 14.0 }, { -1.0, 0.0, 7.0, 9.0 } };
        const Matrix< DDUMat > tIndices = { { 3, 1, 1 }, { 0, 2, 3 } };
        check_matrix( ::reindex_mat( tIndices, 1, tSource ), Matrix< DDRMat >{ { 9.0, 0.0, 0.0 }, { -1.0, 7.0, 9.0 } } );
        CHECK( tSource( 1, 3 ) == 9.0 );
        const Matrix< DDUMat > tEmpty;
        CHECK( ::reindex_mat( tEmpty, 0, tSource ).numel() == 0 );
    }

    TEST_CASE( "LINALG rank and condition distinguish singular and rectangular systems", "[linalgebra],[function_edges]" )
    {
        const Matrix< DDRMat > tEmpty;
        CHECK( rank( tEmpty ) == 0 );
        CHECK( rank( Matrix< DDRMat >( 3, 2, 0.0 ) ) == 0 );
        const Matrix< DDRMat > tDependent = { { 1.0, 2.0, 0.0 }, { 2.0, 4.0, 0.0 }, { 0.0, 0.0, 3.0 } };
        CHECK( rank( tDependent ) == 2 );
        const Matrix< DDRMat > tScaled = 1e-8 * tDependent;
        CHECK( rank( tScaled ) == 2 );
        const Matrix< DDRMat > tRectangular = { { 3.0, 0.0 }, { 0.0, 1.0 }, { 0.0, 0.0 } };
        CHECK( cond( tRectangular ) == Approx( 3.0 ) );
        const Matrix< DDRMat > tSingular = { { 3.0, 0.0 }, { 0.0, 0.0 } };
        CHECK( std::isinf( cond( tSingular ) ) );
    }

    TEST_CASE( "LINALG trace and extrema include rectangular and complex values", "[linalgebra],[function_edges]" )
    {
        const Matrix< DDRMat > tA = { { 1.0, -2.0, 8.0 }, { 3.0, 4.0, -5.0 } };
        const Matrix< DDRMat > tB = { { 2.0, -3.0, 7.0 }, { -1.0, 6.0, -4.0 } };
        CHECK( trace( tA ) == 5.0 );
        check_matrix( Matrix< DDRMat >( min( tA, tB ) ), Matrix< DDRMat >{ { 1.0, -3.0, 7.0 }, { -1.0, 4.0, -5.0 } } );
        check_matrix( Matrix< DDRMat >( max( tA, tB ) ), Matrix< DDRMat >{ { 2.0, -2.0, 8.0 }, { 3.0, 6.0, -4.0 } } );
        const Matrix< DDCMat > tComplex = { { cplx( 1, 2 ), cplx( 3 ) }, { cplx( 4 ), cplx( 5, -1 ) } };
        CHECK( trace( tComplex ) == cplx( 6, 1 ) );
        CHECK( trace( Matrix< DDRMat >() ) == 0.0 );
    }

    TEST_CASE( "LINALG reductions preserve vector orientation and empty shapes", "[linalgebra],[function_edges]" )
    {
        const Matrix< DDRMat > tRow    = { { 2.0, 3.0, 4.0 } };
        const Matrix< DDRMat > tColumn = { { 2.0 }, { 3.0 }, { 4.0 } };
        check_matrix( Matrix< DDRMat >( prod( tRow ) ), tRow );
        check_matrix( Matrix< DDRMat >( prod( tColumn ) ), Matrix< DDRMat >( 1, 1, 24.0 ) );
        const Matrix< DDRMat > tEmpty;
        CHECK( sum( tEmpty ) == 0.0 );
        CHECK( Matrix< DDRMat >( prod( tEmpty ) ).numel() == 0 );
        check_matrix( Matrix< DDRMat >( sum( tEmpty, 0 ) ), Matrix< DDRMat >( 1, 0 ) );
        check_matrix( Matrix< DDRMat >( sum( tEmpty, 1 ) ), Matrix< DDRMat >( 0, 1 ) );
        CHECK( all_true( Matrix< DDBMat >() ) );
        CHECK_FALSE( all_true( Matrix< DDBMat >{ { true, true }, { true, false } } ) );
    }

    TEST_CASE( "LINALG concatenation retains widths and handles destination aliases", "[linalgebra],[function_edges]" )
    {
        Matrix< DDRMat >       tA = { { 1.0 }, { 2.0 } };
        const Matrix< DDRMat > tB = { { 3.0, 4.0 }, { 5.0, 6.0 } };
        tA                        = join_horiz( tA, tB );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 3.0, 4.0 }, { 2.0, 5.0, 6.0 } } );
        const Matrix< DDRMat > tDefaultEmpty;
        check_matrix( Matrix< DDRMat >( join_horiz( tDefaultEmpty, tB ) ), tB );
        check_matrix( Matrix< DDRMat >( join_horiz( tB, tDefaultEmpty ) ), tB );
        check_matrix( Matrix< DDRMat >( join_horiz( tDefaultEmpty, tDefaultEmpty ) ), tDefaultEmpty );
        const Matrix< DDRMat > tEmpty( 2, 0 );
        check_matrix( Matrix< DDRMat >( join_horiz( tEmpty, tB ) ), tB );
        check_matrix( Matrix< DDRMat >( join_horiz( tB, tEmpty ) ), tB );
    }

    TEST_CASE( "LINALG sorting indices form ordered permutations with duplicates", "[linalgebra],[function_edges]" )
    {
        const Matrix< DDRMat > tA          = { { 3.0 }, { -2.0 }, { 3.0 }, { 1.0 } };
        const Matrix< DDUMat > tAscending  = sort_index( tA );
        const Matrix< DDUMat > tDescending = sort_index( tA, "descend" );
        for ( bool tDescendingOrder : { false, true } )
        {
            const auto& tOrder = tDescendingOrder ? tDescending : tAscending;
            REQUIRE( tOrder.n_rows() == 4 );
            REQUIRE( tOrder.n_cols() == 1 );
            std::vector< bool > tSeen( 4, false );
            for ( uint i = 0; i < 4; ++i )
            {
                REQUIRE( tOrder( i ) < 4 );
                CHECK_FALSE( tSeen[ tOrder( i ) ] );
                tSeen[ tOrder( i ) ] = true;
                if ( i > 0 )
                {
                    CHECK( ( tDescendingOrder ? tA( tOrder( i - 1 ) ) >= tA( tOrder( i ) ) : tA( tOrder( i - 1 ) ) <= tA( tOrder( i ) ) ) );
                }
            }
        }
        CHECK( sort_index( Matrix< DDRMat >( 0, 1 ) ).numel() == 0 );
    }

    TEST_CASE( "LINALG find uses exact nonzero values and first occurrence indices", "[linalgebra],[function_edges]" )
    {
        const Matrix< DDRMat > tA = { { 0.0, 1e-30, 0.0 }, { 2.0, 0.0, -3.0 } };
        check_matrix( Matrix< DDNIMat >( find( tA ) ), Matrix< DDNIMat >{ { 1 }, { 2 }, { 5 } } );
        check_matrix( Matrix< DDNIMat >( find( tA, 2 ) ), Matrix< DDNIMat >{ { 1 }, { 2 } } );
        CHECK( Matrix< DDNIMat >( find( tA, 20 ) ).numel() == 3 );
        CHECK( Matrix< DDNIMat >( find( Matrix< DDRMat >( 2, 3, 0.0 ) ) ).numel() == 0 );
        Matrix< DDRMat > tDuplicates = { { 3.0, 1.0, 3.0, 2.0, 1.0 } };
        check_matrix( Matrix< DDNIMat >( find_unique( tDuplicates ) ), Matrix< DDNIMat >{ { 0 }, { 1 }, { 3 } } );
    }

    TEST_CASE( "LINALG intersection indices identify each common value", "[linalgebra],[function_edges]" )
    {
        const Matrix< DDRMat > tA = { { 4.0, 2.0, 2.0, 1.0 } };
        const Matrix< DDRMat > tB = { { 2.0, 4.0, 4.0, 3.0 } };
        Matrix< DDRMat >       tCommon;
        Matrix< DDUMat >       tFirst, tSecond;
        intersect( tA, tB, tCommon, tFirst, tSecond );
        check_matrix( tCommon, Matrix< DDRMat >{ { 2.0, 4.0 } } );
        check_matrix( Matrix< DDRMat >( intersect( tA, tB ) ), tCommon );
        REQUIRE( tFirst.numel() == 2 );
        REQUIRE( tSecond.numel() == 2 );
        for ( uint i = 0; i < 2; ++i )
        {
            REQUIRE( tFirst( i ) < tA.numel() );
            REQUIRE( tSecond( i ) < tB.numel() );
            CHECK( tA( tFirst( i ) ) == tCommon( i ) );
            CHECK( tB( tSecond( i ) ) == tCommon( i ) );
        }
    }

    TEST_CASE( "LINALG interpolation sorts sample locations and preserves query shape", "[linalgebra],[function_edges]" )
    {
        const Matrix< DDRMat > tX     = { { 2.0 }, { 0.0 }, { 1.0 } };
        const Matrix< DDRMat > tY     = { { 4.0 }, { 0.0 }, { 2.0 } };
        Matrix< DDRMat >       tQuery = { { 0.0, 0.5, 2.0 } };
        interp1( tX, tY, tQuery, tQuery );
        check_matrix( tQuery, Matrix< DDRMat >{ { 0.0, 1.0, 4.0 } } );
    }

    TEST_CASE( "LINALG identity and linspace include rectangular and endpoint cases", "[linalgebra],[function_edges]" )
    {
        Matrix< DDRMat > tA( 1, 1, 99.0 );
        eye( 2, 4, tA );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0, 0.0 } } );
        check_matrix( Matrix< DDRMat >( eye( 3, 2 ) ), Matrix< DDRMat >{ { 1.0, 0.0 }, { 0.0, 1.0 }, { 0.0, 0.0 } } );
        check_matrix( Matrix< DDRMat >( linspace( 3.0, -1.0, 3 ) ), Matrix< DDRMat >{ { 3.0 }, { 1.0 }, { -1.0 } } );
        check_matrix( Matrix< DDRMat >( linspace( 1.0, 3.0, 1 ) ), Matrix< DDRMat >( 1, 1, 3.0 ) );
        CHECK( Matrix< DDRMat >( linspace( 1.0, 3.0, 0 ) ).numel() == 0 );
    }

    TEST_CASE( "LINALG diagonal offsets select subdiagonals as documented", "[linalgebra],[function_edges]" )
    {
        Matrix< DDRMat >        tA           = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }, { 7.0, 8.0, 9.0 } };
        const Matrix< DDRMat >& tConst       = tA;
        const size_t            tSubdiagonal = static_cast< size_t >( -1 );
        check_matrix( Matrix< DDRMat >( diag_vec( tConst, tSubdiagonal ) ), Matrix< DDRMat >{ { 4.0 }, { 8.0 } } );
        check_matrix( Matrix< DDRMat >( diag_mat( tA, tSubdiagonal ) ), Matrix< DDRMat >{ { 0.0, 0.0, 0.0 }, { 4.0, 0.0, 0.0 }, { 0.0, 8.0, 0.0 } } );
    }

    TEST_CASE( "LINALG rectangular solves and wide factorizations satisfy their equations", "[linalgebra],[function_edges]" )
    {
        Matrix< DDRMat > tA( 12, 2 ), tB( 12, 1 );
        for ( uint i = 0; i < 12; ++i )
        {
            tA( i, 0 ) = 1.0;
            tA( i, 1 ) = i;
            tB( i )    = 3.0 + 2.0 * i;
        }
        check_matrix( Matrix< DDRMat >( solve( tA, tB ) ), Matrix< DDRMat >{ { 3.0 }, { 2.0 } } );
        Matrix< DDRMat > tWide( 12, 14, 0.0 );
        Matrix< DDRMat > tExpected( 14, 1, 0.0 );
        for ( uint i = 0; i < 12; ++i )
        {
            tWide( i, i )  = 1.0;
            tExpected( i ) = tB( i );
        }
        check_matrix( Matrix< DDRMat >( solve( tWide, tB ) ), tExpected );

        const Matrix< DDRMat > tSmall = { { 1.0, 2.0, 3.0 }, { 4.0, 6.0, 8.0 } };
        Matrix< DDRMat >       tQ, tR, tL, tU, tP;
        qr( tQ, tR, tSmall );
        CHECK( norm( tQ * tR - tSmall ) < 1e-12 );
        const Matrix< DDRMat > tI = eye( 2, 2 );
        CHECK( norm( trans( tQ ) * tQ - tI ) < 1e-12 );
        lu( tL, tU, tP, tSmall );
        CHECK( norm( tP * tSmall - tL * tU ) < 1e-12 );
    }
}    // namespace moris
