/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_linsolve.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_linsolve.hpp"
#include "fn_norm.hpp"
#include "UT_LINALG_Test_Helpers.hpp"

namespace moris
{
    TEST_CASE( "moris::linsolve", "[linalgebra],[linsolve]" )
    {
        SECTION("linsolve_real" )
        {
            // matrix size
            uint imat = 20;

            for (uint iz=0;iz<3;++iz)
            {
                Matrix< DDRMat > tA( imat, imat, 0.0 );
                Matrix< DDRMat > tB( imat, 1 , 0.0 );

                // build matrix and RHS
                tA(0,0) = 2.0*imat;
                tA(0,1) =-1.0*imat;

                for (uint ir=1;ir<imat-1;++ir)
                {
                    tA(ir,ir-1)=-1.0*imat;
                    tA(ir,ir  )= 2.0*imat;
                    tA(ir,ir+1)=-1.0*imat;

                    tB(ir,0) = (real)ir/imat;
                }

                tA(imat-1,imat-2) =-1.0*imat;
                tA(imat-1,imat-1) = 2.0*imat;

                // solve system
                Matrix< DDRMat > tSol = moris::solve( tA, tB );

                // check for accuracy
                real tError = norm( tB - tA*tSol);

                REQUIRE( tError < 1e-8 * norm(tB) );

                imat += imat;
            }
       }

        //    SECTION( "linsolve_complex" )
        //    {
        //        moris::Matrix< moris::DDCMat > A( 2, 2 );
        //        moris::Matrix< moris::DDCMat > B( 2, 1 );
        //
        //        A( 0, 0 ) = { 1.0, 1.0 }; A( 0, 1 ) = { 2.0, -1.0 };
        //        A( 1, 0 ) = { 7.0, 0.0 }; A( 1, 1 ) = { 8.0, -2.0 };
        //
        //        B( 0, 0 ) = { 2.0,  7.0 };
        //        B( 1, 0 ) = { 4.0, -9.0 };
        //
        //        moris::Matrix< moris::DDCMat > C = moris::solve( A,B );
        //
        //        REQUIRE( moris::equal_to( C( 0,0 ).real(),  4.529729729729730 ) );
        //        REQUIRE( moris::equal_to( C( 0,0 ).imag(), -3.778378378378378 ) );
        //        REQUIRE( moris::equal_to( C( 1,0 ).real(), -3.772972972972973 ) );
        //        REQUIRE( moris::equal_to( C( 1,0 ).imag(),  1.237837837837838 ) );
        //    }
            }
}

namespace moris
{
    TEST_CASE( "LINALG square solves preserve scaled inputs expressions and aliases", "[linalgebra],[linsolve]" )
    {
        // An indefinite, nonsymmetric system exercises the general square path.
        // Two independently known solutions check multiple RHS columns.
        Matrix< DDRMat > tA( 12, 12, 0.0 ), tExpected( 12, 2 );
        for ( uint i = 0; i < 12; ++i )
        {
            tA( i, i ) = i % 2 == 0 ? 4.0 : -4.0;
            if ( i + 1 < 12 )
            {
                tA( i, i + 1 ) = 0.5;
                tA( i + 1, i ) = -0.25;
            }
            tExpected( i, 0 ) = 1.0 + i / 12.0;
            tExpected( i, 1 ) = -2.0 + i / 12.0;
        }
        const Matrix< DDRMat > tB = tA * tExpected;

        for ( real tScale : { 1e-150, 1.0, 1e150 } )
        {
            CAPTURE( tScale );
            linalg_test::check_matrix( Matrix< DDRMat >( moris::solve( tScale * tA, tScale * tB ) ), tExpected );
        }

        Matrix< DDRMat > tAliasedA = tA;
        tAliasedA                  = solve( tAliasedA, tB );
        linalg_test::check_matrix( tAliasedA, tExpected );

        Matrix< DDRMat > tAliasedB = tB;
        tAliasedB                  = solve( tA, tAliasedB );
        linalg_test::check_matrix( tAliasedB, tExpected );

        // Widely different RHS magnitudes must each retain relative accuracy.
        Matrix< DDRMat > tScaledB = tB;
        for ( uint i = 0; i < 12; ++i )
        {
            tScaledB( i, 0 ) *= 1e150;
            tScaledB( i, 1 ) *= 1e-150;
        }
        const Matrix< DDRMat > tScaledX = solve( tA, tScaledB );
        for ( uint i = 0; i < 12; ++i )
        {
            CHECK( tScaledX( i, 0 ) / 1e150 == Approx( tExpected( i, 0 ) ).margin( 1e-12 ) );
            CHECK( tScaledX( i, 1 ) / 1e-150 == Approx( tExpected( i, 1 ) ).margin( 1e-12 ) );
        }
    }

    TEST_CASE( "LINALG difficult solves retain least squares and minimum norm solutions", "[linalgebra],[linsolve]" )
    {
        SECTION( "Ill conditioned full rank system" )
        {
            Matrix< DDRMat > tA( 12, 12, 0.0 ), tB( 12, 1, 1.0 );
            for ( uint i = 0; i < 12; ++i )
            {
                tA( i, i ) = 1.0;
            }
            tA( 11, 11 ) = 1e-10;
            tB( 11 )     = 1e-10;
            linalg_test::check_matrix( Matrix< DDRMat >( solve( tA, tB ) ), Matrix< DDRMat >( 12, 1, 1.0 ) );
        }

        SECTION( "Singular inconsistent square system" )
        {
            // The last block has identical columns. Least squares gives their
            // sum as 2.6, and minimum norm splits that sum equally.
            Matrix< DDRMat > tA( 12, 12, 0.0 ), tB( 12, 2, 1.0 ), tExpected( 12, 2, 1.0 );
            for ( uint i = 0; i < 10; ++i )
            {
                tA( i, i ) = 1.0;
            }
            tA( 10, 10 ) = tA( 10, 11 ) = 1.0;
            tA( 11, 10 ) = tA( 11, 11 ) = 2.0;
            tB( 10, 0 )                 = 3.0;
            tB( 11, 0 )                 = 5.0;
            tB( 10, 1 )                 = -2.0;
            tB( 11, 1 )                 = 1.0;
            tExpected( 10, 0 ) = tExpected( 11, 0 ) = 1.3;
            tExpected( 10, 1 ) = tExpected( 11, 1 ) = 0.0;
            linalg_test::check_matrix( Matrix< DDRMat >( solve( tA, tB ) ), tExpected );
        }

        SECTION( "Rectangular inconsistent system" )
        {
            Matrix< DDRMat > tA( 12, 2, 0.0 ), tB( 12, 1, 0.0 );
            tA( 0, 0 ) = tA( 2, 0 ) = 1.0;
            tA( 1, 1 ) = tA( 2, 1 ) = 1.0;
            tB( 0 )                 = 1.0;
            tB( 1 )                 = 2.0;
            tB( 2 )                 = 4.0;
            linalg_test::check_matrix( Matrix< DDRMat >( solve( tA, tB ) ), Matrix< DDRMat >{ { 4.0 / 3.0 }, { 7.0 / 3.0 } } );
        }
    }

    TEST_CASE( "LINALG complex square solves preserve multiple right hand sides", "[linalgebra],[linsolve]" )
    {
        // A non-Hermitian system prevents an implicit SPD assumption.
        Matrix< DDCMat > tA( 12, 12, cplx( 0.0 ) ), tExpected( 12, 2 );
        for ( uint i = 0; i < 12; ++i )
        {
            tA( i, i ) = cplx( 4.0, 1.0 );
            if ( i + 1 < 12 )
            {
                tA( i, i + 1 ) = cplx( 0.5, -0.25 );
                tA( i + 1, i ) = cplx( -0.1, 0.3 );
            }
            tExpected( i, 0 ) = cplx( 1.0 + i, -0.5 );
            tExpected( i, 1 ) = cplx( -2.0, i / 12.0 );
        }
        const Matrix< DDCMat > tB = tA * tExpected;
        linalg_test::check_matrix( Matrix< DDCMat >( solve( tA, tB ) ), tExpected );
    }

#ifdef MORIS_USE_EIGEN
    TEST_CASE( "LINALG Eigen square solves retain robust fallbacks", "[linalgebra],[linsolve]" )
    {
        SECTION( "Nearly singular rank decision" )
        {
            // A tiny residual alone cannot detect this unstable direction.
            // Preserve COD's rank truncation rather than accepting LU's value.
            Matrix< DDRMat > tA( 12, 12, 0.0 ), tB( 12, 1, 1.0 ), tExpected( 12, 1, 1.0 );
            for ( uint i = 0; i < 12; ++i )
            {
                tA( i, i ) = 1.0;
            }
            tA( 11, 11 )    = 1e-18;
            tB( 11 )        = 1e-18;
            tExpected( 11 ) = 0.0;
            linalg_test::check_matrix( Matrix< DDRMat >( solve( tA, tB ) ), tExpected );
        }

        SECTION( "Wilkinson pivot growth" )
        {
            // This well-conditioned matrix has exponential LU pivot growth.
            // The residual guard must catch lost accuracy even when rcond is
            // acceptable; the independently known solution remains attainable.
            Matrix< DDRMat > tA( 48, 48, 0.0 ), tExpected( 48, 1 );
            for ( uint i = 0; i < 48; ++i )
            {
                tA( i, i ) = 1.0;
                for ( uint j = 0; j < i; ++j )
                {
                    tA( i, j ) = -1.0;
                }
                tA( i, 47 )    = 1.0;
                tExpected( i ) = 1.0 / ( i + 1.0 );
            }
            const Matrix< DDRMat > tB = tA * tExpected;
            linalg_test::check_matrix( Matrix< DDRMat >( solve( tA, tB ) ), tExpected );
        }
    }
#endif
}    // namespace moris
