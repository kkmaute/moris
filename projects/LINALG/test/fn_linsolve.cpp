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
#include "typedefs.hpp"
#include "fn_linsolve.hpp"
#include "fn_norm.hpp"

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

