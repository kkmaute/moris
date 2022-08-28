/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sylvester.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_sylvester.hpp"
#include "fn_norm.hpp"

namespace moris
{
    TEST_CASE( "moris::sylvester", "[linalgebra],[sylvester]" )
    {
        SECTION("sylvester_real_symmetric" )
        {
            // matrix size
            uint imat = 10;

            Matrix< DDRMat > tA( imat, imat, 0.0 );
            Matrix< DDRMat > tB( imat, imat, 0.0 );
            Matrix< DDRMat > tC( imat, imat, 0.0 );

            // build symmetric matrix
            tA(0,0) = 2.0*imat;
            tB(0,0) = 2.0*imat/3.0;
            tC(0,0) = 2.0*imat/2.0;

            for (uint ir=1;ir<imat-1;++ir)
            {
                tA(ir,ir-1)=-1.0*imat;
                tA(ir,ir  )= 2.0*imat;
                tA(ir,ir+1)=-1.0*imat;

                tB(ir,ir-1)=-1.0*imat/6.0;
                tB(ir,ir  )= 2.0*imat/3.0;
                tB(ir,ir+1)=-1.0*imat/6.0;

                tC(ir,ir-1)=-1.0*imat/4.0;
                tC(ir,ir  )= 2.0*imat/2.0;
                tC(ir,ir+1)=-1.0*imat/4.0;
            }

            tA(imat-1,imat-2) =-1.0*imat;
            tA(imat-1,imat-1) = 2.0*imat;

            tB(imat-1,imat-2) =-1.0*imat/6.0;
            tB(imat-1,imat-1) = 2.0*imat/3.0;

            tC(imat-1,imat-2) =-1.0*imat/4.0;
            tC(imat-1,imat-1) = 2.0*imat/2.0;

            // compute solution to tA*tX + tX*tB + tC = 0 - Option 1
             Matrix< DDRMat > tX = moris::sylvester(tA,tB,tC);

             // check for accuracy
             real tError = norm( tA*tX + tX*tB + tC);

             REQUIRE( tError < 1e-8 * norm(tC) );

             // compute solution to tA*tX + tX*tB + tC = 0 - Option 2
             tX.fill(0.0);
             moris::sylvester(tA,tB,tC,tX);

             // check for accuracy
             tError = norm( tA*tX + tX*tB + tC);

             REQUIRE( tError < 1e-8 * norm(tC) );
        }

        SECTION("sylvester_real_nonsymmetric" )
        {
            // matrix size
            uint imat = 20;

            Matrix< DDRMat > tA( imat, imat, 0.0 );
            Matrix< DDRMat > tB( imat, imat, 0.0 );
            Matrix< DDRMat > tC( imat, imat, 0.0 );

            // build symmetric matrix
            tA(0,0) = 2.0*imat;
            tB(0,0) = 2.0*imat/3.0;
            tC(0,0) = 2.0*imat/2.0;

            for (uint ir=1;ir<imat-1;++ir)
            {
                tA(ir,ir-1)=-1.0*imat;
                tA(ir,ir  )= 2.0*imat;
                tA(ir,ir+1)=-1.0*imat;

                tB(ir,ir-1)=-0.9*imat/6.0;
                tB(ir,ir  )= 2.0*imat/3.0;
                tB(ir,ir+1)=-1.1*imat/6.0;

                tC(ir,ir-1)=-0.8*imat/4.0;
                tC(ir,ir  )= 2.0*imat/2.0;
                tC(ir,ir+1)=-1.3*imat/4.0;
            }

            tA(imat-1,imat-2) =-1.0*imat;
            tA(imat-1,imat-1) = 2.0*imat;

            tB(imat-1,imat-2) =-0.9*imat/6.0;
            tB(imat-1,imat-1) = 2.0*imat/3.0;

            tC(imat-1,imat-2) =-0.8*imat/4.0;
            tC(imat-1,imat-1) = 2.0*imat/2.0;

            // compute solution to tA*tX + tX*tB + tC = 0 - Option 1
            Matrix< DDRMat > tX = moris::sylvester(tA,tB,tC);

            // check for accuracy
            real tError = norm( tA*tX + tX*tB + tC);

            REQUIRE( tError < 1e-8 * norm(tC) );

            // compute solution to tA*tX + tX*tB + tC = 0 - Option 2
            tX.fill(0.0);
            moris::sylvester(tA,tB,tC,tX);

            // check for accuracy
            tError = norm( tA*tX + tX*tB + tC);

            REQUIRE( tError < 1e-8 * norm(tC) );
        }
    }
}

