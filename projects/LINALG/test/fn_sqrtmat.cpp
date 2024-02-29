/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sqrtmat.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_sqrtmat.hpp"
#include "fn_norm.hpp"

namespace moris
{
    TEST_CASE( "moris::sqrtmat", "[linalgebra],[sqrtmat]" )
    {
        SECTION("sqrtmat_real_symmetric" )
        {
            // matrix size
            uint imat = 20;

            Matrix< DDRMat > tA( imat, imat, 0.0 );

            // build symmetric matrix
            tA(0,0) = 2.0*imat;

            for (uint ir=1;ir<imat-1;++ir)
            {
                tA(ir,ir-1)=-1.0*imat;
                tA(ir,ir  )= 2.0*imat;
                tA(ir,ir+1)=-1.0*imat;
            }

            tA(imat-1,imat-2) =-1.0*imat;
            tA(imat-1,imat-1) = 2.0*imat;

            // compute square root - Option 1
            Matrix< DDRMat > tSqrt = moris::sqrtmat(tA);

            // check for accuracy
            real tError = norm( tA - tSqrt*tSqrt);

            REQUIRE( tError < 1e-8 * norm(tA) );

            // compute square root - Option 2
            tSqrt.fill(0.0);
            moris::sqrtmat(tA, tSqrt);

            // check for accuracy
            tError = norm( tA - tSqrt*tSqrt);

            REQUIRE( tError < 1e-8 * norm(tA) );
        }

        SECTION("sqrtmat_real_symmetric" )
        {
            // matrix size
            uint imat = 20;

            Matrix< DDRMat > tA( imat, imat, 0.0 );

            // build non-symmetric matrix
            tA(0,0) = 2.0*imat;

            for (uint ir=1;ir<imat-1;++ir)
            {
                tA(ir,ir-1)=-0.9*imat;
                tA(ir,ir  )= 2.0*imat;
                tA(ir,ir+1)=-1.1*imat;
            }

            tA(imat-1,imat-2) =-1.0*imat;
            tA(imat-1,imat-1) = 2.0*imat;

            // compute square root - Option 1
            Matrix< DDRMat > tSqrt = moris::sqrtmat(tA);

            // check for accuracy
            real tError = norm( tA - tSqrt*tSqrt);

            REQUIRE( tError < 1e-8 * norm(tA) );

            // compute square root - Option 2
            tSqrt.fill(0.0);
            moris::sqrtmat(tA, tSqrt);

            // check for accuracy
            tError = norm( tA - tSqrt*tSqrt);

            REQUIRE( tError < 1e-8 * norm(tA) );
        }
    }
}

