/*
 * UT_FEM_Geometry_Interpolator.cpp
 *
 *  Created on: May 18, 2020
 *      Author: noel
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"

#include "cl_FEM_Geometry_Interpolator.hpp"

namespace moris
{
    namespace fem
    {

        TEST_CASE( "GI", "[moris],[fem],[GI_update]" )
        {
            // create evaluation point xi, tau
            //------------------------------------------------------------------------------
            Matrix< DDRMat > tParamPoint = {{ 1.0 }, { 0.25 }, { 0.1 }};

            // space and time geometry interpolators
            //------------------------------------------------------------------------------
            // create a space geometry interpolation rule
            Interpolation_Rule tGIRule(
                    mtk::Geometry_Type::QUAD,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // create a space time geometry interpolator
            Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

            // create space coeff xHat
            Matrix< DDRMat > tXHat = {{ 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0}};

            // create time coeff tHat
            Matrix< DDRMat > tTHat = {{ 0.0 }, {1.0 }};

            // set the coefficients xHat, tHat
            tGI.set_coeff( tXHat, tTHat );

            // set the evaluation point
            tGI.set_space_time( tParamPoint );

            // update local coordinates
            Matrix< DDRMat > tPhysCoordinates  = {{ 0.75, 0.25 }};
            Matrix< DDRMat > tParamCoordinates = {{ -0.1,  0.34  }};
            tGI.update_local_coordinates(
                    tPhysCoordinates,
                    tParamCoordinates );

            print( tParamCoordinates, "tParamCoordinates" );

        }/* END_TEST_CASE */
    }
}
