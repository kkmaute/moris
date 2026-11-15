/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Geometry_Interpolator.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"

#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "fn_norm.hpp"

namespace moris::fem
{

    TEST_CASE( "GI", "[moris],[fem],[GI_update]" )
    {
        // create evaluation point xi, tau
        //------------------------------------------------------------------------------
        Matrix< DDRMat > tParamPoint = { { 1.0 }, { 0.25 }, { 0.1 } };

        // space and time geometry interpolators
        //------------------------------------------------------------------------------
        // create a space geometry interpolation rule
        mtk::Interpolation_Rule tGIRule(
                mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create space coeff xHat
        Matrix< DDRMat > tXHat = { { 0.0, 0.0 },
            { 1.0, 0.0 },
            { 1.0, 1.0 },
            { 0.0, 1.0 } };

        // create time coeff tHat
        Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set the evaluation point
        tGI.set_space_time( tParamPoint );

        // update local coordinates
        Matrix< DDRMat > tPhysCoordinates  = { { 0.75, 0.25 } };
        Matrix< DDRMat > tParamCoordinates = { { -0.1, 0.34 } };
        tGI.update_parametric_coordinates(
                tPhysCoordinates,
                tParamCoordinates );

        // print( tParamCoordinates, "tParamCoordinates" );

    } /* END_TEST_CASE */

    TEST_CASE( "GI_Current_Coordinates_Preserve_Time", "[moris],[fem],[GI_update]" )
    {
        // Time-dependent displacements must be sampled at the current field time,
        // and evaluating the nodal coordinates must restore the field's evaluation point.
        mtk::Interpolation_Rule tRule(
                mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        Geometry_Interpolator tGI( tRule );
        Matrix< DDRMat >      tXHat = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };
        tGI.set_coeff( tXHat, Matrix< DDRMat >{ { 0.0 }, { 1.0 } } );
        tGI.set_param_coeff();

        Field_Interpolator tFI( 2, tRule, &tGI, { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        Matrix< DDRMat >   tCoefficients( 8, 2 );
        for ( uint iNode = 0; iNode < 4; ++iNode )
        {
            real const tOffset            = static_cast< real >( iNode );
            tCoefficients( iNode, 0 )     = 1.0 + tOffset;
            tCoefficients( iNode, 1 )     = -2.0 - tOffset;
            tCoefficients( iNode + 4, 0 ) = 3.0 + tOffset;
            tCoefficients( iNode + 4, 1 ) = 4.0 - tOffset;
        }
        tFI.set_coeff( tCoefficients );

        for ( real tTau : { -0.5, 0.5 } )
        {
            Matrix< DDRMat > tPoint = { { 0.2 }, { -0.3 }, { tTau } };
            tGI.set_space_time( tPoint );
            tGI.reset_eval_flags_coordinates();
            tFI.set_space_time( tPoint );
            Matrix< DDRMat > const tValueBefore = tFI.val();

            Matrix< DDRMat > const tCurrentCoordinates = tGI.get_space_coeff_current( &tFI );
            for ( uint iNode = 0; iNode < 4; ++iNode )
            {
                real const tOffset = static_cast< real >( iNode );
                REQUIRE( tCurrentCoordinates( iNode, 0 ) == Approx( tXHat( iNode, 0 ) + 2.0 + tTau + tOffset ) );
                REQUIRE( tCurrentCoordinates( iNode, 1 ) == Approx( tXHat( iNode, 1 ) + 1.0 + 3.0 * tTau - tOffset ) );
            }

            REQUIRE( norm( tFI.get_space_time() - tPoint ) < 1.0e-12 );
            REQUIRE( norm( tFI.val() - tValueBefore ) < 1.0e-12 );
        }
    }
}    // namespace moris::fem
