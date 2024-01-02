/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Integration_Rule.cpp
 *
 */

#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/sr
#include "cl_MTK_Integrator.hpp" //MTK/sr

using namespace moris;
using namespace fem;

TEST_CASE("Integration_Rule_LINE", "[moris],[fem],[Integration_Rule_LINE]")
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // list of integration rules
    Vector< mtk::Integration_Order > tIntegrationOrderList = {
            mtk::Integration_Order::BAR_1,
            mtk::Integration_Order::BAR_2,
            mtk::Integration_Order::BAR_3,
            mtk::Integration_Order::BAR_4,
            mtk::Integration_Order::BAR_5,
            mtk::Integration_Order::BAR_6 };

    // create a space quad element
    Matrix<DDRMat > tXHat  = {{  0.0 }, { 1.0 }};
    Matrix<DDRMat > tXiHat = {{ -1.0 }, { 1.0 }};

    // reference volume
    real tSurfaceExact = 1.0;

    //create a line time element
    Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 }};
    Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

    // create a space and time geometry interpolation rule
    mtk::Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::LINE,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space and time geometry interpolator
    Geometry_Interpolator tGI( tGeoInterpRule );

    //set the coefficients xHat, tHat
    tGI.set_space_coeff( tXHat );
    tGI.set_time_coeff(  tTHat );

    // set the param coefficients xiHat, tauHat
    tGI.set_space_param_coeff( tXiHat );
    tGI.set_time_param_coeff(  tTauHat );

    // loop over integration rules
    for( uint iRule = 0; iRule < tIntegrationOrderList.size(); iRule++ )
    {
        // booleans for checks
        bool tSurfaceCheck = true;

        // get integration order
        mtk::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        mtk::Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::LINE,
                mtk::Integration_Type::GAUSS,
                tIntegrationOrder,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::BAR_1 );

        // create an integrator
        mtk::Integrator tIntegrator( tIntegrationRule );

        //get number of integration points, integration points and weights
        uint tNumOfIntegPoints =
                tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegrationPoints;
        tIntegrator.get_points( tIntegrationPoints );
        Matrix< DDRMat > tIntegrationWeights;
        tIntegrator.get_weights( tIntegrationWeights );

        // init the side surface
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in the surface ref space
            Matrix< DDRMat > tLocalIntegrationPoint =
                    tIntegrationPoints.get_column( iGP );

            // set local integration point
            tGI.set_space_time( tLocalIntegrationPoint );

            // evaluate surfDetJ
            real tSurfDetJ = tGI.det_J();

            // add contribution to the surface
            tSurface += tSurfDetJ * tIntegrationWeights( iGP );
        }

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - tSurfaceExact ) < tEpsilon );

        // check surfDetJ values
        REQUIRE( tSurfaceCheck );
    }
}

TEST_CASE("Integration_Rule_QUAD", "[moris],[fem],[Integration_Rule_QUAD]")
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // list of integration rules
    Vector< mtk::Integration_Order > tIntegrationOrderList = {
            //Integration_Order::QUAD_1x1,
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::QUAD_3x3,
            mtk::Integration_Order::QUAD_4x4,
            mtk::Integration_Order::QUAD_5x5 };

    // create a space quad element
    Matrix<DDRMat > tXHat  = {{  0.0,  0.0 }, { 1.0,  0.0 }, { 1.0, 1.0 }, {  0.0, 1.0 }};
    Matrix<DDRMat > tXiHat = {{ -1.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 }};

    // reference volume
    real tSurfaceExact = 1.0;

    //create a line time element
    Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 }};
    Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

    // create a space and time geometry interpolation rule
    mtk::Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space and time geometry interpolator
    Geometry_Interpolator tGI( tGeoInterpRule );

    //set the coefficients xHat, tHat
    tGI.set_space_coeff( tXHat );
    tGI.set_time_coeff(  tTHat );

    // set the param coefficients xiHat, tauHat
    tGI.set_space_param_coeff( tXiHat );
    tGI.set_time_param_coeff(  tTauHat );

    // loop over integration rules
    for( uint iRule = 0; iRule < tIntegrationOrderList.size(); iRule++ )
    {
        // booleans for checks
        bool tSurfaceCheck = true;

        // get integration order
        mtk::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        mtk::Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::QUAD,
                mtk::Integration_Type::GAUSS,
                tIntegrationOrder,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::BAR_1 );

        // create an integrator
        mtk::Integrator tIntegrator( tIntegrationRule );

        //get number of integration points, integration points and weights
        uint tNumOfIntegPoints =
                tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegrationPoints;
        tIntegrator.get_points( tIntegrationPoints );
        Matrix< DDRMat > tIntegrationWeights;
        tIntegrator.get_weights( tIntegrationWeights );

        // init the side surface
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in the surface ref space
            Matrix< DDRMat > tLocalIntegrationPoint =
                    tIntegrationPoints.get_column( iGP );

            // set local integration point
            tGI.set_space_time( tLocalIntegrationPoint );

            // evaluate surfDetJ
            real tSurfDetJ = tGI.det_J();

            // add contribution to the surface
            tSurface += tSurfDetJ * tIntegrationWeights( iGP );
        }

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - tSurfaceExact ) < tEpsilon );

        // check surfDetJ values
        REQUIRE( tSurfaceCheck );
    }
}

TEST_CASE("Integration_Rule_HEX", "[moris],[fem],[Integration_Rule_HEX]")
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // list of integration rules
    Vector< mtk::Integration_Order > tIntegrationOrderList = {
            //Integration_Order::HEX_1x1x1,
            mtk::Integration_Order::HEX_2x2x2,
            mtk::Integration_Order::HEX_3x3x3,
            mtk::Integration_Order::HEX_4x4x4,
            mtk::Integration_Order::HEX_5x5x5 };

    // create a space hex element
    Matrix< DDRMat > tXHat = {
            { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 }};
    Matrix< DDRMat > tXiHat = {
            { -1.0, -1.0, -1.0 }, {  1.0, -1.0, -1.0 }, {  1.0,  1.0, -1.0 }, { -1.0,  1.0, -1.0 },
            { -1.0, -1.0,  1.0 }, {  1.0, -1.0,  1.0 }, {  1.0,  1.0,  1.0 }, { -1.0,  1.0,  1.0 }};

    // reference volume
    real tSurfaceExact = 1.0;

    //create a line time element
    Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 }};
    Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

    // create a space and time geometry interpolation rule
    mtk::Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space and time geometry interpolator
    Geometry_Interpolator tGI( tGeoInterpRule );

    //set the coefficients xHat, tHat
    tGI.set_space_coeff( tXHat );
    tGI.set_time_coeff(  tTHat );

    // set the param coefficients xiHat, tauHat
    tGI.set_space_param_coeff( tXiHat );
    tGI.set_time_param_coeff(  tTauHat );

    // loop over integration rules
    for( uint iRule = 0; iRule < tIntegrationOrderList.size(); iRule++ )
    {
        // booleans for checks
        bool tSurfaceCheck = true;

        // get integration order
        mtk::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        mtk::Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::HEX,
                mtk::Integration_Type::GAUSS,
                tIntegrationOrder,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::BAR_1 );

        // create an integrator
        mtk::Integrator tIntegrator( tIntegrationRule );

        //get number of integration points, integration points and weights
        uint tNumOfIntegPoints =
                tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegrationPoints;
        tIntegrator.get_points( tIntegrationPoints );
        Matrix< DDRMat > tIntegrationWeights;
        tIntegrator.get_weights( tIntegrationWeights );

        // init the side surface
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in the surface ref space
            Matrix< DDRMat > tLocalIntegrationPoint =
                    tIntegrationPoints.get_column( iGP );

            // set local integration point
            tGI.set_space_time( tLocalIntegrationPoint );

            // evaluate surfDetJ
            real tSurfDetJ = tGI.det_J();

            // add contribution to the surface
            tSurface += tSurfDetJ * tIntegrationWeights( iGP );
        }

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - tSurfaceExact ) < tEpsilon );

        // check surfDetJ values
        REQUIRE( tSurfaceCheck );
    }
}

TEST_CASE("Integration_Rule_TRI", "[moris],[fem],[Integration_Rule_TRI]")
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // list of integration rules
    Vector< mtk::Integration_Order > tIntegrationOrderList = {
            mtk::Integration_Order::TRI_1,
            mtk::Integration_Order::TRI_3,
            mtk::Integration_Order::TRI_4,
            mtk::Integration_Order::TRI_6,
            mtk::Integration_Order::TRI_7,
            mtk::Integration_Order::TRI_12,
            mtk::Integration_Order::TRI_13,
            mtk::Integration_Order::TRI_16,
            mtk::Integration_Order::TRI_19,
            mtk::Integration_Order::TRI_25 };

    // create a space quad element
    Matrix<DDRMat > tXHat  = {{ 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }};
    Matrix<DDRMat > tXiHat = {{ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }};

    // reference volume
    real tSurfaceExact = 0.5;

    //create a line time element
    Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 }};
    Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

    // create a space and time geometry interpolation rule
    mtk::Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::TRI,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space and time geometry interpolator
    Geometry_Interpolator tGI( tGeoInterpRule );

    //set the coefficients xHat, tHat
    tGI.set_space_coeff( tXHat );
    tGI.set_time_coeff(  tTHat );

    // set the param coefficients xiHat, tauHat
    tGI.set_space_param_coeff( tXiHat );
    tGI.set_time_param_coeff(  tTauHat );

    // loop over integration rules
    for( uint iRule = 0; iRule < tIntegrationOrderList.size(); iRule++ )
    {
        // booleans for checks
        bool tSurfaceCheck = true;

        // get integration order
        mtk::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        mtk::Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::TRI,
                mtk::Integration_Type::GAUSS,
                tIntegrationOrder,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::BAR_1 );

        // create an integrator
        mtk::Integrator tIntegrator( tIntegrationRule );

        //get number of integration points, integration points and weights
        uint tNumOfIntegPoints =
                tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegrationPoints;
        tIntegrator.get_points( tIntegrationPoints );
        Matrix< DDRMat > tIntegrationWeights;
        tIntegrator.get_weights( tIntegrationWeights );

        // init the side surface
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in the surface ref space
            Matrix< DDRMat > tLocalIntegrationPoint =
                    tIntegrationPoints.get_column( iGP );

            // set local integration point
            tGI.set_space_time( tLocalIntegrationPoint );

            // evaluate surfDetJ
            real tSurfDetJ = tGI.det_J();

            // add contribution to the surface
            tSurface += tSurfDetJ * tIntegrationWeights( iGP );
        }

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - tSurfaceExact ) < tEpsilon );

        // check surfDetJ values
        REQUIRE( tSurfaceCheck );
    }
}

TEST_CASE("Integration_Rule_TET", "[moris],[fem],[Integration_Rule_TET]")
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // list of integration rules
    Vector< mtk::Integration_Order > tIntegrationOrderList = {
            //Integration_Order::TET_1,
            mtk::Integration_Order::TET_4,
            mtk::Integration_Order::TET_5,
            mtk::Integration_Order::TET_11,
            mtk::Integration_Order::TET_15,
            mtk::Integration_Order::TET_20,
            mtk::Integration_Order::TET_35,
            mtk::Integration_Order::TET_56 };

    // create a space quad element
    Matrix< DDRMat > tXHat = {{ 0.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0,  0.0, 1.0 }};
    Matrix< DDRMat > tXiHat = {
            { 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000 }};

    // reference volume
    real tSurfaceExact = 1.666666666666667e-01;

    //create a line time element
    Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 }};
    Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

    // create a space and time geometry interpolation rule
    mtk::Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::TET,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space and time geometry interpolator
    Geometry_Interpolator tGI( tGeoInterpRule );

    //set the coefficients xHat, tHat
    tGI.set_space_coeff( tXHat );
    tGI.set_time_coeff(  tTHat );

    // set the param coefficients xiHat, tauHat
    tGI.set_space_param_coeff( tXiHat );
    tGI.set_time_param_coeff(  tTauHat );

    // loop over integration rules
    for( uint iRule = 0; iRule < tIntegrationOrderList.size(); iRule++ )
    {
        // booleans for checks
        bool tSurfaceCheck = true;

        // get integration order
        mtk::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        mtk::Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::TET,
                mtk::Integration_Type::GAUSS,
                tIntegrationOrder,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::BAR_1 );

        // create an integrator
        mtk::Integrator tIntegrator( tIntegrationRule );

        //get number of integration points, integration points and weights
        uint tNumOfIntegPoints =
                tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegrationPoints;
        tIntegrator.get_points( tIntegrationPoints );
        Matrix< DDRMat > tIntegrationWeights;
        tIntegrator.get_weights( tIntegrationWeights );

        // init the side surface
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in the surface ref space
            Matrix< DDRMat > tLocalIntegrationPoint =
                    tIntegrationPoints.get_column( iGP );

            // set local integration point
            tGI.set_space_time( tLocalIntegrationPoint );

            // evaluate surfDetJ
            real tSurfDetJ = tGI.det_J();

            // add contribution to the surface
            tSurface += tSurfDetJ * tIntegrationWeights( iGP );
        }

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - tSurfaceExact ) < tEpsilon );

        // check surfDetJ values
        REQUIRE( tSurfaceCheck );
    }
}
