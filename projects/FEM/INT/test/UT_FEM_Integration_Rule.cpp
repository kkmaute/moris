#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/sr
#include "cl_FEM_Integrator.hpp" //FEM/INT/sr

using namespace moris;
using namespace fem;

TEST_CASE("Integration_Rule_LINE", "[moris],[fem],[Integration_Rule_LINE]")
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // list of integration rules
    moris::Cell< fem::Integration_Order > tIntegrationOrderList = {
            Integration_Order::BAR_1,
            Integration_Order::BAR_2,
            Integration_Order::BAR_3,
            Integration_Order::BAR_4,
            Integration_Order::BAR_5,
            Integration_Order::BAR_6 };

    // create a space quad element
    Matrix<DDRMat > tXHat  = {{  0.0 }, { 1.0 }};
    Matrix<DDRMat > tXiHat = {{ -1.0 }, { 1.0 }};

    // reference volume
    real tSurfaceExact = 1.0;

    //create a line time element
    Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 }};
    Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

    // create a space and time geometry interpolation rule
    Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::LINE,
            Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            Interpolation_Type::LAGRANGE,
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
        fem::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::LINE,
                Integration_Type::GAUSS,
                tIntegrationOrder,
                Integration_Type::GAUSS,
                Integration_Order::BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

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

            // map integration point
            const Matrix< DDRMat > & tGlobalIntegrationPoint = tGI.map_integration_point();

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
    moris::Cell< fem::Integration_Order > tIntegrationOrderList = {
            //Integration_Order::QUAD_1x1,
            Integration_Order::QUAD_2x2,
            Integration_Order::QUAD_3x3,
            Integration_Order::QUAD_4x4,
            Integration_Order::QUAD_5x5 };

    // create a space quad element
    Matrix<DDRMat > tXHat  = {{  0.0,  0.0 }, { 1.0,  0.0 }, { 1.0, 1.0 }, {  0.0, 1.0 }};
    Matrix<DDRMat > tXiHat = {{ -1.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 }};

    // reference volume
    real tSurfaceExact = 1.0;

    //create a line time element
    Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 }};
    Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

    // create a space and time geometry interpolation rule
    Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::QUAD,
            Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            Interpolation_Type::LAGRANGE,
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
        fem::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::QUAD,
                Integration_Type::GAUSS,
                tIntegrationOrder,
                Integration_Type::GAUSS,
                Integration_Order::BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

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

            // map integration point
            const Matrix< DDRMat > & tGlobalIntegrationPoint = tGI.map_integration_point();

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
    moris::Cell< fem::Integration_Order > tIntegrationOrderList = {
            //Integration_Order::HEX_1x1x1,
            Integration_Order::HEX_2x2x2,
            Integration_Order::HEX_3x3x3,
            Integration_Order::HEX_4x4x4,
            Integration_Order::HEX_5x5x5 };

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
    Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::HEX,
            Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            Interpolation_Type::LAGRANGE,
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
        fem::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::HEX,
                Integration_Type::GAUSS,
                tIntegrationOrder,
                Integration_Type::GAUSS,
                Integration_Order::BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

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

            // map integration point
            const Matrix< DDRMat > & tGlobalIntegrationPoint = tGI.map_integration_point();

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
    moris::Cell< fem::Integration_Order > tIntegrationOrderList = {
            Integration_Order::TRI_1,
            Integration_Order::TRI_3,
            Integration_Order::TRI_4,
            Integration_Order::TRI_6,
            Integration_Order::TRI_7,
            Integration_Order::TRI_12,
            Integration_Order::TRI_13,
            Integration_Order::TRI_16,
            Integration_Order::TRI_19,
            Integration_Order::TRI_25 };

    // create a space quad element
    Matrix<DDRMat > tXHat  = {{ 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }};
    Matrix<DDRMat > tXiHat = {{ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }};

    // reference volume
    real tSurfaceExact = 0.5;

    //create a line time element
    Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 }};
    Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

    // create a space and time geometry interpolation rule
    Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::TRI,
            Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            Interpolation_Type::LAGRANGE,
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
        fem::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::TRI,
                Integration_Type::GAUSS,
                tIntegrationOrder,
                Integration_Type::GAUSS,
                Integration_Order::BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

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

            // map integration point
            const Matrix< DDRMat > & tGlobalIntegrationPoint = tGI.map_integration_point();

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
    moris::Cell< fem::Integration_Order > tIntegrationOrderList = {
            //Integration_Order::TET_1,
            Integration_Order::TET_4,
            Integration_Order::TET_5,
            Integration_Order::TET_11,
            Integration_Order::TET_15,
            Integration_Order::TET_20,
            Integration_Order::TET_35,
            Integration_Order::TET_56 };

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
    Interpolation_Rule tGeoInterpRule(
            mtk::Geometry_Type::TET,
            Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            Interpolation_Type::LAGRANGE,
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
        fem::Integration_Order tIntegrationOrder =
                tIntegrationOrderList( iRule );

        // create an integration rule
        Integration_Rule tIntegrationRule(
                mtk::Geometry_Type::TET,
                Integration_Type::GAUSS,
                tIntegrationOrder,
                Integration_Type::GAUSS,
                Integration_Order::BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

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

            // map integration point
            const Matrix< DDRMat > & tGlobalIntegrationPoint = tGI.map_integration_point();

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
