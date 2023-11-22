/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Integrator.cpp
 *
 */

#include <string>

#include <catch.hpp>

#include "moris_typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_MTK_Integration_Rule.hpp" //MTK/src
#include "cl_MTK_Integrator.hpp" //MTK/src
#include "cl_MTK_Integrator_Test_Polynomial.hpp" //MTK/src

#include "paths.hpp"
using namespace moris;
using namespace fem;

TEST_CASE( "Integrator", "[moris],[fem]" )
{
//------------------------------------------------------------------------------

    // define epsilon environment
    const real tEpsilon = 1e-12;

    // define path of files
    std::string tPrefix = moris::get_base_moris_dir();
    tPrefix = tPrefix + "/projects/FEM/INT/test/data/" ;

//------------------------------------------------------------------------------

    SECTION( "GAUSS QUAD_2x2" )
    {

        // initialize polynomial
        mtk::Integrator_Test_Polynomial tPoly(
                tPrefix+ "integrate_quad_2x2.bin" );

        // create integration rule
        mtk::Integration_Rule tRule(
                mtk::Geometry_Type::QUAD,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::QUAD_2x2 );

        // create integrator
        mtk::Integrator tIntegrator( tRule );

        // get number of points
        uint tNumberOfPoints = tIntegrator.get_number_of_points();

        // get points
        Matrix< DDRMat > tPoints;
        tIntegrator.get_points( tPoints );

        // get weights
        Matrix< DDRMat > tWeights;
        tIntegrator.get_weights( tWeights );

        // initialize value
        real tValue = 0.0;

        // loop over all points
        for( uint k=0; k<tNumberOfPoints; ++k )
        {
            tValue += tWeights( k ) * tPoly.eval( tPoints.get_column( k ) );
        }

        // calculate error
        real tError = std::abs( tPoly.get_integral() - tValue );

        // perform test
        REQUIRE( tError < tEpsilon );
    }

//------------------------------------------------------------------------------

    SECTION( "GAUSS QUAD_3x3" )
    {
        // initialize polynomial
        mtk::Integrator_Test_Polynomial tPoly(
                tPrefix+ "integrate_quad_3x3.bin");

        // create integration rule
        mtk::Integration_Rule tRule(
                mtk::Geometry_Type::QUAD,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::QUAD_3x3 );

        // create integrator
        mtk::Integrator tIntegrator( tRule );

        // get number of points
        uint tNumberOfPoints = tIntegrator.get_number_of_points();

        // get points
        Matrix< DDRMat > tPoints;
        tIntegrator.get_points( tPoints );

        // get weights
        Matrix< DDRMat > tWeights;
        tIntegrator.get_weights( tWeights );

        // initialize value
        real tValue = 0.0;

        // loop over all points
        for( uint k=0; k<tNumberOfPoints; ++k )
        {
            tValue += tWeights( k ) * tPoly.eval( tPoints.get_column( k ) );
        }

        // calculate error
        real tError = std::abs( tPoly.get_integral() - tValue );

        // perform test
        REQUIRE( tError < tEpsilon );
    }

//------------------------------------------------------------------------------

    SECTION( "GAUSS QUAD_4x4" )
    {
        // initialize polynomial
        mtk::Integrator_Test_Polynomial tPoly(
                tPrefix+ "integrate_quad_4x4.bin");

        // create integration rule
        mtk::Integration_Rule tRule(
                mtk::Geometry_Type::QUAD,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::QUAD_4x4 );

        // create integrator
        mtk::Integrator tIntegrator( tRule );

        // get number of points
        uint tNumberOfPoints = tIntegrator.get_number_of_points();

        // get points
        Matrix< DDRMat > tPoints;
        tIntegrator.get_points( tPoints );

        // get weights
        Matrix< DDRMat > tWeights;
        tIntegrator.get_weights( tWeights );

        // initialize value
        real tValue = 0.0;

        // loop over all points
        for( uint k=0; k<tNumberOfPoints; ++k )
        {
            tValue += tWeights( k ) * tPoly.eval( tPoints.get_column( k ) );
        }

        // calculate error
        real tError = std::abs( tPoly.get_integral() - tValue );

        // perform test
        REQUIRE( tError < tEpsilon );
    }

//------------------------------------------------------------------------------

    SECTION( "GAUSS QUAD_5x5" )
    {
        // initialize polynomial
        mtk::Integrator_Test_Polynomial tPoly(
                tPrefix+ "integrate_quad_5x5.bin");

        // create integration rule
        mtk::Integration_Rule tRule(
                mtk::Geometry_Type::QUAD,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::QUAD_5x5 );

        // create integrator
        mtk::Integrator tIntegrator( tRule );

        // get number of points
        uint tNumberOfPoints = tIntegrator.get_number_of_points();

        // get points
        Matrix< DDRMat > tPoints;
        tIntegrator.get_points( tPoints );

        // get weights
        Matrix< DDRMat > tWeights;
        tIntegrator.get_weights( tWeights );

        // initialize value
        real tValue = 0.0;

        // loop over all points
        for( uint k=0; k<tNumberOfPoints; ++k )
        {
            tValue += tWeights( k ) * tPoly.eval( tPoints.get_column( k ) );
        }

        // calculate error
        real tError = std::abs( tPoly.get_integral() - tValue );

        // perform test
        REQUIRE( tError < tEpsilon );
    }

//------------------------------------------------------------------------------

    SECTION( "GAUSS HEX_2x2x2" )
    {
        // initialize polynomial
        mtk::Integrator_Test_Polynomial tPoly(
                tPrefix+ "integrate_hex_2x2x2.bin");

        // create integration rule
        mtk::Integration_Rule tRule(
                mtk::Geometry_Type::HEX,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::HEX_2x2x2 );

        // create integrator
        mtk::Integrator tIntegrator( tRule );

        // get number of points
        uint tNumberOfPoints = tIntegrator.get_number_of_points();

        // get points
        Matrix< DDRMat > tPoints;
        tIntegrator.get_points( tPoints );

        // get weights
        Matrix< DDRMat > tWeights;
        tIntegrator.get_weights( tWeights );

        // initialize value
        real tValue = 0.0;

        // loop over all points
        for( uint k=0; k<tNumberOfPoints; ++k )
        {
            tValue += tWeights( k ) * tPoly.eval( tPoints.get_column( k ) );
        }

        // calculate error
        real tError = std::abs( tPoly.get_integral() - tValue );

        // perform test
        REQUIRE( tError < tEpsilon );
    }

//------------------------------------------------------------------------------

    SECTION( "GAUSS HEX_3x3x3" )
    {
        // initialize polynomial
        mtk::Integrator_Test_Polynomial tPoly(
                tPrefix+ "integrate_hex_3x3x3.bin");

        // create integration rule
        mtk::Integration_Rule tRule(
                mtk::Geometry_Type::HEX,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::HEX_3x3x3 );

        // create integrator
        mtk::Integrator tIntegrator( tRule );

        // get number of points
        uint tNumberOfPoints = tIntegrator.get_number_of_points();

        // get points
        Matrix< DDRMat > tPoints;
        tIntegrator.get_points( tPoints );

        // get weights
        Matrix< DDRMat > tWeights;
        tIntegrator.get_weights( tWeights );

        // initialize value
        real tValue = 0.0;

        // loop over all points
        for( uint k=0; k<tNumberOfPoints; ++k )
        {
            tValue += tWeights( k ) * tPoly.eval( tPoints.get_column( k ) );
        }

        // calculate error
        real tError = std::abs( tPoly.get_integral() - tValue );

        // perform test
        REQUIRE( tError < tEpsilon );
    }

//------------------------------------------------------------------------------

    SECTION( "GAUSS HEX_4x4x4" )
    {
        // initialize polynomial
        mtk::Integrator_Test_Polynomial tPoly(
                tPrefix+ "integrate_hex_4x4x4.bin");

        // create integration rule
        mtk::Integration_Rule tRule(
                mtk::Geometry_Type::HEX,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::HEX_4x4x4 );

        // create integrator
        mtk::Integrator tIntegrator( tRule );

        // get number of points
        uint tNumberOfPoints = tIntegrator.get_number_of_points();

        // get points
        Matrix< DDRMat > tPoints;
        tIntegrator.get_points( tPoints );

        // get weights
        Matrix< DDRMat > tWeights;
        tIntegrator.get_weights( tWeights );

        // initialize value
        real tValue = 0.0;

        // loop over all points
        for( uint k=0; k<tNumberOfPoints; ++k )
        {
            tValue += tWeights( k ) * tPoly.eval( tPoints.get_column( k ) );
        }

        // calculate error
        real tError = std::abs( tPoly.get_integral() - tValue );

        // perform test
        REQUIRE( tError < tEpsilon );
    }

//------------------------------------------------------------------------------

    SECTION( "GAUSS HEX_5x5x5" )
    {
        // initialize polynomial
        mtk::Integrator_Test_Polynomial tPoly(
                tPrefix+ "integrate_hex_5x5x5.bin");

        // create integration rule
        mtk::Integration_Rule tRule(
                mtk::Geometry_Type::HEX,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::HEX_5x5x5 );

        // create integrator
        mtk::Integrator tIntegrator( tRule );

        // get number of points
        uint tNumberOfPoints = tIntegrator.get_number_of_points();

        // get points
        Matrix< DDRMat > tPoints;
        tIntegrator.get_points( tPoints );

        // get weights
        Matrix< DDRMat > tWeights;
        tIntegrator.get_weights( tWeights );

        // initialize value
        real tValue = 0.0;

        // loop over all points
        for( uint k=0; k<tNumberOfPoints; ++k )
        {
            tValue += tWeights( k ) * tPoly.eval( tPoints.get_column( k ) );
        }

        // calculate error
        real tError = std::abs( tPoly.get_integral() - tValue );

        // perform test
        REQUIRE( tError < tEpsilon );
    }

//------------------------------------------------------------------------------
}
