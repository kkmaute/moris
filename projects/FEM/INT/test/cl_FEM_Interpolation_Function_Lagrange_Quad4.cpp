/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Interpolation_Function_Lagrange_Quad4.cpp
 *
 */

#include <catch.hpp>

#include "moris_typedefs.hpp" //MRS/COR/src
#include "paths.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_save_matrix_to_binary_file.hpp" //LNA/src
#include "fn_load_matrix_from_binary_file.hpp" //LNA/src
#include "op_times.hpp" //LNA/src
#include "op_minus.hpp" //LNA/src
#include "fn_trans.hpp" //LNA/src
#include "fn_norm.hpp" //LNA/src
#include "fn_sum.hpp" //LNA/src

#include "cl_MTK_Interpolation_Rule.hpp" //MTK/src
#include "fn_FEM_Check.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "Lagrange QUAD4", "[moris],[fem],[Quad4LagInterpolation]" )
{

    //------------------------------------------------------------------------------

    // step 1: load MATLAB precomputed data from binary files
    std::string tPrefix = moris::get_base_moris_dir();
    tPrefix = tPrefix + "/projects/FEM/INT/test/data/" ;

    // load point coordinates from file
    Matrix< DDRMat > tXi;
    load_matrix_from_binary_file( tXi,
            tPrefix + "points_2d.bin" );

    // load values from nodes from file
    Matrix< DDRMat > tPhiHat;
    load_matrix_from_binary_file( tPhiHat,
            tPrefix + "lagrange_quad4_phihat.bin" );

    // load solutions for N*tPhiHat
    Matrix< DDRMat > tPhi;
    load_matrix_from_binary_file( tPhi,
            tPrefix + "lagrange_quad4_phi.bin" );

    // load solutions for dNdXi*tPhiHat
    Matrix< DDRMat > tdPhidXi;
    load_matrix_from_binary_file( tdPhidXi,
            tPrefix + "lagrange_quad4_dphidxi.bin" );

    // load solutions for d2NdXi2*tPhiHat
    Matrix< DDRMat > td2PhidXi2;
    load_matrix_from_binary_file( td2PhidXi2,
            tPrefix + "lagrange_quad4_d2phidxi2.bin" );

    //------------------------------------------------------------------------------

    // step 2 create function and interpolation matrices

    // create rule
    mtk::Interpolation_Rule tRule( mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::CONSTANT,
            mtk::Interpolation_Order::CONSTANT );

    // create shape function object
    auto tFunction = tRule.create_space_interpolation_function();

    // create matrix that contains the shape function
    Matrix< DDRMat > tN;

    // create matrix that contains the first derivative
    Matrix< DDRMat > tdNdXi;

    // create matrix that contains the second derivative
    Matrix< DDRMat > td2NdXi2;

    // create matrix that contains the third derivative
    Matrix< DDRMat > td3NdXi3;

    //------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // get number of points to test
    auto tNumberOfTestPoints = tXi.n_cols();

    //------------------------------------------------------------------------------

    SECTION( "QUAD4: test for unity" )
    {
        bool tCheck = true;
        for( uint k=0; k<tNumberOfTestPoints; ++k )
        {
            // evaluate shape function at point k
            tFunction->eval_N( tXi.get_column(k ), tN );

            // test unity
            tCheck = tCheck && ( std::abs( sum(tN) - 1.0 ) < tEpsilon );
        }

        REQUIRE( tCheck );
    }

    //------------------------------------------------------------------------------

    SECTION( "QUAD4: test N" )
    {
        bool tCheck = true;
        for( uint k=0; k<tNumberOfTestPoints; ++k )
        {
            // evaluate shape function at point k
            tFunction->eval_N( tXi.get_column( k ), tN );

            // test evaluated value
            Matrix< DDRMat > tError  = tN * tPhiHat ;
            tError( 0 ) -= tPhi( k );

            // test error
            tCheck = tCheck && ( norm(tError) < tEpsilon );
        }

        REQUIRE( tCheck );
    }

    //------------------------------------------------------------------------------

    SECTION( "QUAD4: test dNdXi" )
    {
        bool tCheck = true;
        for( uint k=0; k<tNumberOfTestPoints; ++k )
        {
            // td2NdXi2evaluate shape function at point k
            tFunction->eval_dNdXi( tXi.get_column(k ), tdNdXi );

            // test evaluated value
            Matrix< DDRMat > tError = tdPhidXi.get_column( k );
            tError = tError - tdNdXi*tPhiHat;

            // test error
            tCheck = tCheck && ( norm(tError) < tEpsilon );
        }

        REQUIRE( tCheck );
    }

    //------------------------------------------------------------------------------

    SECTION( "QUAD4: test d2NdXi2" )
    {
        bool tCheck = true;
        for( uint k=0; k<tNumberOfTestPoints; ++k )
        {
            // evaluate shape function at point k
            tFunction->eval_d2NdXi2( tXi.get_column( k ), td2NdXi2 );

            // test evaluated valueN
            Matrix< DDRMat > tError = td2PhidXi2.get_column( k );
            tError = tError - td2NdXi2*tPhiHat;

            // test error
            tCheck = tCheck && ( norm(tError) < tEpsilon );
        }

        REQUIRE( tCheck );
    }

    //------------------------------------------------------------------------------

    SECTION( "QUAD4: test d3NdXi3" )
    {
        Matrix< DDRMat > td3NdXi3FD( 4, 4, 0.0 );

        bool tCheck = true;
        for( uint k=0; k<tNumberOfTestPoints; ++k )
        {
            // get evaluation point
            Matrix< DDRMat > tEvalPoint = tXi.get_column( k );

            // evaluate shape function at point k
            tFunction->eval_d3NdXi3( tEvalPoint, td3NdXi3 );

            for ( uint l=0; l<tEvalPoint.numel(); l++ )
            {
                // perturbation size
                real tDeltaXi = 1e-6;

                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tEvalPoint;
                tPertEvalPoint( l ) = tPertEvalPoint( l ) - tDeltaXi;

                // evaluate td2NdXi2 at point k - delta xi
                Matrix< DDRMat > td2NdXi2Minus;
                tFunction->eval_d2NdXi2( tPertEvalPoint, td2NdXi2Minus );

                // perturbed evaluation point
                tPertEvalPoint = tEvalPoint;
                tPertEvalPoint( l ) = tPertEvalPoint( l ) + tDeltaXi;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > td2NdXi2Plus;
                tFunction->eval_d2NdXi2( tPertEvalPoint, td2NdXi2Plus );

                // evaluate FD
                Matrix< DDRMat > td3NdXi3FDTemp = ( td2NdXi2Plus - td2NdXi2Minus ) / ( 2.0 * tDeltaXi );

                // evaluate td3NdXi3_FD
                if( l == 0 )
                {
                    td3NdXi3FD.get_row( 0 ) = td3NdXi3FDTemp.get_row( 0 );
                    td3NdXi3FD.get_row( 2 ) = td3NdXi3FDTemp.get_row( 2 );
                    td3NdXi3FD.get_row( 3 ) = td3NdXi3FDTemp.get_row( 1 );
                }
                else
                {
                    td3NdXi3FD.get_row( 1 ) = td3NdXi3FDTemp.get_row( 1 );
                }
            }

            // test error
            tCheck = tCheck && fem::check( td3NdXi3, td3NdXi3FD, 1e-9 );
        }

        REQUIRE( tCheck );
    }

    //------------------------------------------------------------------------------

    // tidy up
    delete tFunction;

    //------------------------------------------------------------------------------
}

