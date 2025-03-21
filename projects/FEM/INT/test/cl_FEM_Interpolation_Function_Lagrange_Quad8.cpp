/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Interpolation_Function_Lagrange_Quad8.cpp
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
#include "fn_sum.hpp" //LNA/src
#include "fn_trans.hpp" //LNA/src
#include "fn_norm.hpp" //LNA/src
#include "cl_FEM_Enums.hpp"

#include "cl_MTK_Interpolation_Rule.hpp" //MTK/src

using namespace moris;
using namespace fem;

TEST_CASE( "Lagrange QUAD8", "[moris],[fem],[Quad8LagInterpolation]" )
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
            tPrefix + "lagrange_quad8_phihat.bin" );

    // load solutions for N*tPhiHat
    Matrix< DDRMat > tPhi;
    load_matrix_from_binary_file( tPhi,
            tPrefix + "lagrange_quad8_phi.bin" );

    // load solutions for dNdXi*tPhiHat
    Matrix< DDRMat > tdPhidXi;
    load_matrix_from_binary_file( tdPhidXi,
            tPrefix + "lagrange_quad8_dphidxi.bin" );

    // load solutions for d2NdXi2*tPhiHat
    Matrix< DDRMat > td2PhidXi2;
    load_matrix_from_binary_file( td2PhidXi2,
            tPrefix + "lagrange_quad8_d2phidxi2.bin" );

    //------------------------------------------------------------------------------

    // step 2 create function and interpolation matrices

    // create rule
    mtk::Interpolation_Rule tRule( mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::SERENDIPITY,
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

    //------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // get number of points to test
    auto tNumberOfTestPoints = tXi.n_cols();

    //------------------------------------------------------------------------------

    SECTION( "QUAD8: test for unity" )
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

    SECTION( "QUAD8: test N" )
    {
        bool tCheck = true;
        for( uint k=0; k<tNumberOfTestPoints; ++k )
        {
            // evaluate shape function at point k
            tFunction->eval_N( tXi.get_column(k ), tN );

            // test evaluated value
            Matrix< DDRMat > tError  = tN * tPhiHat ;
            tError( 0 ) -= tPhi( k );

            // test error
            tCheck = tCheck && ( norm(tError) < tEpsilon );
        }

        REQUIRE( tCheck );
    }

    //------------------------------------------------------------------------------

    SECTION( "QUAD8: test dNdXi" )
    {
        bool tCheck = true;
        for( uint k=0; k<tNumberOfTestPoints; ++k )
        {
            // evaluate shape function at point k
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

    SECTION( "QUAD8: test d2NdXi2" )
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

    // tidy up
    delete tFunction;

    //------------------------------------------------------------------------------
}

