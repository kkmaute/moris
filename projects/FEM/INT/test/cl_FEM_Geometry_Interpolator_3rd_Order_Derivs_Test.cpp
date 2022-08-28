/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Geometry_Interpolator_3rd_Order_Derivs_Test.cpp
 *
 */

#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/sr

using namespace moris;
using namespace fem;

TEST_CASE( "Geometry_Interpolator_Derivatives", "[moris],[fem],[GeoInterpolator_Derivatives]" )
{
    // define an epsilon environment
    double tEpsilon = 1E-12;

    SECTION( "Geometry Interpolator : 2D space - 3rd derivatives" )
    {
        // space and time geometry interpolator
        //------------------------------------------------------------------------------
        // create a distorted QUAD9 space element
        Matrix< DDRMat > tXHat( 9, 2 );

        tXHat( 0, 0 ) =  1.0;   tXHat( 0, 1 ) =  1.0;
        tXHat( 1, 0 ) = 11.0;   tXHat( 1, 1 ) =  2.0;
        tXHat( 2, 0 ) = 11.0;   tXHat( 2, 1 ) = 10.0;
        tXHat( 3, 0 ) =  1.0;   tXHat( 3, 1 ) =  9.0;
        tXHat( 4, 0 ) =  5.0;   tXHat( 4, 1 ) =  0.0;
        tXHat( 5, 0 ) = 12.0;   tXHat( 5, 1 ) =  7.0;
        tXHat( 6, 0 ) =  6.0;   tXHat( 6, 1 ) = 10.0;
        tXHat( 7, 0 ) =  0.0;   tXHat( 7, 1 ) =  5.0;
        tXHat( 8, 0 ) =  6.0;   tXHat( 8, 1 ) =  6.0;

        //create a line time element
        Matrix< DDRMat > tTHat( 2, 1 );
        tTHat( 0 ) = 0.0;
        tTHat( 1 ) = 1.0;

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tGeomInterpRule(
                mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeomInterpolator( tGeomInterpRule );

        //set nodal points xHat, tHat
        tGeomInterpolator.set_coeff( tXHat, tTHat );

        // set random test point in element: xi = (0.2, -0.6)
        Matrix< DDRMat > tXi(2,1);
        tXi( 0, 0 ) =   0.2;
        tXi( 1, 0 ) = - 0.6;

        tGeomInterpolator.set_space( tXi );

        // set nominal values --------------------------------------------------------------

        // first help matrix
        Matrix< DDRMat > tJ3a_nominal = {
                { 198.359290368,   1.450571968, 115.505512704,  22.419794304 },
                {   2.176782336, 244.844425216,  31.523033088, 152.166739968 },
                {  44.079842304,   8.016588544, 229.892401152,  84.262747392 },
                {   9.795520512,  44.303690752,  96.470424576, 246.606114816 }
        };

        // second help matrix
        Matrix< DDRMat > tJ3b_nominal = {
                {  16.79616,   5.29776,  30.55392 },
                { - 5.28768, -40.53888, -33.92256 },
                {  10.10880,   9.66880,   9.28160 },
                { - 5.96160, - 2.94560, - 4.73120 }
        };

        // third help matrix
        Matrix< DDRMat > tJ3c_nominal = {
                {   0.0,   0.0 },
                {   0.0,   0.0 },
                { - 2.2, - 3.2 },
                { - 1.6, - 0.6 }
        };

        // ---------------------------------------------------------------------------------

        // construct jacobian matrices

        Matrix< DDRMat > tdNdXi   = tGeomInterpolator.dNdXi();
        Matrix< DDRMat > td2NdXi2 = tGeomInterpolator.d2NdXi2();
        Matrix< DDRMat > td3NdXi3 = tGeomInterpolator.d3NdXi3();

        Matrix< DDRMat > tJ3a;
        Matrix< DDRMat > tJ3b;
        Matrix< DDRMat > tJ3c;
        Matrix< DDRMat > tJ = tGeomInterpolator.space_jacobian();
        Matrix< DDRMat > tJ2b;
        tGeomInterpolator.second_space_jacobian( tJ2b );

        tGeomInterpolator.space_jacobian_and_matrices_for_third_derivatives(tJ,
                tJ2b,
                tJ3a, tJ3b, tJ3c,
                tdNdXi, td2NdXi2, td3NdXi3 );

        // check if all entries in the matrices are equal to the expected values
        bool tJ3aCheck = true;
        for ( uint i = 0; i < 4; i++)
        {
            for ( uint j = 0; j < 4; j++)
            {
                tJ3aCheck = tJ3aCheck && ( std::abs( tJ3a_nominal( i, j ) - tJ3a( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tJ3aCheck );

        bool tJ3bCheck = true;
        for ( uint i = 0; i < 4; i++)
        {
            for ( uint j = 0; j < 3; j++)
            {
                tJ3bCheck = tJ3bCheck && ( std::abs( tJ3b_nominal( i, j ) - tJ3b( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tJ3bCheck );

        bool tJ3cCheck = true;
        for ( uint i = 0; i < 4; i++)
        {
            for ( uint j = 0; j < 2; j++)
            {
                tJ3cCheck = tJ3cCheck && ( std::abs( tJ3c_nominal( i, j ) - tJ3c( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tJ3cCheck );
    }

    //___________________________________________________________________________________________________________
    SECTION( "Geometry Interpolator : 3D space - 3rd derivatives" )
    {
        // space and time geometry interpolator
        //------------------------------------------------------------------------------
        // create a distorted HEX27 space element
        Matrix< DDRMat > tXHat = {
                { 1.0,  1.0,  4.0},
                {11.0,  2.0,  4.0},
                {11.0, 10.0,  4.0},
                { 1.0,  9.0,  4.0},
                { 1.0,  2.0,  0.0},
                {11.0,  1.0,  0.0},
                {11.0,  9.0,  0.0},
                { 1.0, 10.0,  0.0},
                { 5.0,  0.0,  5.0},
                {12.0,  7.0,  6.0},
                { 6.0, 10.0,  5.0},
                { 0.0,  5.0,  6.0},
                { 1.0,  2.0,  2.0},
                {11.0,  1.0,  2.0},
                {12.0,  9.0,  2.0},
                { 2.0, 10.0,  2.0},
                { 5.0,  0.0, -1.0},
                {12.0,  7.0, -1.0},
                { 6.0, 10.0, -1.0},
                { 0.0,  5.0, -1.0},
                { 6.0,  6.0,  3.0},
                { 6.0,  6.0,  7.0},
                { 6.0,  6.0, -2.0},
                { 0.0,  6.0,  2.0},
                {12.0,  6.0,  2.0},
                { 6.0,  0.0,  2.0},
                { 5.0, 10.0,  2.0}
        };

        //create a line time element
        Matrix< DDRMat > tTHat( 2, 1 );
        tTHat( 0 ) = 0.0;
        tTHat( 1 ) = 1.0;

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeomInterpolator( tGeomInterpRule );

        //set nodal points xHat, tHat
        tGeomInterpolator.set_coeff( tXHat, tTHat );

        // set random test point in element: xi = (0.2, -0.6, 0.3)
        Matrix< DDRMat > tXi(3,1);
        tXi( 0, 0 ) =   0.2;
        tXi( 1, 0 ) = - 0.6;
        tXi( 2, 0 ) =   0.3;

        tGeomInterpolator.set_space( tXi );

        // set nominal values --------------------------------------------------------------

        // first help matrix
        Matrix< DDRMat > tJ3a_nominal = {
                {172.80124713094300,   0.00349796383200,  -0.00144136526234,  14.12833361043460,  -10.51341610431280,   0.38504701002240,   -0.00780889317120,   0.21321589594522,   0.00581088651264,   -0.57305547104256},
                {  0.07343372746938, 250.38054419795200,   0.17780771091917,   3.31577584411584,    0.29582520033370,  49.90608460763520,   67.01497157456640,   0.39724048775347,   5.97890759913216,    8.90499127262976},
                { -0.04025486249165,   0.00024297062400, -68.75169402414690,   0.02198795010048,   -1.44353711377613,  -0.00400340828160,   -0.04785399889920, -17.25505334319510,   3.14167637901312,    0.52565777252352},
                { 12.99161392776810,   0.14523693307200,   0.00717518160691, 196.24625746965500,   16.91847133677160,  10.66785905037600,   -0.20319468312000,  -0.70225556794163,   0.06113891501056,   -6.99463078116864},
                {-10.63254829812940,   0.00143789817600,  -0.05227296908902,   1.35634913206272, -126.66293155646700,   0.09762275001600,   -0.09653984524800,   5.15065976694375,   0.14128908115968,   -6.99429973635072},
                {  0.97674082363699,   6.03029869411200,  -0.03571837925990,  29.42867362276320,    2.60336881755699, 222.06798510207400,   -3.41135227982080,   1.70803077999206,  -0.75270378268800,   38.95692157907840},
                { -0.06009935781427,   2.47885796121600,  -1.29537306313114,  -1.79818082020608,   -0.87979197409587, -13.28525845463040, -162.29788820295700,  -2.03769936427622, -29.01881399104510,  -24.02496221545470},
                {  0.65422608452813,   0.00059107276800,  -1.89574660135322,  -0.22040361584640,   15.62707255069900,   0.01519522467840,   -0.07804920913920,  93.15976386222490,   2.60532932567040,   -2.41659503861760},
                {  0.04918629264998,   0.02454159052800,   9.43711250772787,   0.72239699460096,    1.24192688526950,  -0.26794938193920,   -3.22018422435840,   8.60681773034701, 105.48899230941200,   17.46012144451580},
                { -0.79938058810982,   0.05970213849600,   0.26021720943821, -11.90780023368960,  -10.61246154501730,   1.86668424103680,   -3.95862203005440, -12.61544026574030,   2.56300784215040, -143.67019439769600}
        };

        // second help matrix
        Matrix< DDRMat > tJ3b_nominal = {
                {-5.855099904,   0.710424000,  0.191399424,  -0.785862720,  -9.318928896,  25.907653440},
                {-1.662812208, -42.316999200, -1.581018912, -21.496368960,  -3.410245872, -27.838712880},
                { 1.174569984,   0.061401600,  7.236304896,  -4.141317120,  14.645385216,  -0.551093760},
                {11.968956480,   9.891266400, -0.218916480,  -2.859043200,  -5.545145280,   0.936088800},
                { 2.045053440,   0.192067200,  2.188712960,  -6.328416000,   7.757007360,   2.971584000},
                {-6.461459360,   2.105758000, -0.385376640,  -5.177408400,  -4.213285440,   1.205926000},
                { 0.938045760,   2.280624000,  1.442538240, -17.539324800,   4.609079040,   8.136297600},
                {-6.481520640,   0.088728000, -4.469288960,  -2.613644800,  -4.945781760,   1.461225600},
                {-0.874849920,   2.091280000, 17.130424320,  -5.365731200,  -4.167143680,  -7.122385600},
                { 2.921261184,   2.007724800,  2.341453056,   2.497464960, -15.912760704,   2.378381280},

        };

        // third help matrix
        Matrix< DDRMat > tJ3c_nominal = {
                { 0.000,  0.000,  0.000},
                { 0.000,  0.000,  0.000},
                { 0.000,  0.000,  0.000},
                {-0.562, -3.200, -2.184},
                { 0.864,  0.000,  2.768},
                {-1.236, -0.990,  0.728},
                {-0.024, -0.320,  3.552},
                { 0.576,  1.640,  0.512},
                { 1.928,  0.240, -1.104},
                {-0.216,  0.960,  0.288},
        };

        // third shape function derivatives
        Matrix< DDRMat > td3NdXi3_nominal = {
                { 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.000},
                { 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.000},
                { 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.000},
                { 0.1155,  0.1155,  0.0105,  0.0105, -0.2145, -0.2145, -0.0195, -0.0195, -0.2310, -0.1260, -0.0210, -0.1260, -1.0010, -1.0010, -0.0910, -0.0910,  0.4290,  0.2340,  0.0390,  0.2340, -2.1840,  0.2520, -0.4680,  1.0920,  1.0920,  2.0020,  0.182},
                {-0.0960, -0.0960,  0.0240,  0.0240,  0.3840,  0.3840, -0.0960, -0.0960,  0.1920, -0.1280, -0.0480, -0.1280, -0.2880, -0.2880,  0.0720,  0.0720, -0.7680,  0.5120,  0.1920,  0.5120,  0.7680,  0.2560, -1.0240, -0.3840, -0.3840,  0.5760, -0.144},
                { 0.0315, -0.0735, -0.0735,  0.0315, -0.0585,  0.1365,  0.1365, -0.0585,  0.0420,  0.1470,  0.0420, -0.0630, -0.2730,  0.6370,  0.6370, -0.2730, -0.0780, -0.2730, -0.0780,  0.1170,  0.7280, -0.0840,  0.1560,  0.5460, -1.2740, -0.3640, -0.364},
                { 0.0160, -0.0240, -0.0240,  0.0160, -0.0640,  0.0960,  0.0960, -0.0640, -0.1920,  0.0480, -0.1920, -0.0320,  0.0480, -0.0720, -0.0720,  0.0480,  0.7680, -0.1920,  0.7680,  0.1280,  1.1520,  0.3840, -1.5360, -0.0960,  0.1440, -0.5760, -0.576},
                {-0.1440,  0.3360, -0.0840,  0.0360, -0.1440,  0.3360, -0.0840,  0.0360, -0.1920,  0.4480,  0.0480, -0.1920,  0.2880, -0.6720,  0.1680, -0.0720, -0.1920,  0.4480,  0.0480, -0.1920,  0.5120, -0.2560, -0.2560,  0.3840, -0.8960,  0.3840, -0.096},
                { 0.0880, -0.1320, -0.0120,  0.0080,  0.0880, -0.1320, -0.0120,  0.0080, -1.0560,  0.1440, -0.0960, -0.0960, -0.1760,  0.2640,  0.0240, -0.0160, -1.0560,  0.1440, -0.0960, -0.0960, -2.3040,  1.1520,  1.1520,  0.1920, -0.2880,  2.1120,  0.192},
                {-0.0660,  0.1540,  0.0140, -0.0060,  0.2640, -0.6160, -0.0560,  0.0240, -0.0880, -0.1680, -0.0080,  0.0720, -0.1980,  0.4620,  0.0420, -0.0180,  0.3520,  0.6720,  0.0320, -0.2880,  0.2880,  0.0960, -0.3840,  0.2160, -0.5040, -0.2640, -0.024}
        };

        // ---------------------------------------------------------------------------------

        // construct jacobian matrices

        Matrix< DDRMat > tdNdXi   = tGeomInterpolator.dNdXi();
        Matrix< DDRMat > td2NdXi2 = tGeomInterpolator.d2NdXi2();
        Matrix< DDRMat > td3NdXi3 = tGeomInterpolator.d3NdXi3();

        Matrix< DDRMat > tJ3a;
        Matrix< DDRMat > tJ3b;
        Matrix< DDRMat > tJ3c;
        Matrix< DDRMat > tJ = tGeomInterpolator.space_jacobian();
        Matrix< DDRMat > tJ2b;
        tGeomInterpolator.second_space_jacobian( tJ2b );

        tGeomInterpolator.space_jacobian_and_matrices_for_third_derivatives(tJ,
                tJ2b,
                tJ3a, tJ3b, tJ3c,
                tdNdXi, td2NdXi2, td3NdXi3 );

        // for debugging
        //        Matrix< DDRMat > t_Delta_d3NdXi3 = td3NdXi3_nominal - td3NdXi3;
        //        print( t_Delta_d3NdXi3 , "delta_d3NdXi3" );
        //        Matrix< DDRMat > t_Delta_J3a = tJ3a_nominal - tJ3a;
        //        print( t_Delta_J3a , "delta_J3a" );
        //        Matrix< DDRMat > t_Delta_J3b = tJ3b_nominal - tJ3b;
        //        print( t_Delta_J3b , "delta_J3b" );
        //        Matrix< DDRMat > t_Delta_J3c = tJ3c_nominal - tJ3c;
        //        print( t_Delta_J3c , "delta_J3c" );

        // check if all entries in the matrices are equal to the expected values
        bool tJ3aCheck = true;
        for ( uint i = 0; i < 10; i++)
        {
            for ( uint j = 0; j < 10; j++)
            {
                tJ3aCheck = tJ3aCheck && ( std::abs( tJ3a_nominal( i, j ) - tJ3a( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tJ3aCheck );

        bool tJ3bCheck = true;
        for ( uint i = 0; i < 10; i++)
        {
            for ( uint j = 0; j < 6; j++)
            {
                tJ3bCheck = tJ3bCheck && ( std::abs( tJ3b_nominal( i, j ) - tJ3b( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tJ3bCheck );

        bool tJ3cCheck = true;
        for ( uint i = 0; i < 10; i++)
        {
            for ( uint j = 0; j < 3; j++)
            {
                tJ3cCheck = tJ3cCheck && ( std::abs( tJ3c_nominal( i, j ) - tJ3c( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tJ3cCheck );

        bool td3NdXi3Check = true;
        for ( uint i = 0; i < 10; i++)
        {
            for ( uint j = 0; j < 27; j++)
            {
                td3NdXi3Check = td3NdXi3Check && ( std::abs( td3NdXi3_nominal( i, j ) - td3NdXi3( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( td3NdXi3Check );

    }

    //------------------------------------------------------------------------------
}
