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
        Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::QUADRATIC,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

        //set the coefficients xHat, tHat
        tGeomInterpolator.set_coeff( tXHat, tTHat );

        // specify test point in element: xi = (0.2, -0.6)
        Matrix< DDRMat > tXi(2,1);
        tXi( 0, 0 ) =   0.2;
        tXi( 1, 0 ) = - 0.6;

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

		Matrix< DDRMat > tdNdXi   = tGeomInterpolator.dNdXi( tXi );
		Matrix< DDRMat > td2NdXi2 = tGeomInterpolator.d2NdXi2( tXi );
		Matrix< DDRMat > td3NdXi3 = tGeomInterpolator.d3NdXi3( tXi );

		Matrix< DDRMat > tJ3a;
		Matrix< DDRMat > tJ3b;
		Matrix< DDRMat > tJ3c;
		Matrix< DDRMat > tJ = tGeomInterpolator.space_jacobian( tdNdXi );
		Matrix< DDRMat > tJ2b = tGeomInterpolator.second_space_jacobian( td2NdXi2 );

		tGeomInterpolator.space_jacobian_and_matrices_for_third_derivatives(tJ,
                                                                            tJ2b,
                                                                            tJ3a,
					                                                        tJ3b,
					                                                        tJ3c,
                                                                            tdNdXi,
                                                                            td2NdXi2,
				                                                            td3NdXi3);

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


//------------------------------------------------------------------------------
}


