/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Bspline.cpp
 *
 */

#include <catch.hpp>
#include "cl_Bspline.hpp" // MOD/src

TEST_CASE("moris::model::Bspline",
        "[moris],[model],[Bspline]")
{
//    SECTION( "Building 1D b-spline basis functions and its derivative")
//    {
//        // Creates with the natural coordinate, a knot vector and a polynomial degree the basis functions and its derivative in 1D
//        moris::real Xi1 = 0.5; // Natural coordinate
//        moris::Mat<moris::real> KnotVector = {{-2, -1, 0, 1, 2, 3}}; // Knot vector for 1D
//        moris::uint PolynomialDegree1 = 2; // Polynomial degree in 1D
//        moris::uint n = KnotVector.n_cols()-PolynomialDegree1-1;   // Number of basis functions
//        moris::real comp_value = 0; // Variable for comparison purposes
//
//        // Creates 1D b-spline basis functions
//        moris::Mat<moris::real> BSPLINEd = moris::Bspline::build_spline_1d(Xi1,KnotVector,PolynomialDegree1);
//
//        for (moris::uint  i = 0; i<n;i++)
//        {
//            comp_value = comp_value + BSPLINEd(i);
//        }
//        // The sum of the shape functions must be one
//        REQUIRE( moris::equal_to( comp_value, 1.0 ) );
//
//        // Creates the derivative of 1D b-spline basis functions
//        moris::Mat<moris::real> bspline_deriv_1d = moris::Bspline::build_spline_deriv_1d(Xi1,KnotVector,PolynomialDegree1);
//        comp_value = 0;
//        for (moris::uint  i = 0; i<n;i++)
//        {
//            comp_value = comp_value + bspline_deriv_1d(i);
//        }
//        // The sum of the derivative of the shape functions must be zero
//        REQUIRE( moris::equal_to( comp_value, 0 ) );
//    }
//
//    SECTION( "Building 2D b-spline basis functions and its derivatives")
//    {
//        moris::Mat<moris::real> Xi = {{0, 1, 0}};// Natural coordinates in the 2 different directions, the third direction must be zero
//        moris::Mat<moris::real> KnotVector1 = {{-2, -1, 0, 1, 2, 3}}; // Knot vector for the direction Xi_1
//        moris::Mat<moris::real> KnotVector2 = {{-2, -1, 0, 1, 2, 3}}; // Knot vector for the direction Xi_2
//        moris::Mat<moris::real> KnotVector3 = {{0}}; // Knot vector for the third direction is zero
//        moris::Mat<moris::uint> PolynomialDegree = {{2, 2, 0}}; // Polynomial degree in each direction, the third direction must be zero
//        moris::Mat<moris::uint> Element_flag = {{0, 0, 0}}; // Gives the first basis functions, which starts in a specific element. A flag for each direction.
//        // In this example, the first basis function starts in this element
//        moris::uint dim = 2; // Dimensions of the problem
//        moris::real comp_value = 0; // Variable for comparison purposes
//        moris::real comp_value2 = 0; // Variable for comparison purposes
//
//        // Creates 2D b-spline basis functions with respect to the knot vectors and polynomial degrees
//        moris::Mat<moris::real> bspline_2d = moris::Bspline::build_spline_nd(Xi,KnotVector1,KnotVector2,KnotVector3,PolynomialDegree,Element_flag,dim);
//
//        for (moris::uint  i = 0; i<bspline_2d.n_cols();i++)
//        {
//            comp_value = comp_value + bspline_2d(i);
//        }
//        // The sum of the shape functions must be one
//        REQUIRE( moris::equal_to( comp_value, 1.0 ) );
//
//        //Creates the derivative of the 2D b-spline basis functions with respect to the knot vectors and polynomial degrees
//        moris::Mat<moris::real> bspline_deriv_2d = moris::Bspline::build_spline_deriv_nd(Xi,KnotVector1,KnotVector2,KnotVector3,PolynomialDegree,Element_flag,dim);
//
//        comp_value = 0;
//        for (moris::uint  i = 0; i<bspline_deriv_2d.n_cols();i++)
//        {
//            comp_value  = comp_value  + bspline_deriv_2d(0,i);
//            comp_value2 = comp_value2 + bspline_deriv_2d(1,i);
//        }
//        // The sum of the derivative of the shape functions must be zero in each direction
//        REQUIRE( moris::equal_to( comp_value, 0.0 ) );
//        REQUIRE( moris::equal_to( comp_value2, 0.0 ) );
//    }
//
//    SECTION( "Building uniform 2D b-spline basis functions and its derivatives")
//    {
//        moris::Mat<moris::real> Xi = {{0, 1, 0}};// Natural coordinates in the 2 different directions, the third direction must be zero
//        moris::Mat<moris::real> KnotVector = {{-2, -1, 0, 1, 2, 3}}; // Knot vector for each direction
//        moris::uint PolynomialDegree = 2; // Polynomial degree in each direction
//        moris::Mat<moris::uint> Element_flag = {{0, 0, 0}}; // Gives the first basis functions, which starts in a specific element. A flag for each direction.
//        // In this example, the first basis function starts in this element
//        moris::uint dim = 2; // Dimensions of the problem
//        moris::real comp_value = 0; // Variable for comparison purposes
//        moris::real comp_value2 = 0; // Variable for comparison purposes
//
//        // Creates uniform 2D b-spline basis functions with respect to the knot vectors and polynomial degrees
//        moris::Mat<moris::real> bspline_uniform_2d = moris::Bspline::build_spline_uniform_nd(Xi,KnotVector,PolynomialDegree,Element_flag,dim);
//
//        for (moris::uint  i = 0; i<bspline_uniform_2d.n_cols();i++)
//        {
//            comp_value = comp_value + bspline_uniform_2d(i);
//        }
//        // The sum of the shape functions must be one
//        REQUIRE( moris::equal_to( comp_value, 1.0 ) );
//
//        //Creates the derivative of the uniform 2D b-spline basis functions with respect to the knot vectors and polynomial degrees
//        moris::Mat<moris::real> bspline_deriv_uniform_2d = moris::Bspline::build_spline_deriv_uniform_nd(Xi,KnotVector,PolynomialDegree,Element_flag,dim);
//
//        comp_value = 0;
//        for (moris::uint  i = 0; i<bspline_deriv_uniform_2d.n_cols();i++)
//        {
//            comp_value  = comp_value  + bspline_deriv_uniform_2d(0,i);
//            comp_value2 = comp_value2 + bspline_deriv_uniform_2d(1,i);
//        }
//        // The sum of the derivative of the shape functions must be zero in each direction
//        REQUIRE( moris::equal_to( comp_value, 0.0 ) );
//        REQUIRE( moris::equal_to( comp_value2, 0.0 ) );
//    }
//
//
//    SECTION( "Building 3D b-spline basis functions and its derivatives")
//    {
//        moris::Mat<moris::real> Xi = {{0, 1, 0}};// Natural coordinates in the 2 different directions, the third direction must be zero
//        moris::Mat<moris::real> KnotVector1 = {{-2, -1, 0, 1, 2, 3}}; // Knot vector for the direction Xi_1
//        moris::Mat<moris::real> KnotVector2 = {{-2, -1, 0, 1, 2, 3}}; // Knot vector for the direction Xi_2
//        moris::Mat<moris::real> KnotVector3 = {{-2, -1, 0, 1, 2, 3}}; // Knot vector for the third direction is zero
//        moris::Mat<moris::uint> PolynomialDegree = {{2, 2, 2}}; // Polynomial degree in each direction, the third direction must be zero
//        moris::Mat<moris::uint> Element_flag = {{0, 0, 0}}; // Gives the first basis functions, which starts in a specific element. A flag for each direction.
//        // In this example, the first basis function starts in this element
//        moris::uint dim = 3; // Dimensions of the problem
//        moris::real comp_value = 0; // Variable for comparison purposes
//        moris::real comp_value2 = 0; // Variable for comparison purposes
//        moris::real comp_value3 = 0; // Variable for comparison purposes
//
//        // Creates 2D b-spline basis functions with respect to the knot vectors and polynomial degrees
//        moris::Mat<moris::real> bspline_3d = moris::Bspline::build_spline_nd(Xi,KnotVector1,KnotVector2,KnotVector3,PolynomialDegree,Element_flag,dim);
//
//        for (moris::uint  i = 0; i<bspline_3d.n_cols();i++)
//        {
//            comp_value = comp_value + bspline_3d(i);
//        }
//        // The sum of the shape functions must be one
//        REQUIRE( moris::equal_to( comp_value, 1.0 ) );
//
//        //Creates the derivative of the 2D b-spline basis functions with respect to the knot vectors and polynomial degrees
//        moris::Mat<moris::real> bspline_deriv_3d = moris::Bspline::build_spline_deriv_nd(Xi,KnotVector1,KnotVector2,KnotVector3,PolynomialDegree,Element_flag,dim);
//
//        comp_value = 0;
//        for (moris::uint  i = 0; i<bspline_deriv_3d.n_cols();i++)
//        {
//            comp_value  = comp_value   + bspline_deriv_3d(0,i);
//            comp_value2 = comp_value2 + bspline_deriv_3d(1,i);
//            comp_value3 = comp_value3 + bspline_deriv_3d(2,i);
//        }
//        // The sum of the derivative of the shape functions must be zero in each direction
//        REQUIRE( moris::equal_to( comp_value, 0.0 ) );
//        REQUIRE( moris::equal_to( comp_value2, 0.0 ) );
//        REQUIRE( moris::equal_to( comp_value3, 0.0 ) );
//    }
}

