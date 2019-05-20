/*
 * cl_FEM_Interpolation_Function_Lagrange_Hex27.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX27_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX27_HPP_

#include "assert.h"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 27 >::get_number_of_param_dimensions() const
        {
            return 3;
        }

//------------------------------------------------------------------------------

    template<>
    mtk::Interpolation_Order
    Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 27 >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::QUADRATIC;
    }

//------------------------------------------------------------------------------

    template<>
    Matrix< DDRMat >
    Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 27 >::get_param_coords() const
    {
        Matrix< DDRMat > tXiHat( 3, 27 );
        tXiHat( 0,  0 ) = -1.000000; tXiHat( 1,  0 ) = -1.000000; tXiHat( 2,  0 ) = -1.000000;
        tXiHat( 0,  1 ) =  1.000000; tXiHat( 1,  1 ) = -1.000000; tXiHat( 2,  1 ) = -1.000000;
        tXiHat( 0,  2 ) =  1.000000; tXiHat( 1,  2 ) =  1.000000; tXiHat( 2,  2 ) = -1.000000;
        tXiHat( 0,  3 ) = -1.000000; tXiHat( 1,  3 ) =  1.000000; tXiHat( 2,  3 ) = -1.000000;
        tXiHat( 0,  4 ) = -1.000000; tXiHat( 1,  4 ) = -1.000000; tXiHat( 2,  4 ) =  1.000000;
        tXiHat( 0,  5 ) =  1.000000; tXiHat( 1,  5 ) = -1.000000; tXiHat( 2,  5 ) =  1.000000;
        tXiHat( 0,  6 ) =  1.000000; tXiHat( 1,  6 ) =  1.000000; tXiHat( 2,  6 ) =  1.000000;
        tXiHat( 0,  7 ) = -1.000000; tXiHat( 1,  7 ) =  1.000000; tXiHat( 2,  7 ) =  1.000000;

        tXiHat( 0,  8 ) =  0.000000; tXiHat( 1,  8 ) = -1.000000; tXiHat( 2,  8 ) = -1.000000;
        tXiHat( 0,  9 ) =  1.000000; tXiHat( 1,  9 ) =  0.000000; tXiHat( 2,  9 ) = -1.000000;
        tXiHat( 0, 10 ) =  0.000000; tXiHat( 1, 10 ) =  1.000000; tXiHat( 2, 10 ) = -1.000000;
        tXiHat( 0, 11 ) = -1.000000; tXiHat( 1, 11 ) =  0.000000; tXiHat( 2, 11 ) = -1.000000;

        tXiHat( 0, 12 ) = -1.000000; tXiHat( 1, 12 ) = -1.000000; tXiHat( 2, 12 ) =  0.000000;
        tXiHat( 0, 13 ) =  1.000000; tXiHat( 1, 13 ) = -1.000000; tXiHat( 2, 13 ) =  0.000000;
        tXiHat( 0, 14 ) =  1.000000; tXiHat( 1, 14 ) =  1.000000; tXiHat( 2, 14 ) =  0.000000;
        tXiHat( 0, 15 ) = -1.000000; tXiHat( 1, 15 ) =  1.000000; tXiHat( 2, 15 ) =  0.000000;

        tXiHat( 0, 16 ) =  0.000000; tXiHat( 1, 16 ) = -1.000000; tXiHat( 2, 16 ) =  1.000000;
        tXiHat( 0, 17 ) =  1.000000; tXiHat( 1, 17 ) =  0.000000; tXiHat( 2, 17 ) =  1.000000;
        tXiHat( 0, 18 ) =  0.000000; tXiHat( 1, 18 ) =  1.000000; tXiHat( 2, 18 ) =  1.000000;
        tXiHat( 0, 19 ) = -1.000000; tXiHat( 1, 19 ) =  0.000000; tXiHat( 2, 19 ) =  1.000000;

        tXiHat( 0, 20 ) =  0.000000; tXiHat( 1, 20 ) =  0.000000; tXiHat( 2, 20 ) =  0.000000;
        tXiHat( 0, 21 ) =  0.000000; tXiHat( 1, 21 ) =  0.000000; tXiHat( 2, 21 ) = -1.000000;
        tXiHat( 0, 22 ) =  0.000000; tXiHat( 1, 22 ) =  0.000000; tXiHat( 2, 22 ) =  1.000000;
        tXiHat( 0, 23 ) = -1.000000; tXiHat( 1, 23 ) =  0.000000; tXiHat( 2, 23 ) =  0.000000;
        tXiHat( 0, 24 ) =  1.000000; tXiHat( 1, 24 ) =  0.000000; tXiHat( 2, 24 ) =  0.000000;
        tXiHat( 0, 25 ) =  0.000000; tXiHat( 1, 25 ) = -1.000000; tXiHat( 2, 25 ) =  0.000000;
        tXiHat( 0, 26 ) =  0.000000; tXiHat( 1, 26 ) =  1.000000; tXiHat( 2, 26 ) =  0.000000;
        return tXiHat;
    }

//------------------------------------------------------------------------------

    template<>
    Matrix< DDRMat >
    Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 27  >::eval_N( const Matrix< DDRMat > & aXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 3,
                "HEX27 - eval_N: aXi not allocated or hat wrong size." );

        // unpack xi and eta from input vector
        auto    xi = aXi( 0 );
        auto   eta = aXi( 1 );
        auto  zeta = aXi( 2 );

        // often used constants
        auto    xi2 = std::pow(   xi, 2 );
        auto   eta2 = std::pow(  eta, 2 );
        auto  zeta2 = std::pow( zeta, 2 );

        auto a = -0.25 * eta * zeta;
        auto b = -0.25 * xi * zeta;
        auto c = -0.25 * xi * eta;
        auto d = 0.125 * xi * eta * zeta;

        // populate output matrix
        Matrix< DDRMat > tN(1,27);
        tN(  0 ) = d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        tN(  1 ) = d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        tN(  2 ) = d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        tN(  3 ) = d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        tN(  4 ) = d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        tN(  5 ) = d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        tN(  6 ) = d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        tN(  7 ) = d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        tN(  8 ) = a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
        tN(  9 ) = b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        tN( 10 ) = a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
        tN( 11 ) = b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        tN( 12 ) = c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );
        tN( 13 ) = c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );
        tN( 14 ) = c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );
        tN( 15 ) = c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );
        tN( 16 ) = a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
        tN( 17 ) = b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        tN( 18 ) = a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
        tN( 19 ) = b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        tN( 20 ) = -( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 );
        tN( 21 ) = ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
        tN( 22 ) = ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
        tN( 23 ) = ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.5;
        tN( 24 ) = ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.5;
        tN( 25 ) = ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.5;
        tN( 26 ) = ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.5;
        return tN;
    }

//------------------------------------------------------------------------------

    template<>
    Matrix< DDRMat >
    Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 27  >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 3,
                "HEX27 - eval_dNdXi: aXi not allocated or hat wrong size." );

        // unpack xi and eta from input vector
        auto   xi = aXi( 0 );
        auto  eta = aXi( 1 );
        auto zeta = aXi( 2 );

        // often used constants
        auto   xi2 = std::pow(   xi, 2 );
        auto  eta2 = std::pow(  eta, 2 );
        auto zeta2 = std::pow( zeta, 2 );

        auto a = 0.125*eta*zeta;
        auto b = 0.125*xi*zeta;
        auto c = 0.125*xi*eta;
        auto d = -0.5*xi*eta*zeta;

        Matrix< DDRMat > tdNdXi(3,27);
        tdNdXi( 0,  0 ) = a * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 1,  0 ) = b * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 2,  0 ) = c * ( 2.0 * zeta - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );

        tdNdXi( 0,  1 ) = a * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 1,  1 ) = b * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 2,  1 ) = c * ( 2.0 * zeta - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );

        tdNdXi( 0,  2 ) = a * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 1,  2 ) = b * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 2,  2 ) = c * ( 2.0 * zeta - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );

        tdNdXi( 0,  3 ) = a * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 1,  3 ) = b * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 2,  3 ) = c * ( 2.0 * zeta - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );

        tdNdXi( 0,  4 ) = a * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 1,  4 ) = b * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 2,  4 ) = c * ( 2.0 * zeta + 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );

        tdNdXi( 0,  5 ) = a * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 1,  5 ) = b * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 2,  5 ) = c * ( 2.0 * zeta + 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );

        tdNdXi( 0,  6 ) = a * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 1,  6 ) = b * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 2,  6 ) = c * ( 2.0 * zeta + 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );

        tdNdXi( 0,  7 ) = a * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 1,  7 ) = b * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 2,  7 ) = c * ( 2.0 * zeta + 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );

        tdNdXi( 0,  8 ) = d * ( eta - 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 1,  8 ) = - ( zeta * ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
        tdNdXi( 2,  8 ) = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( eta - 1.0 ) ) * 0.25;

        tdNdXi( 0,  9 ) = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
        tdNdXi( 1,  9 ) = d * ( xi + 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 2,  9 ) = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( xi + 1.0 ) ) * 0.25;


        tdNdXi( 0, 10 ) = d * ( eta + 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 1, 10 ) =  - ( zeta * ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
        tdNdXi( 2, 10 ) = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( eta + 1.0 ) ) * 0.25;

        tdNdXi( 0, 11 ) = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
        tdNdXi( 1, 11 ) = d * ( xi - 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 2, 11 ) = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( xi - 1.0 ) ) * 0.25;

        tdNdXi( 0, 12 ) = - ( eta * ( 2.0 * xi - 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
        tdNdXi( 1, 12 ) = - ( xi * ( 2.0 * eta - 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
        tdNdXi( 2, 12 ) = d * ( eta - 1.0 ) * ( xi - 1.0 );

        tdNdXi( 0, 13 ) = - ( eta * ( 2.0 * xi + 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
        tdNdXi( 1, 13 ) = - ( xi * ( 2.0 * eta - 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
        tdNdXi( 2, 13 ) = d * ( eta - 1.0 ) * ( xi + 1.0 );

        tdNdXi( 0, 14 ) = - ( eta * ( 2.0 * xi + 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
        tdNdXi( 1, 14 ) = - ( xi * ( 2.0 * eta + 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
        tdNdXi( 2, 14 ) = d * ( eta + 1.0 ) * ( xi + 1.0 );

        tdNdXi( 0, 15 ) = - ( eta * ( 2.0 * xi - 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
        tdNdXi( 1, 15 ) = - ( xi * ( 2.0 * eta + 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
        tdNdXi( 2, 15 ) = d * ( eta + 1.0 ) * ( xi - 1.0 );

        tdNdXi( 0, 16 ) = d * ( eta - 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 1, 16 ) = - ( zeta * ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
        tdNdXi( 2, 16 ) = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( eta - 1.0 ) ) * 0.25;


        tdNdXi( 0, 17 ) = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
        tdNdXi( 1, 17 ) = d * ( xi + 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 2, 17 ) = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( xi + 1.0 ) ) * 0.25;

        tdNdXi( 0, 18 ) = d * ( eta + 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 1, 18 ) = - ( zeta * ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
        tdNdXi( 2, 18 ) = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( eta + 1.0 ) ) * 0.25;

        tdNdXi( 0, 19 ) = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
        tdNdXi( 1, 19 ) = d * ( xi - 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 2, 19 ) = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( xi - 1.0 ) ) * 0.25;

        tdNdXi( 0, 20 ) = - 2.0 * xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 );
        tdNdXi( 1, 20 ) = - 2.0 * eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 );
        tdNdXi( 2, 20 ) = - 2.0 * zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 );

        tdNdXi( 0, 21 ) = 8.0 * b * ( eta2 - 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 1, 21 ) = 8.0 * a * ( xi2 - 1.0 ) * ( zeta - 1.0 );
        tdNdXi( 2, 21 ) = ( ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.5;

        tdNdXi( 0, 22 ) = 8.0 * b * ( eta2 - 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 1, 22 ) = 8.0 * a * ( xi2 - 1.0 ) * ( zeta + 1.0 );
        tdNdXi( 2, 22 ) = ( ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.5;

        tdNdXi( 0, 23 ) =( ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
        tdNdXi( 1, 23 ) = 8.0 * c * ( zeta2 - 1.0 ) * ( xi - 1.0 );
        tdNdXi( 2, 23 ) = 8.0 * b * ( eta2 - 1.0 ) * ( xi - 1.0 );

        tdNdXi( 0, 24 ) =( ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
        tdNdXi( 1, 24 ) = 8.0 * c * ( zeta2 - 1.0 ) * ( xi + 1.0 );
        tdNdXi( 2, 24 ) = 8.0 * b * ( eta2 - 1.0 ) * ( xi + 1.0 );

        tdNdXi( 0, 25 ) = 8.0 * c * ( zeta2 - 1.0 ) * ( eta - 1.0 );
        tdNdXi( 1, 25 ) = ( ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
        tdNdXi( 2, 25 ) = 8.0 * a * ( xi2 - 1.0 ) * ( eta - 1.0 );

        tdNdXi( 0, 26 ) = 8.0 * c * ( zeta2 - 1.0 ) * ( eta + 1.0 );
        tdNdXi( 1, 26 ) = ( ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
        tdNdXi( 2, 26 ) = 8.0 * a * ( xi2 - 1.0 ) * ( eta + 1.0 );
        return tdNdXi;

    }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 27  >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "HEX27 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // often used constants
            auto   xi2 = std::pow(   xi, 2 );
            auto  eta2 = std::pow(  eta, 2 );
            auto zeta2 = std::pow( zeta, 2 );

            auto a = eta * zeta;
            auto b = xi * zeta;
            auto c = xi * eta;
            auto d = 2 * xi* eta * zeta;

            Matrix< DDRMat > td2NdXi2(6,27);
            td2NdXi2( 0,  0 ) = ( a * (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 1,  0 ) = ( b * (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 2,  0 ) = ( c * (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            td2NdXi2( 3,  0 ) = (   xi * ( 2.0 *  eta - 1.0 ) * ( 2.0 * zeta - 1.0 ) * (   xi - 1.0 ) ) * 0.125;
            td2NdXi2( 4,  0 ) = (  eta * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta - 1.0 ) * (  eta - 1.0 ) ) * 0.125;
            td2NdXi2( 5,  0 ) = ( zeta * ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.125;

            td2NdXi2( 0,  1 ) = ( a * (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 1,  1 ) = ( b * (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 2,  1 ) = ( c * (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            td2NdXi2( 3,  1 ) = (   xi * ( 2.0 *  eta - 1.0 ) * ( 2.0 * zeta - 1.0 ) * (   xi + 1.0 ) ) * 0.125;
            td2NdXi2( 4,  1 ) = (  eta * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta - 1.0 ) * (  eta - 1.0 ) ) * 0.125;
            td2NdXi2( 5,  1 ) = ( zeta * ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.125;

            td2NdXi2( 0,  2 ) = ( a * (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 1,  2 ) = ( b * (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 2,  2 ) = ( c * (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            td2NdXi2( 3,  2 ) = (   xi * ( 2.0 *  eta + 1.0 ) * ( 2.0 * zeta - 1.0 ) * (   xi + 1.0 ) ) * 0.125;
            td2NdXi2( 4,  2 ) = (  eta * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta - 1.0 ) * (  eta + 1.0 ) ) * 0.125;
            td2NdXi2( 5,  2 ) = ( zeta * ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.125;

            td2NdXi2( 0,  3 ) = ( a * (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 1,  3 ) = ( b * (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 2,  3 ) = ( c * (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            td2NdXi2( 3,  3 ) = (   xi * ( 2.0 *  eta + 1.0 ) * ( 2.0 * zeta - 1.0 ) * (   xi - 1.0 ) ) * 0.125;
            td2NdXi2( 4,  3 ) = (  eta * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta - 1.0 ) * (  eta + 1.0 ) ) * 0.125;
            td2NdXi2( 5,  3 ) = ( zeta * ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.125;

            td2NdXi2( 0,  4 ) = ( a * (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 1,  4 ) = ( b * (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 2,  4 ) = ( c * (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            td2NdXi2( 3,  4 ) = (   xi * ( 2.0 *  eta - 1.0 ) * ( 2.0 * zeta + 1.0 ) * (   xi - 1.0 ) ) * 0.125;
            td2NdXi2( 4,  4 ) = (  eta * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta + 1.0 ) * (  eta - 1.0 ) ) * 0.125;
            td2NdXi2( 5,  4 ) = ( zeta * ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.125;

            td2NdXi2( 0,  5 ) = ( a * (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 1,  5 ) = ( b * (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 2,  5 ) = ( c * (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            td2NdXi2( 3,  5 ) = (   xi * ( 2.0 *  eta - 1.0 ) * ( 2.0 * zeta + 1.0 ) * (   xi + 1.0 ) ) * 0.125;
            td2NdXi2( 4,  5 ) = (  eta * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta + 1.0 ) * (  eta - 1.0 ) ) * 0.125;
            td2NdXi2( 5,  5 ) = ( zeta * ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.125;

            td2NdXi2( 0,  6 ) = ( a * (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 1,  6 ) = ( b * (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 2,  6 ) = ( c * (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            td2NdXi2( 3,  6 ) = (   xi * ( 2.0 *  eta + 1.0 ) * ( 2.0 * zeta + 1.0 ) * (   xi + 1.0 ) ) * 0.125;
            td2NdXi2( 4,  6 ) = (  eta * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta + 1.0 ) * (  eta + 1.0 ) ) * 0.125;
            td2NdXi2( 5,  6 ) = ( zeta * ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.125;

            td2NdXi2( 0,  7 ) = ( a * (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 1,  7 ) = ( b * (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 2,  7 ) = ( c * (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            td2NdXi2( 3,  7 ) = (   xi * ( 2.0 *  eta + 1.0 ) * ( 2.0 * zeta + 1.0 ) * (   xi - 1.0 ) ) * 0.125;
            td2NdXi2( 4,  7 ) = (  eta * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta + 1.0 ) * (  eta + 1.0 ) ) * 0.125;
            td2NdXi2( 5,  7 ) = ( zeta * ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.125;

            td2NdXi2( 0,  8 ) = - ( a * (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 1,  8 ) = - ( zeta * (   xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 2,  8 ) = - (  eta * (   xi2 - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 3,  8 ) = - ( ( 2.0 *  eta - 1.0 ) * (   xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 4,  8 ) = - ( c * ( 2.0 * zeta - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 5,  8 ) = - ( b * ( 2.0 *  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;

            td2NdXi2( 0,  9 ) = - ( zeta * (  eta2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 1,  9 ) = - ( b * (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 2,  9 ) = - (   xi * (  eta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 3,  9 ) = - ( c * ( 2.0 * zeta - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 4,  9 ) = - ( (  eta2 - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 5,  9 ) = - ( a * ( 2.0 *   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;

            td2NdXi2( 0, 10 ) = - ( a * (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 1, 10 ) = - ( zeta * (   xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 2, 10 ) = - (  eta * (   xi2 - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 3, 10 ) = - ( ( 2.0 *  eta + 1.0 ) * (   xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 4, 10 ) = - ( c * ( 2.0 * zeta - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 5, 10 ) = - ( b * ( 2.0 *  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;

            td2NdXi2( 0, 11 ) = - ( zeta * (  eta2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 1, 11 ) = - ( b * (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 2, 11 ) = - (   xi * (  eta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 3, 11 ) = - ( c * ( 2.0 * zeta - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 4, 11 ) = - ( (  eta2 - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 5, 11 ) = - ( a * ( 2.0 *   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;

            td2NdXi2( 0, 12 ) = - (  eta * ( zeta2 - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 1, 12 ) = - (   xi * ( zeta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 2, 12 ) = - ( c * (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 3, 12 ) = - ( b * ( 2.0 *  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 4, 12 ) = - ( a * ( 2.0 *   xi - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 5, 12 ) = - ( ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.25;

            td2NdXi2( 0, 13 ) = - (  eta * ( zeta2 - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 1, 13 ) = - (   xi * ( zeta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 13 ) = - ( c * (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 3, 13 ) = - ( b * ( 2.0 *  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 4, 13 ) = - ( a * ( 2.0 *   xi + 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 5, 13 ) = - ( ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta2 - 1.0 ) ) * 0.25;

            td2NdXi2( 0, 14 ) = - (  eta * ( zeta2 - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 14 ) = - (   xi * ( zeta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 14 ) = - ( c * (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 3, 14 ) = - ( b * ( 2.0 *  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 4, 14 ) = - ( a * ( 2.0 *   xi + 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 5, 14 ) = - ( ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta2 - 1.0 ) ) * 0.25;

            td2NdXi2( 0, 15 ) = - (  eta * ( zeta2 - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 15 ) = - (   xi * ( zeta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 2, 15 ) = - ( c * (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 3, 15 ) = - ( b * ( 2.0 *  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 4, 15 ) = - ( a * ( 2.0 *   xi - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 5, 15 ) = - ( ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.25;

            td2NdXi2( 0, 16 ) = - ( a * (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 16 ) = - ( zeta * (   xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 16 ) = - (  eta * (   xi2 - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 3, 16 ) = - ( ( 2.0 *  eta - 1.0 ) * (   xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 4, 16 ) =  - ( c * ( 2.0 * zeta + 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 5, 16 ) =  - ( b * ( 2.0 *  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;

            td2NdXi2( 0, 17 ) = - ( zeta * (  eta2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 17 ) = - ( b * (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 17 ) = - (   xi * (  eta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 3, 17 ) = - ( c * ( 2.0 * zeta + 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 4, 17 ) =  - ( (  eta2 - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 5, 17 ) =  - ( a * ( 2.0 *   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;

            td2NdXi2( 0, 18 ) = - ( a * (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 18 ) = - ( zeta * (   xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 18 ) = - (  eta * (   xi2 - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 3, 18 ) = - ( ( 2.0 *  eta + 1.0 ) * (   xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 4, 18 ) = - ( c * ( 2.0 * zeta + 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 5, 18 ) = - ( b * ( 2.0 *  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;

            td2NdXi2( 0, 19 ) = - ( zeta * (  eta2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 19 ) = - ( b * (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 19 ) = - (   xi * (  eta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 3, 19 ) = - ( c * ( 2.0 * zeta + 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 4, 19 ) =  - ( (  eta2 - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 5, 19 ) =  - ( a * ( 2.0 *   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;

            td2NdXi2( 0, 20 ) = - 2.0 * (  eta2 - 1.0 ) * ( zeta2 - 1.0 );
            td2NdXi2( 1, 20 ) = - 2.0 * (   xi2 - 1.0 ) * ( zeta2 - 1.0 );
            td2NdXi2( 2, 20 ) = - 2.0 * (  eta2 - 1.0 ) * (   xi2 - 1.0 );
            td2NdXi2( 3, 20 ) = - 4.0 * a * (   xi2 - 1.0 );
            td2NdXi2( 4, 20 ) = - 4.0 * b * (  eta2 - 1.0 );
            td2NdXi2( 5, 20 ) = - 4.0 * c * ( zeta2 - 1.0 );

            td2NdXi2( 0, 21 ) = zeta * (  eta2 - 1.0 ) * ( zeta - 1.0 );
            td2NdXi2( 1, 21 ) = zeta * (   xi2 - 1.0 ) * ( zeta - 1.0 );
            td2NdXi2( 2, 21 ) = (  eta2 - 1.0 ) * (   xi2 - 1.0 );
            td2NdXi2( 3, 21 ) =  eta * (   xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 );
            td2NdXi2( 4, 21 ) =   xi * (  eta2 - 1.0 ) * ( 2.0 * zeta - 1.0 );
            td2NdXi2( 5, 21 ) = d * ( zeta - 1.0 );

            td2NdXi2( 0, 22 ) = zeta * (  eta2 - 1.0 ) * ( zeta + 1.0 );
            td2NdXi2( 1, 22 ) = zeta * (   xi2 - 1.0 ) * ( zeta + 1.0 );
            td2NdXi2( 2, 22 ) = (  eta2 - 1.0 ) * (   xi2 - 1.0 );
            td2NdXi2( 3, 22 ) =  eta * (   xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 );
            td2NdXi2( 4, 22 ) =   xi * (  eta2 - 1.0 ) * ( 2.0 * zeta + 1.0 );
            td2NdXi2( 5, 22 ) = d * ( zeta + 1.0 );

            td2NdXi2( 0, 23 ) = (  eta2 - 1.0 ) * ( zeta2 - 1.0 );
            td2NdXi2( 1, 23 ) =   xi * ( zeta2 - 1.0 ) * (   xi - 1.0 );
            td2NdXi2( 2, 23 ) =   xi * (  eta2 - 1.0 ) * (   xi - 1.0 );
            td2NdXi2( 3, 23 ) = d * (   xi - 1.0 );
            td2NdXi2( 4, 23 ) = zeta * (  eta2 - 1.0 ) * ( 2.0 *   xi - 1.0 );
            td2NdXi2( 5, 23 ) =  eta * ( 2.0 *   xi - 1.0 ) * ( zeta2 - 1.0 );

            td2NdXi2( 0, 24 ) = (  eta2 - 1.0 ) * ( zeta2 - 1.0 );
            td2NdXi2( 1, 24 ) =   xi * ( zeta2 - 1.0 ) * (   xi + 1.0 );
            td2NdXi2( 2, 24 ) =   xi * (  eta2 - 1.0 ) * (   xi + 1.0 );
            td2NdXi2( 3, 24 ) = d * (   xi + 1.0 );
            td2NdXi2( 4, 24 ) = zeta * (  eta2 - 1.0 ) * ( 2.0 *   xi + 1.0 );
            td2NdXi2( 5, 24 ) =  eta * ( 2.0 *   xi + 1.0 ) * ( zeta2 - 1.0 );

            td2NdXi2( 0, 25 ) =  eta * ( zeta2 - 1.0 ) * (  eta - 1.0 );
            td2NdXi2( 1, 25 ) = (   xi2 - 1.0 ) * ( zeta2 - 1.0 );
            td2NdXi2( 2, 25 ) =  eta * (   xi2 - 1.0 ) * (  eta - 1.0 );
            td2NdXi2( 3, 25 ) = zeta * ( 2.0 *  eta - 1.0 ) * (   xi2 - 1.0 );
            td2NdXi2( 4, 25 ) = d * (  eta - 1.0 );
            td2NdXi2( 5, 25 ) =   xi * ( 2.0 *  eta - 1.0 ) * ( zeta2 - 1.0 );

            td2NdXi2( 0, 26 ) =  eta * ( zeta2 - 1.0 ) * (  eta + 1.0 );
            td2NdXi2( 1, 26 ) = (   xi2 - 1.0 ) * ( zeta2 - 1.0 );
            td2NdXi2( 2, 26 ) =  eta * (   xi2 - 1.0 ) * (  eta + 1.0 );
            td2NdXi2( 3, 26 ) = zeta * ( 2.0 *  eta + 1.0 ) * (   xi2 - 1.0 );
            td2NdXi2( 4, 26 ) = d * (  eta + 1.0 );
            td2NdXi2( 5, 26 ) =   xi * ( 2.0 *  eta + 1.0 ) * ( zeta2 - 1.0 );
            return td2NdXi2;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 27  >::eval_d3NdXi3( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "HEX27 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // often used parameters
            // 1st dimension
            real a0 = 0.5 * ( std::pow(xi,2) - xi );
            real a1 = 1.0 - std::pow(xi,2);
            real a2 = 0.5 * ( std::pow(xi,2) + xi );

            real da0 =   xi - 0.5;
            real da1 = - xi * 2.0;
            real da2 =   xi + 0.5;

            real dda0 =   1.0;
            real dda1 = - 2.0;
            real dda2 =   1.0;

            // 2nd dimension
            real b0 = 0.5 * ( std::pow(xi,2) - eta );
            real b1 = 1.0 - std::pow(eta,2);
            real b2 = 0.5 * ( std::pow(xi,2) + eta );

            real db0 =   eta - 0.5;
            real db1 = - eta * 2.0;
            real db2 =   eta + 0.5;

            real ddb0 =   1.0;
            real ddb1 = - 2.0;
            real ddb2 =   1.0;

            // 3rd dimension
            real c0 = 0.5 * ( std::pow(xi,2) - zeta );
            real c1 = 1.0 - std::pow(zeta,2);
            real c2 = 0.5 * ( std::pow(xi,2) + zeta );

            real dc0 =   zeta - 0.5;
            real dc1 = - zeta * 2.0;
            real dc2 =   zeta + 0.5;

            real ddc0 =   1.0;
            real ddc1 = - 2.0;
            real ddc2 =   1.0;

            // 3rd derivatives are = 0 for all dimensions


            Matrix< DDRMat > td3NdXi3(10,27,0.0);

            // 0th node: (0,0,0)
            td3NdXi3( 3,  0 ) =  dda0*  db0*   c0;
            td3NdXi3( 4,  0 ) =  dda0*   b0*  dc0;
            td3NdXi3( 5,  0 ) =   da0* ddb0*   c0;
            td3NdXi3( 6,  0 ) =    a0* ddb0*  dc0;
            td3NdXi3( 7,  0 ) =   da0*   b0* ddc0;
            td3NdXi3( 8,  0 ) =    a0*  db0* ddc0;
            td3NdXi3( 9,  0 ) =   da0*  db0*  dc0;

            // 1th node: (2,0,0)
            td3NdXi3( 3,  1 ) =  dda2*  db0*   c0;
            td3NdXi3( 4,  1 ) =  dda2*   b0*  dc0;
            td3NdXi3( 5,  1 ) =   da2* ddb0*   c0;
            td3NdXi3( 6,  1 ) =    a2* ddb0*  dc0;
            td3NdXi3( 7,  1 ) =   da2*   b0* ddc0;
            td3NdXi3( 8,  1 ) =    a2*  db0* ddc0;
            td3NdXi3( 9,  1 ) =   da2*  db0*  dc0;

            // 2th node: (2,2,0)
            td3NdXi3( 3,  2 ) =  dda2*  db2*   c0;
            td3NdXi3( 4,  2 ) =  dda2*   b2*  dc0;
            td3NdXi3( 5,  2 ) =   da2* ddb2*   c0;
            td3NdXi3( 6,  2 ) =    a2* ddb2*  dc0;
            td3NdXi3( 7,  2 ) =   da2*   b2* ddc0;
            td3NdXi3( 8,  2 ) =    a2*  db2* ddc0;
            td3NdXi3( 9,  2 ) =   da2*  db2*  dc0;

            // 3th node: (0,2,0)
            td3NdXi3( 3,  3 ) =  dda0*  db2*   c0;
            td3NdXi3( 4,  3 ) =  dda0*   b2*  dc0;
            td3NdXi3( 5,  3 ) =   da0* ddb2*   c0;
            td3NdXi3( 6,  3 ) =    a0* ddb2*  dc0;
            td3NdXi3( 7,  3 ) =   da0*   b2* ddc0;
            td3NdXi3( 8,  3 ) =    a0*  db2* ddc0;
            td3NdXi3( 9,  3 ) =   da0*  db2*  dc0;

            // 4th node: (0,0,2)
            td3NdXi3( 3,  4 ) =  dda0*  db0*   c2;
            td3NdXi3( 4,  4 ) =  dda0*   b0*  dc2;
            td3NdXi3( 5,  4 ) =   da0* ddb0*   c2;
            td3NdXi3( 6,  4 ) =    a0* ddb0*  dc2;
            td3NdXi3( 7,  4 ) =   da0*   b0* ddc2;
            td3NdXi3( 8,  4 ) =    a0*  db0* ddc2;
            td3NdXi3( 9,  4 ) =   da0*  db0*  dc2;

            // 5th node: (2,0,2)
            td3NdXi3( 3,  5 ) =  dda2*  db0*   c2;
            td3NdXi3( 4,  5 ) =  dda2*   b0*  dc2;
            td3NdXi3( 5,  5 ) =   da2* ddb0*   c2;
            td3NdXi3( 6,  5 ) =    a2* ddb0*  dc2;
            td3NdXi3( 7,  5 ) =   da2*   b0* ddc2;
            td3NdXi3( 8,  5 ) =    a2*  db0* ddc2;
            td3NdXi3( 9,  5 ) =   da2*  db0*  dc2;

            // 6th node: (2,2,2)
            td3NdXi3( 3,  6 ) =  dda2*  db2*   c2;
            td3NdXi3( 4,  6 ) =  dda2*   b2*  dc2;
            td3NdXi3( 5,  6 ) =   da2* ddb2*   c2;
            td3NdXi3( 6,  6 ) =    a2* ddb2*  dc2;
            td3NdXi3( 7,  6 ) =   da2*   b2* ddc2;
            td3NdXi3( 8,  6 ) =    a2*  db2* ddc2;
            td3NdXi3( 9,  6 ) =   da2*  db2*  dc2;

            // 7th node: (0,2,2)
            td3NdXi3( 3,  7 ) =  dda0*  db2*   c2;
            td3NdXi3( 4,  7 ) =  dda0*   b2*  dc2;
            td3NdXi3( 5,  7 ) =   da0* ddb2*   c2;
            td3NdXi3( 6,  7 ) =    a0* ddb2*  dc2;
            td3NdXi3( 7,  7 ) =   da0*   b2* ddc2;
            td3NdXi3( 8,  7 ) =    a0*  db2* ddc2;
            td3NdXi3( 9,  7 ) =   da0*  db2*  dc2;

            // 8th node: (1,0,0)
            td3NdXi3( 3,  8 ) =  dda1*  db0*   c0;
            td3NdXi3( 4,  8 ) =  dda1*   b0*  dc0;
            td3NdXi3( 5,  8 ) =   da1* ddb0*   c0;
            td3NdXi3( 6,  8 ) =    a1* ddb0*  dc0;
            td3NdXi3( 7,  8 ) =   da1*   b0* ddc0;
            td3NdXi3( 8,  8 ) =    a1*  db0* ddc0;
            td3NdXi3( 9,  8 ) =   da1*  db0*  dc0;

            // 9th node: (2,1,0)
            td3NdXi3( 3,  9 ) =  dda2*  db1*   c0;
            td3NdXi3( 4,  9 ) =  dda2*   b1*  dc0;
            td3NdXi3( 5,  9 ) =   da2* ddb1*   c0;
            td3NdXi3( 6,  9 ) =    a2* ddb1*  dc0;
            td3NdXi3( 7,  9 ) =   da2*   b1* ddc0;
            td3NdXi3( 8,  9 ) =    a2*  db1* ddc0;
            td3NdXi3( 9,  9 ) =   da2*  db1*  dc0;

            // 10th node: (1,2,0)
            td3NdXi3( 3, 10 ) =  dda1*  db2*   c0;
            td3NdXi3( 4, 10 ) =  dda1*   b2*  dc0;
            td3NdXi3( 5, 10 ) =   da1* ddb2*   c0;
            td3NdXi3( 6, 10 ) =    a1* ddb2*  dc0;
            td3NdXi3( 7, 10 ) =   da1*   b2* ddc0;
            td3NdXi3( 8, 10 ) =    a1*  db2* ddc0;
            td3NdXi3( 9, 10 ) =   da1*  db2*  dc0;

            // 11th node: (0,1,0)
            td3NdXi3( 3, 11 ) =  dda0*  db1*   c0;
            td3NdXi3( 4, 11 ) =  dda0*   b1*  dc0;
            td3NdXi3( 5, 11 ) =   da0* ddb1*   c0;
            td3NdXi3( 6, 11 ) =    a0* ddb1*  dc0;
            td3NdXi3( 7, 11 ) =   da0*   b1* ddc0;
            td3NdXi3( 8, 11 ) =    a0*  db1* ddc0;
            td3NdXi3( 9, 11 ) =   da0*  db1*  dc0;

            // 12th node: (0,0,1)
            td3NdXi3( 3, 12 ) =  dda0*  db0*   c1;
            td3NdXi3( 4, 12 ) =  dda0*   b0*  dc1;
            td3NdXi3( 5, 12 ) =   da0* ddb0*   c1;
            td3NdXi3( 6, 12 ) =    a0* ddb0*  dc1;
            td3NdXi3( 7, 12 ) =   da0*   b0* ddc1;
            td3NdXi3( 8, 12 ) =    a0*  db0* ddc1;
            td3NdXi3( 9, 12 ) =   da0*  db0*  dc1;

            // 13th node: (2,0,1)
            td3NdXi3( 3, 13 ) =  dda2*  db0*   c1;
            td3NdXi3( 4, 13 ) =  dda2*   b0*  dc1;
            td3NdXi3( 5, 13 ) =   da2* ddb0*   c1;
            td3NdXi3( 6, 13 ) =    a2* ddb0*  dc1;
            td3NdXi3( 7, 13 ) =   da2*   b0* ddc1;
            td3NdXi3( 8, 13 ) =    a2*  db0* ddc1;
            td3NdXi3( 9, 13 ) =   da2*  db0*  dc1;

            // 14th node: (2,2,1)
            td3NdXi3( 3, 14 ) =  dda2*  db2*   c1;
            td3NdXi3( 4, 14 ) =  dda2*   b2*  dc1;
            td3NdXi3( 5, 14 ) =   da2* ddb2*   c1;
            td3NdXi3( 6, 14 ) =    a2* ddb2*  dc1;
            td3NdXi3( 7, 14 ) =   da2*   b2* ddc1;
            td3NdXi3( 8, 14 ) =    a2*  db2* ddc1;
            td3NdXi3( 9, 14 ) =   da2*  db2*  dc1;

            // 15th node: (0,2,1)
            td3NdXi3( 3, 15 ) =  dda0*  db2*   c1;
            td3NdXi3( 4, 15 ) =  dda0*   b2*  dc1;
            td3NdXi3( 5, 15 ) =   da0* ddb2*   c1;
            td3NdXi3( 6, 15 ) =    a0* ddb2*  dc1;
            td3NdXi3( 7, 15 ) =   da0*   b2* ddc1;
            td3NdXi3( 8, 15 ) =    a0*  db2* ddc1;
            td3NdXi3( 9, 15 ) =   da0*  db2*  dc1;

            // 16th node: (1,0,2)
            td3NdXi3( 3, 16 ) =  dda1*  db0*   c2;
            td3NdXi3( 4, 16 ) =  dda1*   b0*  dc2;
            td3NdXi3( 5, 16 ) =   da1* ddb0*   c2;
            td3NdXi3( 6, 16 ) =    a1* ddb0*  dc2;
            td3NdXi3( 7, 16 ) =   da1*   b0* ddc2;
            td3NdXi3( 8, 16 ) =    a1*  db0* ddc2;
            td3NdXi3( 9, 16 ) =   da1*  db0*  dc2;

            // 17th node: (2,1,2)
            td3NdXi3( 3, 17 ) =  dda2*  db1*   c2;
            td3NdXi3( 4, 17 ) =  dda2*   b1*  dc2;
            td3NdXi3( 5, 17 ) =   da2* ddb1*   c2;
            td3NdXi3( 6, 17 ) =    a2* ddb1*  dc2;
            td3NdXi3( 7, 17 ) =   da2*   b1* ddc2;
            td3NdXi3( 8, 17 ) =    a2*  db1* ddc2;
            td3NdXi3( 9, 17 ) =   da2*  db1*  dc2;

            // 18th node: (1,2,2)
            td3NdXi3( 3, 18 ) =  dda1*  db2*   c2;
            td3NdXi3( 4, 18 ) =  dda1*   b2*  dc2;
            td3NdXi3( 5, 18 ) =   da1* ddb2*   c2;
            td3NdXi3( 6, 18 ) =    a1* ddb2*  dc2;
            td3NdXi3( 7, 18 ) =   da1*   b2* ddc2;
            td3NdXi3( 8, 18 ) =    a1*  db2* ddc2;
            td3NdXi3( 9, 18 ) =   da1*  db2*  dc2;

            // 19th node: (0,1,2)
            td3NdXi3( 3, 19 ) =  dda0*  db1*   c2;
            td3NdXi3( 4, 19 ) =  dda0*   b1*  dc2;
            td3NdXi3( 5, 19 ) =   da0* ddb1*   c2;
            td3NdXi3( 6, 19 ) =    a0* ddb1*  dc2;
            td3NdXi3( 7, 19 ) =   da0*   b1* ddc2;
            td3NdXi3( 8, 19 ) =    a0*  db1* ddc2;
            td3NdXi3( 9, 19 ) =   da0*  db1*  dc2;

            // 20th node: (1,1,1)
            td3NdXi3( 3, 20 ) =  dda1*  db1*   c1;
            td3NdXi3( 4, 20 ) =  dda1*   b1*  dc1;
            td3NdXi3( 5, 20 ) =   da1* ddb1*   c1;
            td3NdXi3( 6, 20 ) =    a1* ddb1*  dc1;
            td3NdXi3( 7, 20 ) =   da1*   b1* ddc1;
            td3NdXi3( 8, 20 ) =    a1*  db1* ddc1;
            td3NdXi3( 9, 20 ) =   da1*  db1*  dc1;

            // 21th node: (1,1,0)
            td3NdXi3( 3, 21 ) =  dda1*  db1*   c0;
            td3NdXi3( 4, 21 ) =  dda1*   b1*  dc0;
            td3NdXi3( 5, 21 ) =   da1* ddb1*   c0;
            td3NdXi3( 6, 21 ) =    a1* ddb1*  dc0;
            td3NdXi3( 7, 21 ) =   da1*   b1* ddc0;
            td3NdXi3( 8, 21 ) =    a1*  db1* ddc0;
            td3NdXi3( 9, 21 ) =   da1*  db1*  dc0;

            // 22th node: (1,1,2)
            td3NdXi3( 3, 22 ) =  dda1*  db1*   c2;
            td3NdXi3( 4, 22 ) =  dda1*   b1*  dc2;
            td3NdXi3( 5, 22 ) =   da1* ddb1*   c2;
            td3NdXi3( 6, 22 ) =    a1* ddb1*  dc2;
            td3NdXi3( 7, 22 ) =   da1*   b1* ddc2;
            td3NdXi3( 8, 22 ) =    a1*  db1* ddc2;
            td3NdXi3( 9, 22 ) =   da1*  db1*  dc2;

            // 23th node: (0,1,1)
            td3NdXi3( 3, 23 ) =  dda0*  db1*   c1;
            td3NdXi3( 4, 23 ) =  dda0*   b1*  dc1;
            td3NdXi3( 5, 23 ) =   da0* ddb1*   c1;
            td3NdXi3( 6, 23 ) =    a0* ddb1*  dc1;
            td3NdXi3( 7, 23 ) =   da0*   b1* ddc1;
            td3NdXi3( 8, 23 ) =    a0*  db1* ddc1;
            td3NdXi3( 9, 23 ) =   da0*  db1*  dc1;

            // 24th node: (2,1,1)
            td3NdXi3( 3, 24 ) =  dda2*  db1*   c1;
            td3NdXi3( 4, 24 ) =  dda2*   b1*  dc1;
            td3NdXi3( 5, 24 ) =   da2* ddb1*   c1;
            td3NdXi3( 6, 24 ) =    a2* ddb1*  dc1;
            td3NdXi3( 7, 24 ) =   da2*   b1* ddc1;
            td3NdXi3( 8, 24 ) =    a2*  db1* ddc1;
            td3NdXi3( 9, 24 ) =   da2*  db1*  dc1;

            // 25th node: (1,0,1)
            td3NdXi3( 3, 25 ) =  dda1*  db0*   c1;
            td3NdXi3( 4, 25 ) =  dda1*   b0*  dc1;
            td3NdXi3( 5, 25 ) =   da1* ddb0*   c1;
            td3NdXi3( 6, 25 ) =    a1* ddb0*  dc1;
            td3NdXi3( 7, 25 ) =   da1*   b0* ddc1;
            td3NdXi3( 8, 25 ) =    a1*  db0* ddc1;
            td3NdXi3( 9, 25 ) =   da1*  db0*  dc1;

            // 26th node: (1,2,1)
            td3NdXi3( 3, 26 ) =  dda1*  db2*   c1;
            td3NdXi3( 4, 26 ) =  dda1*   b2*  dc1;
            td3NdXi3( 5, 26 ) =   da1* ddb2*   c1;
            td3NdXi3( 6, 26 ) =    a1* ddb2*  dc1;
            td3NdXi3( 7, 26 ) =   da1*   b2* ddc1;
            td3NdXi3( 8, 26 ) =    a1*  db2* ddc1;
            td3NdXi3( 9, 26 ) =   da1*  db2*  dc1;

            return td3NdXi3;
        }
//------------------------------------------------------------------------------


    } /* namespace fem */
} /* namespace moris */



#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX27_HPP_ */
