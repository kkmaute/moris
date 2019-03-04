/*
 * cl_FEM_Interpolation_Function_Lagrange_Hex20.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX20_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX20_HPP_

#include "assert.h"

//#include "cl_FEM_Interpolation_Matrix.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src

// checked against femdoc:
// - N
// - dNdxi
// - d2Ndxi2

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 20 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::SERENDIPITY;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 20 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat(3,20);
            tXiHat( 0,  0 ) = -1.000000;
            tXiHat( 1,  0 ) = -1.000000;
            tXiHat( 2,  0 ) = -1.000000;
            tXiHat( 0,  1 ) =  1.000000;
            tXiHat( 1,  1 ) = -1.000000;
            tXiHat( 2,  1 ) = -1.000000;
            tXiHat( 0,  2 ) =  1.000000;
            tXiHat( 1,  2 ) =  1.000000;
            tXiHat( 2,  2 ) = -1.000000;
            tXiHat( 0,  3 ) = -1.000000;
            tXiHat( 1,  3 ) =  1.000000;
            tXiHat( 2,  3 ) = -1.000000;
            tXiHat( 0,  4 ) = -1.000000;
            tXiHat( 1,  4 ) = -1.000000;
            tXiHat( 2,  4 ) =  1.000000;
            tXiHat( 0,  5 ) =  1.000000;
            tXiHat( 1,  5 ) = -1.000000;
            tXiHat( 2,  5 ) =  1.000000;
            tXiHat( 0,  6 ) =  1.000000;
            tXiHat( 1,  6 ) =  1.000000;
            tXiHat( 2,  6 ) =  1.000000;
            tXiHat( 0,  7 ) = -1.000000;
            tXiHat( 1,  7 ) =  1.000000;
            tXiHat( 2,  7 ) =  1.000000;
            tXiHat( 0,  8 ) =  0.000000;
            tXiHat( 1,  8 ) = -1.000000;
            tXiHat( 2,  8 ) = -1.000000;
            tXiHat( 0,  9 ) =  1.000000;
            tXiHat( 1,  9 ) =  0.000000;
            tXiHat( 2,  9 ) = -1.000000;
            tXiHat( 0, 10 ) =  0.000000;
            tXiHat( 1, 10 ) =  1.000000;
            tXiHat( 2, 10 ) = -1.000000;
            tXiHat( 0, 11 ) = -1.000000;
            tXiHat( 1, 11 ) =  0.000000;
            tXiHat( 2, 11 ) = -1.000000;
            tXiHat( 0, 12 ) = -1.000000;
            tXiHat( 1, 12 ) = -1.000000;
            tXiHat( 2, 12 ) =  0.000000;
            tXiHat( 0, 13 ) =  1.000000;
            tXiHat( 1, 13 ) = -1.000000;
            tXiHat( 2, 13 ) =  0.000000;
            tXiHat( 0, 14 ) =  1.000000;
            tXiHat( 1, 14 ) =  1.000000;
            tXiHat( 2, 14 ) =  0.000000;
            tXiHat( 0, 15 ) = -1.000000;
            tXiHat( 1, 15 ) =  1.000000;
            tXiHat( 2, 15 ) =  0.000000;
            tXiHat( 0, 16 ) =  0.000000;
            tXiHat( 1, 16 ) = -1.000000;
            tXiHat( 2, 16 ) =  1.000000;
            tXiHat( 0, 17 ) =  1.000000;
            tXiHat( 1, 17 ) =  0.000000;
            tXiHat( 2, 17 ) =  1.000000;
            tXiHat( 0, 18 ) =  0.000000;
            tXiHat( 1, 18 ) =  1.000000;
            tXiHat( 2, 18 ) =  1.000000;
            tXiHat( 0, 19 ) = -1.000000;
            tXiHat( 1, 19 ) =  0.000000;
            tXiHat( 2, 19 ) =  1.000000;
            tXiHat( 0, 20 ) =  0.000000;
            tXiHat( 1, 20 ) =  0.000000;
            tXiHat( 2, 20 ) =  0.000000;
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 20 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
             // make sure that input is correct
             MORIS_ASSERT( aXi.length() >= 3,
                           "eval_shape: aXi not allocated or hat wrong size." );

             // unpack xi and eta from input vector
             auto xi = aXi( 0 );
             auto eta = aXi( 1 );
             auto  zeta = aXi( 2 );

             // often used constants
             auto xi2 = std::pow(   xi, 2 );
             auto eta2 = std::pow(  eta, 2 );
             auto zeta2 = std::pow( zeta, 2 );

             // populate output matrix
             Matrix< DDRMat > tN(1,20);
             tN(  0 ) =   ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * ( eta + xi + zeta + 2.0 ) * 0.125;
             tN(  1 ) = - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * ( eta - xi + zeta + 2.0 ) * 0.125;
             tN(  2 ) = - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * ( eta + xi - zeta - 2.0 ) * 0.125;
             tN(  3 ) = - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * (  - eta + xi + zeta + 2.0 ) * 0.125;
             tN(  4 ) = - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * ( eta + xi - zeta + 2.0 ) * 0.125;
             tN(  5 ) =   ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * ( eta - xi - zeta + 2.0 ) * 0.125;
             tN(  6 ) =   ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * ( eta + xi + zeta - 2.0 ) * 0.125;
             tN(  7 ) = - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * ( eta - xi + zeta - 2.0 ) * 0.125;
             tN(  8 ) = - ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 ) * 0.25;
             tN(  9 ) =   ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.25;
             tN( 10 ) =   ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 ) * 0.25;
             tN( 11 ) = - ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.25;
             tN( 12 ) = - ( zeta2- 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 ) * 0.25;
             tN( 13 ) =   ( zeta2- 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 ) * 0.25;
             tN( 14 ) = - ( zeta2- 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 ) * 0.25;
             tN( 15 ) =   ( zeta2- 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 ) * 0.25;
             tN( 16 ) =   ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 ) * 0.25;
             tN( 17 ) = - ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.25;
             tN( 18 ) = - ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 ) * 0.25;
             tN( 19 ) =   ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.25;
             return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 20 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                          "eval_shape: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // often used constants
            auto   xi2 = std::pow(   xi, 2 );
            auto  eta2 = std::pow(  eta, 2 );
            auto zeta2 = std::pow( zeta, 2 );

            Matrix< DDRMat > tdNdXi(3,20);
            tdNdXi( 0,  0 ) =   ( (  eta - 1.0 ) * ( zeta - 1.0 ) * (  eta + 2.0 *   xi + zeta + 1.0 ) ) * 0.125;
            tdNdXi( 1,  0 ) =   ( (   xi - 1.0 ) * ( zeta - 1.0 ) * ( 2.0 *  eta +   xi + zeta + 1.0 ) ) * 0.125;
            tdNdXi( 2,  0 ) =   ( (  eta - 1.0 ) * (   xi - 1.0 ) * (  eta +   xi + 2.0 * zeta + 1.0 ) ) * 0.125;

            tdNdXi( 0,  1 ) = - ( (  eta - 1.0 ) * ( zeta - 1.0 ) * (  eta - 2.0 *   xi + zeta + 1.0 ) ) * 0.125;
            tdNdXi( 1,  1 ) = - ( (   xi + 1.0 ) * ( zeta - 1.0 ) * ( 2.0 *  eta -   xi + zeta + 1.0 ) ) * 0.125;
            tdNdXi( 2,  1 ) = - ( (  eta - 1.0 ) * (   xi + 1.0 ) * (  eta -   xi + 2.0 * zeta + 1.0 ) ) * 0.125;

            tdNdXi( 0,  2 ) = - ( (  eta + 1.0 ) * ( zeta - 1.0 ) * (  eta + 2.0 *   xi - zeta - 1.0 ) ) * 0.125;
            tdNdXi( 1,  2 ) = - ( (   xi + 1.0 ) * ( zeta - 1.0 ) * ( 2.0 *  eta +   xi - zeta - 1.0 ) ) * 0.125;
            tdNdXi( 2,  2 ) = - ( (  eta + 1.0 ) * (   xi + 1.0 ) * (  eta +   xi - 2.0 * zeta - 1.0 ) ) * 0.125;

            tdNdXi( 0,  3 ) = - ( (  eta + 1.0 ) * ( zeta - 1.0 ) * ( 2.0 *   xi -  eta + zeta + 1.0 ) ) * 0.125;
            tdNdXi( 1,  3 ) = - ( (   xi - 1.0 ) * ( zeta - 1.0 ) * (   xi - 2.0 *  eta + zeta + 1.0 ) ) * 0.125;
            tdNdXi( 2,  3 ) = - ( (  eta + 1.0 ) * (   xi - 1.0 ) * (   xi -  eta + 2.0 * zeta + 1.0 ) ) * 0.125;

            tdNdXi( 0,  4 ) = - ( (  eta - 1.0 ) * ( zeta + 1.0 ) * (  eta + 2.0 *   xi - zeta + 1.0 ) ) * 0.125;
            tdNdXi( 1,  4 ) = - ( (   xi - 1.0 ) * ( zeta + 1.0 ) * ( 2.0 *  eta +   xi - zeta + 1.0 ) ) * 0.125;
            tdNdXi( 2,  4 ) = - ( (  eta - 1.0 ) * (   xi - 1.0 ) * (  eta +   xi - 2.0 * zeta + 1.0 ) ) * 0.125;

            tdNdXi( 0,  5 )=   ( (  eta - 1.0 ) * ( zeta + 1.0 ) * (  eta - 2.0 *   xi - zeta + 1.0 ) ) * 0.125;
            tdNdXi( 1,  5 )=   ( (   xi + 1.0 ) * ( zeta + 1.0 ) * ( 2.0 *  eta -   xi - zeta + 1.0 ) ) * 0.125;
            tdNdXi( 2,  5 )=   ( (  eta - 1.0 ) * (   xi + 1.0 ) * (  eta -   xi - 2.0 * zeta + 1.0 ) ) * 0.125;

            tdNdXi( 0,  6 )=   ( (  eta + 1.0 ) * ( zeta + 1.0 ) * (  eta + 2.0 *   xi + zeta - 1.0 ) ) * 0.125;
            tdNdXi( 1,  6 )=   ( (   xi + 1.0 ) * ( zeta + 1.0 ) * ( 2.0 *  eta +   xi + zeta - 1.0 ) ) * 0.125;
            tdNdXi( 2,  6 )=   ( (  eta + 1.0 ) * (   xi + 1.0 ) * (  eta +   xi + 2.0 * zeta - 1.0 ) ) * 0.125;

            tdNdXi( 0,  7 )= - ( (  eta + 1.0 ) * ( zeta + 1.0 ) * (  eta - 2.0 *   xi + zeta - 1.0 ) ) * 0.125;
            tdNdXi( 1,  7 )= - ( (   xi - 1.0 ) * ( zeta + 1.0 ) * ( 2.0 *  eta -   xi + zeta - 1.0 ) ) * 0.125;
            tdNdXi( 2,  7 )= - ( (  eta + 1.0 ) * (   xi - 1.0 ) * (  eta -   xi + 2.0 * zeta - 1.0 ) ) * 0.125;

            tdNdXi( 0,  8 )= - (   xi * (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            tdNdXi( 1,  8 )= - ( (   xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            tdNdXi( 2,  8 )= - ( (   xi2 - 1.0 ) * (  eta - 1.0 ) ) * 0.25;

            tdNdXi( 0,  9 )=   ( (  eta2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            tdNdXi( 1,  9 )=   (  eta * (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            tdNdXi( 2,  9 )=   ( (  eta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.25;

            tdNdXi( 0, 10 )=   (   xi * (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            tdNdXi( 1, 10 )=   ( (   xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            tdNdXi( 2, 10 )=   ( (   xi2 - 1.0 ) * (  eta + 1.0 ) ) * 0.25;

            tdNdXi( 0, 11 )= - ( (  eta2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            tdNdXi( 1, 11 )= - (  eta * (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            tdNdXi( 2, 11 )= - ( (  eta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.25;

            tdNdXi( 0, 12 )= - ( ( zeta2 - 1.0 ) * (  eta - 1.0 ) ) * 0.25;
            tdNdXi( 1, 12 )= - ( ( zeta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            tdNdXi( 2, 12 )= - ( zeta * (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.5;

            tdNdXi( 0, 13 )=   ( ( zeta2 - 1.0 ) * (  eta - 1.0 ) ) * 0.25;
            tdNdXi( 1, 13 )=   ( ( zeta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            tdNdXi( 2, 13 )=   ( zeta * (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.5;

            tdNdXi( 0, 14 )= - ( ( zeta2 - 1.0 ) * (  eta + 1.0 ) ) * 0.25;
            tdNdXi( 1, 14 )= - ( ( zeta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            tdNdXi( 2, 14 )= - ( zeta * (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.5;

            tdNdXi( 0, 15 )=   ( ( zeta2 - 1.0 ) * (  eta + 1.0 ) ) * 0.25;
            tdNdXi( 1, 15 )=   ( ( zeta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            tdNdXi( 2, 15 )=   ( zeta * (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.5;

            tdNdXi( 0, 16 )=   (   xi * (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            tdNdXi( 1, 16 )=   ( (   xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            tdNdXi( 2, 16 )=   ( (   xi2 - 1.0 ) * (  eta - 1.0 ) ) * 0.25;

            tdNdXi( 0, 17 )= - ( (  eta2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            tdNdXi( 1, 17 )= - (  eta * (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            tdNdXi( 2, 17 )= - ( (  eta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.25;

            tdNdXi( 0, 18 )= - (   xi * (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            tdNdXi( 1, 18 )= - ( (   xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            tdNdXi( 2, 18 )= - ( (   xi2 - 1.0 ) * (  eta + 1.0 ) ) * 0.25;

            tdNdXi( 0, 19 )=   ( (  eta2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            tdNdXi( 1, 19 )=   (  eta * (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            tdNdXi( 2, 19 )=   ( (  eta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            return tdNdXi;
        }
//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 20 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                          "ad2NdXi2: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // often used constants
            auto   xi2 = std::pow(   xi, 2 );
            auto  eta2 = std::pow(  eta, 2 );
            auto zeta2 = std::pow( zeta, 2 );

            Matrix< DDRMat > td2NdXi2(6,20);
            td2NdXi2( 0,  0 ) =   ( (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 1,  0 ) =   ( (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 2,  0 ) =   ( (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            td2NdXi2( 3,  0 ) =   ( (   xi - 1.0 ) * ( 2.0 * ( eta + zeta ) + xi ) ) * 0.125;
            td2NdXi2( 4,  0 ) =   ( (  eta - 1.0 ) * ( 2.0 * ( xi + zeta ) + eta ) ) * 0.125;
            td2NdXi2( 5,  0 ) =   ( ( zeta - 1.0 ) * ( 2.0 * ( xi + eta ) + zeta ) ) * 0.125;

            td2NdXi2( 0,  1 ) =   ( (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 1,  1 ) = - ( (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 2,  1 ) = - ( (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            td2NdXi2( 3,  1 ) = - ( (   xi + 1.0 ) * ( 2.0 * ( eta + zeta ) - xi ) ) * 0.125;
            td2NdXi2( 4,  1 ) = - ( (  eta - 1.0 ) * (  eta - 2.0 *  ( xi - zeta ) ) ) * 0.125;
            td2NdXi2( 5,  1 ) = - ( ( zeta - 1.0 ) * ( 2.0 *  ( eta - xi ) + zeta ) ) * 0.125;

            td2NdXi2( 0,  2 ) = - ( (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 1,  2 ) = - ( (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 2,  2 ) =   ( (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            td2NdXi2( 3,  2 ) = - ( (   xi + 1.0 ) * ( 2.0 *  eta +   xi - 2.0 * zeta ) ) * 0.125;
            td2NdXi2( 4,  2 ) = - ( (  eta + 1.0 ) * (  eta + 2.0 *   xi - 2.0 * zeta ) ) * 0.125;
            td2NdXi2( 5,  2 ) = - ( ( zeta - 1.0 ) * ( 2.0 *  eta + 2.0 *   xi - zeta ) ) * 0.125;

            td2NdXi2( 0,  3 ) = - ( (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 1,  3 ) =   ( (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            td2NdXi2( 2,  3 ) = - ( (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            td2NdXi2( 3,  3 ) = - ( (   xi - 1.0 ) * ( xi - 2.0 * ( eta - zeta ) ) ) * 0.125;
            td2NdXi2( 4,  3 ) = - ( (  eta + 1.0 ) * ( 2.0 * ( xi + zeta ) - eta ) ) * 0.125;
            td2NdXi2( 5,  3 ) = - ( ( zeta - 1.0 ) * ( zeta - 2.0 * ( eta - xi ) ) ) * 0.125;

            td2NdXi2( 0,  4 ) = - ( (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 1,  4 ) = - ( (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 2,  4 ) =   ( (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            td2NdXi2( 3,  4 ) = - ( (   xi - 1.0 ) * ( 2.0 *  eta +   xi - 2.0 * zeta ) ) * 0.125;
            td2NdXi2( 4,  4 ) = - ( (  eta - 1.0 ) * (  eta + 2.0 *   xi - 2.0 * zeta ) ) * 0.125;
            td2NdXi2( 5,  4 ) = - ( ( zeta + 1.0 ) * ( 2.0 *  eta + 2.0 *   xi - zeta ) ) * 0.125;

            td2NdXi2( 0,  5 ) = - ( (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 1,  5 ) =   ( (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 2,  5 ) = - ( (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            td2NdXi2( 3,  5 ) = - ( (   xi + 1.0 ) * ( xi - 2.0 * ( eta - zeta ) ) ) * 0.125;
            td2NdXi2( 4,  5 ) = - ( (  eta - 1.0 ) * ( 2.0 * ( xi + zeta ) - eta ) ) * 0.125;
            td2NdXi2( 5,  5 ) = - ( ( zeta + 1.0 ) * ( zeta - 2.0 * ( eta - xi ) ) ) * 0.125;

            td2NdXi2( 0,  6 ) =   ( (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 1,  6 ) =   ( (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 2,  6 ) =   ( (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            td2NdXi2( 3,  6 ) =   ( (   xi + 1.0 ) * ( 2.0 * ( eta + zeta ) + xi ) ) * 0.125;
            td2NdXi2( 4,  6 ) =   ( (  eta + 1.0 ) * ( 2.0 * ( xi + zeta ) + eta ) ) * 0.125;
            td2NdXi2( 5,  6 ) =   ( ( zeta + 1.0 ) * ( 2.0 * ( xi + eta ) + zeta ) ) * 0.125;

            td2NdXi2( 0,  7 ) =   ( (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 1,  7 ) = - ( (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            td2NdXi2( 2,  7 ) = - ( (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            td2NdXi2( 3,  7 ) = - ( (   xi - 1.0 ) * ( 2.0 * ( eta + zeta ) - xi ) ) * 0.125;
            td2NdXi2( 4,  7 ) = - ( (  eta + 1.0 ) * (  eta - 2.0 *  ( xi - zeta ) ) ) * 0.125;
            td2NdXi2( 5,  7 ) = - ( ( zeta + 1.0 ) * ( 2.0 *  ( eta - xi ) + zeta ) ) * 0.125;

            td2NdXi2( 0,  8 ) = - ( (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 1,  8 ) = 0.0;
            td2NdXi2( 2,  8 ) = 0.0;
            td2NdXi2( 3,  8 ) = 0.25 -   xi2 * 0.25;
            td2NdXi2( 4,  8 ) = - (   xi * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 5,  8 ) = - (   xi * ( zeta - 1.0 ) ) * 0.5;

            td2NdXi2( 0,  9 ) = 0.0;
            td2NdXi2( 1,  9 ) =   ( (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 2,  9 ) = 0.0;
            td2NdXi2( 3,  9 ) =   (  eta * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 4,  9 ) =  eta2 * 0.25 - 0.25;
            td2NdXi2( 5,  9 ) =   (  eta * ( zeta - 1.0 ) ) * 0.5;

            td2NdXi2( 0, 10 ) =   ( (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 1, 10 ) = 0.0;
            td2NdXi2( 2, 10 ) = 0.0;
            td2NdXi2( 3, 10 ) =   xi2 * 0.25 - 0.25;
            td2NdXi2( 4, 10 ) =   (   xi * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 5, 10 ) =   (   xi * ( zeta - 1.0 ) ) * 0.5;

            td2NdXi2( 0, 11 ) = 0.0;
            td2NdXi2( 1, 11 ) = - ( (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            td2NdXi2( 2, 11 ) = 0.0;
            td2NdXi2( 3, 11 ) = - (  eta * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 4, 11 ) = 0.25 -  eta2 * 0.25;
            td2NdXi2( 5, 11 ) = - (  eta * ( zeta - 1.0 ) ) * 0.5;

            td2NdXi2( 0, 12 ) = 0.0;
            td2NdXi2( 1, 12 ) = 0.0;
            td2NdXi2( 2, 12 ) = - ( (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 3, 12 ) = - ( zeta * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 4, 12 ) = - ( zeta * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 5, 12 ) = 0.25 - zeta2 * 0.25;

            td2NdXi2( 0, 13 ) = 0.0;
            td2NdXi2( 1, 13 ) = 0.0;
            td2NdXi2( 2, 13 ) =   ( (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 3, 13 ) =   ( zeta * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 4, 13 ) =   ( zeta * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 5, 13 ) = zeta2 * 0.25 - 0.25;

            td2NdXi2( 0, 14 ) = 0.0;
            td2NdXi2( 1, 14 ) = 0.0;
            td2NdXi2( 2, 14 ) = - ( (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 3, 14 ) = - ( zeta * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 4, 14 ) = - ( zeta * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 5, 14 ) = 0.25 - zeta2 * 0.25;

            td2NdXi2( 0, 15 ) = 0.0;
            td2NdXi2( 1, 15 ) = 0.0;
            td2NdXi2( 2, 15 ) =   ( (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 3, 15 ) =   ( zeta * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 4, 15 ) =   ( zeta * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 5, 15 ) = zeta2 * 0.25 - 0.25;

            td2NdXi2( 0, 16 ) =   ( (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 16 ) = 0.0;
            td2NdXi2( 2, 16 ) = 0.0;
            td2NdXi2( 3, 16 ) =   xi2 * 0.25 - 0.25;
            td2NdXi2( 4, 16 ) =   (   xi * (  eta - 1.0 ) ) * 0.5;
            td2NdXi2( 5, 16 ) =   (   xi * ( zeta + 1.0 ) ) * 0.5;

            td2NdXi2( 0, 17 ) = 0.0;
            td2NdXi2( 1, 17 ) = - ( (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 17 ) = 0.0;
            td2NdXi2( 3, 17 ) = - (  eta * (   xi + 1.0 ) ) * 0.5;
            td2NdXi2( 4, 17 ) = 0.25 -  eta2 * 0.25;
            td2NdXi2( 5, 17 ) = - (  eta * ( zeta + 1.0 ) ) * 0.5;

            td2NdXi2( 0, 18 ) = - ( (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 18 ) = 0.0;
            td2NdXi2( 2, 18 ) = 0.0;
            td2NdXi2( 3, 18 ) = 0.25 -   xi2 * 0.25;
            td2NdXi2( 4, 18 ) = - (   xi * (  eta + 1.0 ) ) * 0.5;
            td2NdXi2( 5, 18 ) = - (   xi * ( zeta + 1.0 ) ) * 0.5;

            td2NdXi2( 0, 19 ) = 0.0;
            td2NdXi2( 1, 19 ) =   ( (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 19 ) = 0.0;
            td2NdXi2( 3, 19 ) =   (  eta * (   xi - 1.0 ) ) * 0.5;
            td2NdXi2( 4, 19 ) =  eta2 * 0.25 - 0.25;
            td2NdXi2( 5, 19 ) =   (  eta * ( zeta + 1.0 ) ) * 0.5;
            return td2NdXi2;
        }
//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */




#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX20_HPP_ */
