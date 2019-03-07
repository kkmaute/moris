/*
 * cl_FEM_Interpolation_Function_Hex64.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX64_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX64_HPP_

#include "assert.h"

//#include "cl_FEM_Interpolation_Matrix.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 64 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 64 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat( 3, 64 );
            real c = 1.0/3.0;
            tXiHat( 0,  0 ) = -1.0;
            tXiHat( 1,  0 ) = -1.0;
            tXiHat( 2,  0 ) = -1.0;
            tXiHat( 0,  1 ) =  1.0;
            tXiHat( 1,  1 ) = -1.0;
            tXiHat( 2,  1 ) = -1.0;
            tXiHat( 0,  2 ) =  1.0;
            tXiHat( 1,  2 ) =  1.0;
            tXiHat( 2,  2 ) = -1.0;
            tXiHat( 0,  3 ) = -1.0;
            tXiHat( 1,  3 ) =  1.0;
            tXiHat( 2,  3 ) = -1.0;
            tXiHat( 0,  4 ) = -1.0;
            tXiHat( 1,  4 ) = -1.0;
            tXiHat( 2,  4 ) =  1.0;
            tXiHat( 0,  5 ) =  1.0;
            tXiHat( 1,  5 ) = -1.0;
            tXiHat( 2,  5 ) =  1.0;
            tXiHat( 0,  6 ) =  1.0;
            tXiHat( 1,  6 ) =  1.0;
            tXiHat( 2,  6 ) =  1.0;
            tXiHat( 0,  7 ) = -1.0;
            tXiHat( 1,  7 ) =  1.0;
            tXiHat( 2,  7 ) =  1.0;
            tXiHat( 0,  8 ) = -c;
            tXiHat( 1,  8 ) = -1.0;
            tXiHat( 2,  8 ) = -1.0;
            tXiHat( 0,  9 ) =  c;
            tXiHat( 1,  9 ) = -1.0;
            tXiHat( 2,  9 ) = -1.0;
            tXiHat( 0, 10 ) = -1.0;
            tXiHat( 1, 10 ) = -c;
            tXiHat( 2, 10 ) = -1.0;
            tXiHat( 0, 11 ) = -1.0;
            tXiHat( 1, 11 ) =  c;
            tXiHat( 2, 11 ) = -1.0;
            tXiHat( 0, 12 ) = -1.0;
            tXiHat( 1, 12 ) = -1.0;
            tXiHat( 2, 12 ) = -c;
            tXiHat( 0, 13 ) = -1.0;
            tXiHat( 1, 13 ) = -1.0;
            tXiHat( 2, 13 ) =  c;
            tXiHat( 0, 14 ) =  1.0;
            tXiHat( 1, 14 ) = -c;
            tXiHat( 2, 14 ) = -1.0;
            tXiHat( 0, 15 ) =  1.0;
            tXiHat( 1, 15 ) =  c;
            tXiHat( 2, 15 ) = -1.0;
            tXiHat( 0, 16 ) =  1.0;
            tXiHat( 1, 16 ) = -1.0;
            tXiHat( 2, 16 ) = -c;
            tXiHat( 0, 17 ) =  1.0;
            tXiHat( 1, 17 ) = -1.0;
            tXiHat( 2, 17 ) =  c;
            tXiHat( 0, 18 ) =  c;
            tXiHat( 1, 18 ) =  1.0;
            tXiHat( 2, 18 ) = -1.0;
            tXiHat( 0, 19 ) = -c;
            tXiHat( 1, 19 ) =  1.0;
            tXiHat( 2, 19 ) = -1.0;
            tXiHat( 0, 20 ) =  1.0;
            tXiHat( 1, 20 ) =  1.0;
            tXiHat( 2, 20 ) = -c;
            tXiHat( 0, 21 ) =  1.0;
            tXiHat( 1, 21 ) =  1.0;
            tXiHat( 2, 21 ) =  c;
            tXiHat( 0, 22 ) = -1.0;
            tXiHat( 1, 22 ) =  1.0;
            tXiHat( 2, 22 ) = -c;
            tXiHat( 0, 23 ) = -1.0;
            tXiHat( 1, 23 ) =  1.0;
            tXiHat( 2, 23 ) =  c;
            tXiHat( 0, 24 ) = -c;
            tXiHat( 1, 24 ) = -1.0;
            tXiHat( 2, 24 ) =  1.0;
            tXiHat( 0, 25 ) =  c;
            tXiHat( 1, 25 ) = -1.0;
            tXiHat( 2, 25 ) =  1.0;
            tXiHat( 0, 26 ) = -1.0;
            tXiHat( 1, 26 ) = -c;
            tXiHat( 2, 26 ) =  1.0;
            tXiHat( 0, 27 ) = -1.0;
            tXiHat( 1, 27 ) =  c;
            tXiHat( 2, 27 ) =  1.0;
            tXiHat( 0, 28 ) =  1.0;
            tXiHat( 1, 28 ) = -c;
            tXiHat( 2, 28 ) =  1.0;
            tXiHat( 0, 29 ) =  1.0;
            tXiHat( 1, 29 ) =  c;
            tXiHat( 2, 29 ) =  1.0;
            tXiHat( 0, 30 ) =  c;
            tXiHat( 1, 30 ) =  1.0;
            tXiHat( 2, 30 ) =  1.0;
            tXiHat( 0, 31 ) = -c;
            tXiHat( 1, 31 ) =  1.0;
            tXiHat( 2, 31 ) =  1.0;
            tXiHat( 0, 32 ) = -c;
            tXiHat( 1, 32 ) = -c;
            tXiHat( 2, 32 ) = -1.0;
            tXiHat( 0, 33 ) = -c;
            tXiHat( 1, 33 ) =  c;
            tXiHat( 2, 33 ) = -1.0;
            tXiHat( 0, 34 ) =  c;
            tXiHat( 1, 34 ) =  c;
            tXiHat( 2, 34 ) = -1.0;
            tXiHat( 0, 35 ) =  c;
            tXiHat( 1, 35 ) = -c;
            tXiHat( 2, 35 ) = -1.0;
            tXiHat( 0, 36 ) = -c;
            tXiHat( 1, 36 ) = -1.0;
            tXiHat( 2, 36 ) = -c;
            tXiHat( 0, 37 ) =  c;
            tXiHat( 1, 37 ) = -1.0;
            tXiHat( 2, 37 ) = -c;
            tXiHat( 0, 38 ) =  c;
            tXiHat( 1, 38 ) = -1.0;
            tXiHat( 2, 38 ) =  c;
            tXiHat( 0, 39 ) = -c;
            tXiHat( 1, 39 ) = -1.0;
            tXiHat( 2, 39 ) =  c;
            tXiHat( 0, 40 ) = -1.0;
            tXiHat( 1, 40 ) = -c;
            tXiHat( 2, 40 ) = -c;
            tXiHat( 0, 41 ) = -1.0;
            tXiHat( 1, 41 ) = -c;
            tXiHat( 2, 41 ) =  c;
            tXiHat( 0, 42 ) = -1.0;
            tXiHat( 1, 42 ) =  c;
            tXiHat( 2, 42 ) =  c;
            tXiHat( 0, 43 ) = -1.0;
            tXiHat( 1, 43 ) =  c;
            tXiHat( 2, 43 ) = -c;
            tXiHat( 0, 44 ) =  1.0;
            tXiHat( 1, 44 ) = -c;
            tXiHat( 2, 44 ) = -c;
            tXiHat( 0, 45 ) =  1.0;
            tXiHat( 1, 45 ) =  c;
            tXiHat( 2, 45 ) = -c;
            tXiHat( 0, 46 ) =  1.0;
            tXiHat( 1, 46 ) =  c;
            tXiHat( 2, 46 ) =  c;
            tXiHat( 0, 47 ) =  1.0;
            tXiHat( 1, 47 ) = -c;
            tXiHat( 2, 47 ) =  c;
            tXiHat( 0, 48 ) =  c;
            tXiHat( 1, 48 ) =  1.0;
            tXiHat( 2, 48 ) = -c;
            tXiHat( 0, 49 ) = -c;
            tXiHat( 1, 49 ) =  1.0;
            tXiHat( 2, 49 ) = -c;
            tXiHat( 0, 50 ) = -c;
            tXiHat( 1, 50 ) =  1.0;
            tXiHat( 2, 50 ) =  c;
            tXiHat( 0, 51 ) =  c;
            tXiHat( 1, 51 ) =  1.0;
            tXiHat( 2, 51 ) =  c;
            tXiHat( 0, 52 ) = -c;
            tXiHat( 1, 52 ) = -c;
            tXiHat( 2, 52 ) =  1.0;
            tXiHat( 0, 53 ) =  c;
            tXiHat( 1, 53 ) = -c;
            tXiHat( 2, 53 ) =  1.0;
            tXiHat( 0, 54 ) =  c;
            tXiHat( 1, 54 ) =  c;
            tXiHat( 2, 54 ) =  1.0;
            tXiHat( 0, 55 ) = -c;
            tXiHat( 1, 55 ) =  c;
            tXiHat( 2, 55 ) =  1.0;
            tXiHat( 0, 56 ) = -c;
            tXiHat( 1, 56 ) = -c;
            tXiHat( 2, 56 ) = -c;
            tXiHat( 0, 57 ) =  c;
            tXiHat( 1, 57 ) = -c;
            tXiHat( 2, 57 ) = -c;
            tXiHat( 0, 58 ) =  c;
            tXiHat( 1, 58 ) =  c;
            tXiHat( 2, 58 ) = -c;
            tXiHat( 0, 59 ) = -c;
            tXiHat( 1, 59 ) =  c;
            tXiHat( 2, 59 ) = -c;
            tXiHat( 0, 60 ) = -c;
            tXiHat( 1, 60 ) = -c;
            tXiHat( 2, 60 ) =  c;
            tXiHat( 0, 61 ) =  c;
            tXiHat( 1, 61 ) = -c;
            tXiHat( 2, 61 ) =  c;
            tXiHat( 0, 62 ) =  c;
            tXiHat( 1, 62 ) =  c;
            tXiHat( 2, 62 ) =  c;
            tXiHat( 0, 63 ) = -c;
            tXiHat( 1, 63 ) =  c;
            tXiHat( 2, 63 ) =  c;
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat>
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 64 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // unpack xi and eta and zeta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

            real c0 =  ( zeta*( 1.0 + 9.0 * zeta * ( 1.0 - zeta ) ) - 1.0 )*0.0625;
            real c1 =  ( 9.0 - zeta * ( 27.0 + zeta*( 9.0 - 27.0*zeta ) ) )*0.0625;
            real c2 =  ( 9.0 + zeta * ( 27.0 - zeta*( 9.0 + 27.0*zeta ) ) )*0.0625;
            real c3 = ( -zeta*( 1.0 - 9.0 * zeta * ( 1.0 + zeta ) ) - 1.0 )*0.0625;

            // populate matrix with values
            Matrix< DDRMat > tN(1,64);
            tN(  0 ) = a0 * b0 * c0;
            tN(  1 ) = a3 * b0 * c0;
            tN(  2 ) = a3 * b3 * c0;
            tN(  3 ) = a0 * b3 * c0;
            tN(  4 ) = a0 * b0 * c3;
            tN(  5 ) = a3 * b0 * c3;
            tN(  6 ) = a3 * b3 * c3;
            tN(  7 ) = a0 * b3 * c3;
            tN(  8 ) = a1 * b0 * c0;
            tN(  9 ) = a2 * b0 * c0;
            tN( 10 ) = a0 * b1 * c0;
            tN( 11 ) = a0 * b2 * c0;
            tN( 12 ) = a0 * b0 * c1;
            tN( 13 ) = a0 * b0 * c2;
            tN( 14 ) = a3 * b1 * c0;
            tN( 15 ) = a3 * b2 * c0;
            tN( 16 ) = a3 * b0 * c1;
            tN( 17 ) = a3 * b0 * c2;
            tN( 18 ) = a2 * b3 * c0;
            tN( 19 ) = a1 * b3 * c0;
            tN( 20 ) = a3 * b3 * c1;
            tN( 21 ) = a3 * b3 * c2;
            tN( 22 ) = a0 * b3 * c1;
            tN( 23 ) = a0 * b3 * c2;
            tN( 24 ) = a1 * b0 * c3;
            tN( 25 ) = a2 * b0 * c3;
            tN( 26 ) = a0 * b1 * c3;
            tN( 27 ) = a0 * b2 * c3;
            tN( 28 ) = a3 * b1 * c3;
            tN( 29 ) = a3 * b2 * c3;
            tN( 30 ) = a2 * b3 * c3;
            tN( 31 ) = a1 * b3 * c3;
            tN( 32 ) = a1 * b1 * c0;
            tN( 33 ) = a1 * b2 * c0;
            tN( 34 ) = a2 * b2 * c0;
            tN( 35 ) = a2 * b1 * c0;
            tN( 36 ) = a1 * b0 * c1;
            tN( 37 ) = a2 * b0 * c1;
            tN( 38 ) = a2 * b0 * c2;
            tN( 39 ) = a1 * b0 * c2;
            tN( 40 ) = a0 * b1 * c1;
            tN( 41 ) = a0 * b1 * c2;
            tN( 42 ) = a0 * b2 * c2;
            tN( 43 ) = a0 * b2 * c1;
            tN( 44 ) = a3 * b1 * c1;
            tN( 45 ) = a3 * b2 * c1;
            tN( 46 ) = a3 * b2 * c2;
            tN( 47 ) = a3 * b1 * c2;
            tN( 48 ) = a2 * b3 * c1;
            tN( 49 ) = a1 * b3 * c1;
            tN( 50 ) = a1 * b3 * c2;
            tN( 51 ) = a2 * b3 * c2;
            tN( 52 ) = a1 * b1 * c3;
            tN( 53 ) = a2 * b1 * c3;
            tN( 54 ) = a2 * b2 * c3;
            tN( 55 ) = a1 * b2 * c3;
            tN( 56 ) = a1 * b1 * c1;
            tN( 57 ) = a2 * b1 * c1;
            tN( 58 ) = a2 * b2 * c1;
            tN( 59 ) = a1 * b2 * c1;
            tN( 60 ) = a1 * b1 * c2;
            tN( 61 ) = a2 * b1 * c2;
            tN( 62 ) = a2 * b2 * c2;
            tN( 63 ) = a1 * b2 * c2;
            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat>
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 64 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                          "eval_shape: aXi not allocated or hat wrong size." );

            // unpack xi and eta and zeta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // often used parameters
            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 ) * 0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) ) * 0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) ) * 0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 ) * 0.0625;

            real c0 =  ( zeta*( 1.0 + 9.0 * zeta * ( 1.0 - zeta ) ) - 1.0 )*0.0625;
            real c1 =  ( 9.0 - zeta * ( 27.0 + zeta*( 9.0 - 27.0*zeta ) ) )*0.0625;
            real c2 =  ( 9.0 + zeta * ( 27.0 - zeta*( 9.0 + 27.0*zeta ) ) )*0.0625;
            real c3 = ( -zeta*( 1.0 - 9.0 * zeta * ( 1.0 + zeta ) ) - 1.0 )*0.0625;

            real da0 = (   1.0 + xi*( 18.0 - 27.0*xi )) * 0.0625;
            real da1 = ( -27.0 - xi*( 18.0 - 81.0*xi )) * 0.0625;
            real da2 = (  27.0 - xi*( 18.0 + 81.0*xi )) * 0.0625;
            real da3 = (  -1.0 + xi*( 18.0 + 27.0*xi )) * 0.0625;

            real db0 = (   1.0 + eta*( 18.0 - 27.0*eta )) * 0.0625;
            real db1 = ( -27.0 - eta*( 18.0 - 81.0*eta )) * 0.0625;
            real db2 = (  27.0 - eta*( 18.0 + 81.0*eta )) * 0.0625;
            real db3 = (  -1.0 + eta*( 18.0 + 27.0*eta )) * 0.0625;

            real dc0 = (   1.0 + zeta*( 18.0 - 27.0*zeta )) * 0.0625;
            real dc1 = ( -27.0 - zeta*( 18.0 - 81.0*zeta )) * 0.0625;
            real dc2 = (  27.0 - zeta*( 18.0 + 81.0*zeta )) * 0.0625;
            real dc3 = (  -1.0 + zeta*( 18.0 + 27.0*zeta )) * 0.0625;

            // populate output matrix
            Matrix< DDRMat > tdNdXi(3,64);
            tdNdXi( 0, 0 ) = b0*c0*da0;
            tdNdXi( 1, 0 ) = a0*c0*db0;
            tdNdXi( 2, 0 ) = a0*b0*dc0;

            tdNdXi( 0, 1 ) = b0*c0*da3;
            tdNdXi( 1, 1 ) = a3*c0*db0;
            tdNdXi( 2, 1 ) = a3*b0*dc0;

            tdNdXi( 0, 2 ) = b3*c0*da3;
            tdNdXi( 1, 2 ) = a3*c0*db3;
            tdNdXi( 2, 2 ) = a3*b3*dc0;

            tdNdXi( 0, 3 ) = b3*c0*da0;
            tdNdXi( 1, 3 ) = a0*c0*db3;
            tdNdXi( 2, 3 ) = a0*b3*dc0;

            tdNdXi( 0, 4 ) = b0*c3*da0;
            tdNdXi( 1, 4 ) = a0*c3*db0;
            tdNdXi( 2, 4 ) = a0*b0*dc3;

            tdNdXi( 0, 5 ) = b0*c3*da3;
            tdNdXi( 1, 5 ) = a3*c3*db0;
            tdNdXi( 2, 5 ) = a3*b0*dc3;

            tdNdXi( 0, 6 ) = b3*c3*da3;
            tdNdXi( 1, 6 ) = a3*c3*db3;
            tdNdXi( 2, 6 ) = a3*b3*dc3;

            tdNdXi( 0, 7 ) = b3*c3*da0;
            tdNdXi( 1, 7 ) = a0*c3*db3;
            tdNdXi( 2, 7 ) = a0*b3*dc3;

            tdNdXi( 0, 8 ) = b0*c0*da1;
            tdNdXi( 1, 8 ) = a1*c0*db0;
            tdNdXi( 2, 8 ) = a1*b0*dc0;

            tdNdXi( 0, 9 ) = b0*c0*da2;
            tdNdXi( 1, 9 ) = a2*c0*db0;
            tdNdXi( 2, 9 ) = a2*b0*dc0;

            tdNdXi( 0, 10 ) = b1*c0*da0;
            tdNdXi( 1, 10 ) = a0*c0*db1;
            tdNdXi( 2, 10 ) = a0*b1*dc0;

            tdNdXi( 0, 11 ) = b2*c0*da0;
            tdNdXi( 1, 11 ) = a0*c0*db2;
            tdNdXi( 2, 11 ) = a0*b2*dc0;

            tdNdXi( 0, 12 ) = b0*c1*da0;
            tdNdXi( 1, 12 ) = a0*c1*db0;
            tdNdXi( 2, 12 ) = a0*b0*dc1;

            tdNdXi( 0, 13 ) = b0*c2*da0;
            tdNdXi( 1, 13 ) = a0*c2*db0;
            tdNdXi( 2, 13 ) = a0*b0*dc2;

            tdNdXi( 0, 14 ) = b1*c0*da3;
            tdNdXi( 1, 14 ) = a3*c0*db1;
            tdNdXi( 2, 14 ) = a3*b1*dc0;

            tdNdXi( 0, 15 ) = b2*c0*da3;
            tdNdXi( 1, 15 ) = a3*c0*db2;
            tdNdXi( 2, 15 ) = a3*b2*dc0;

            tdNdXi( 0, 16 ) = b0*c1*da3;
            tdNdXi( 1, 16 ) = a3*c1*db0;
            tdNdXi( 2, 16 ) = a3*b0*dc1;

            tdNdXi( 0, 17 ) = b0*c2*da3;
            tdNdXi( 1, 17 ) = a3*c2*db0;
            tdNdXi( 2, 17 ) = a3*b0*dc2;

            tdNdXi( 0, 18 ) = b3*c0*da2;
            tdNdXi( 1, 18 ) = a2*c0*db3;
            tdNdXi( 2, 18 ) = a2*b3*dc0;

            tdNdXi( 0, 19 ) = b3*c0*da1;
            tdNdXi( 1, 19 ) = a1*c0*db3;
            tdNdXi( 2, 19 ) = a1*b3*dc0;

            tdNdXi( 0, 20 ) = b3*c1*da3;
            tdNdXi( 1, 20 ) = a3*c1*db3;
            tdNdXi( 2, 20 ) = a3*b3*dc1;

            tdNdXi( 0, 21 ) = b3*c2*da3;
            tdNdXi( 1, 21 ) = a3*c2*db3;
            tdNdXi( 2, 21 ) = a3*b3*dc2;

            tdNdXi( 0, 22 ) = b3*c1*da0;
            tdNdXi( 1, 22 ) = a0*c1*db3;
            tdNdXi( 2, 22 ) = a0*b3*dc1;

            tdNdXi( 0, 23 ) = b3*c2*da0;
            tdNdXi( 1, 23 ) = a0*c2*db3;
            tdNdXi( 2, 23 ) = a0*b3*dc2;

            tdNdXi( 0, 24 ) = b0*c3*da1;
            tdNdXi( 1, 24 ) = a1*c3*db0;
            tdNdXi( 2, 24 ) = a1*b0*dc3;

            tdNdXi( 0, 25 ) = b0*c3*da2;
            tdNdXi( 1, 25 ) = a2*c3*db0;
            tdNdXi( 2, 25 ) = a2*b0*dc3;

            tdNdXi( 0, 26 ) = b1*c3*da0;
            tdNdXi( 1, 26 ) = a0*c3*db1;
            tdNdXi( 2, 26 ) = a0*b1*dc3;

            tdNdXi( 0, 27 ) = b2*c3*da0;
            tdNdXi( 1, 27 ) = a0*c3*db2;
            tdNdXi( 2, 27 ) = a0*b2*dc3;

            tdNdXi( 0, 28 ) = b1*c3*da3;
            tdNdXi( 1, 28 ) = a3*c3*db1;
            tdNdXi( 2, 28 ) = a3*b1*dc3;

            tdNdXi( 0, 29 ) = b2*c3*da3;
            tdNdXi( 1, 29 ) = a3*c3*db2;
            tdNdXi( 2, 29 ) = a3*b2*dc3;

            tdNdXi( 0, 30 ) = b3*c3*da2;
            tdNdXi( 1, 30 ) = a2*c3*db3;
            tdNdXi( 2, 30 ) = a2*b3*dc3;

            tdNdXi( 0, 31 ) = b3*c3*da1;
            tdNdXi( 1, 31 ) = a1*c3*db3;
            tdNdXi( 2, 31 ) = a1*b3*dc3;

            tdNdXi( 0, 32 ) = b1*c0*da1;
            tdNdXi( 1, 32 ) = a1*c0*db1;
            tdNdXi( 2, 32 ) = a1*b1*dc0;

            tdNdXi( 0, 33 ) = b2*c0*da1;
            tdNdXi( 1, 33 ) = a1*c0*db2;
            tdNdXi( 2, 33 ) = a1*b2*dc0;

            tdNdXi( 0, 34 ) = b2*c0*da2;
            tdNdXi( 1, 34 ) = a2*c0*db2;
            tdNdXi( 2, 34 ) = a2*b2*dc0;

            tdNdXi( 0, 35 ) = b1*c0*da2;
            tdNdXi( 1, 35 ) = a2*c0*db1;
            tdNdXi( 2, 35 ) = a2*b1*dc0;

            tdNdXi( 0, 36 ) = b0*c1*da1;
            tdNdXi( 1, 36 ) = a1*c1*db0;
            tdNdXi( 2, 36 ) = a1*b0*dc1;

            tdNdXi( 0, 37 ) = b0*c1*da2;
            tdNdXi( 1, 37 ) = a2*c1*db0;
            tdNdXi( 2, 37 ) = a2*b0*dc1;

            tdNdXi( 0, 38 ) = b0*c2*da2;
            tdNdXi( 1, 38 ) = a2*c2*db0;
            tdNdXi( 2, 38 ) = a2*b0*dc2;

            tdNdXi( 0, 39 ) = b0*c2*da1;
            tdNdXi( 1, 39 ) = a1*c2*db0;
            tdNdXi( 2, 39 ) = a1*b0*dc2;

            tdNdXi( 0, 40 ) = b1*c1*da0;
            tdNdXi( 1, 40 ) = a0*c1*db1;
            tdNdXi( 2, 40 ) = a0*b1*dc1;

            tdNdXi( 0, 41 ) = b1*c2*da0;
            tdNdXi( 1, 41 ) = a0*c2*db1;
            tdNdXi( 2, 41 ) = a0*b1*dc2;

            tdNdXi( 0, 42 ) = b2*c2*da0;
            tdNdXi( 1, 42 ) = a0*c2*db2;
            tdNdXi( 2, 42 ) = a0*b2*dc2;

            tdNdXi( 0, 43 ) = b2*c1*da0;
            tdNdXi( 1, 43 ) = a0*c1*db2;
            tdNdXi( 2, 43 ) = a0*b2*dc1;

            tdNdXi( 0, 44 ) = b1*c1*da3;
            tdNdXi( 1, 44 ) = a3*c1*db1;
            tdNdXi( 2, 44 ) = a3*b1*dc1;

            tdNdXi( 0, 45 ) = b2*c1*da3;
            tdNdXi( 1, 45 ) = a3*c1*db2;
            tdNdXi( 2, 45 ) = a3*b2*dc1;

            tdNdXi( 0, 46 ) = b2*c2*da3;
            tdNdXi( 1, 46 ) = a3*c2*db2;
            tdNdXi( 2, 46 ) = a3*b2*dc2;

            tdNdXi( 0, 47 ) = b1*c2*da3;
            tdNdXi( 1, 47 ) = a3*c2*db1;
            tdNdXi( 2, 47 ) = a3*b1*dc2;

            tdNdXi( 0, 48 ) = b3*c1*da2;
            tdNdXi( 1, 48 ) = a2*c1*db3;
            tdNdXi( 2, 48 ) = a2*b3*dc1;

            tdNdXi( 0, 49 ) = b3*c1*da1;
            tdNdXi( 1, 49 ) = a1*c1*db3;
            tdNdXi( 2, 49 ) = a1*b3*dc1;

            tdNdXi( 0, 50 ) = b3*c2*da1;
            tdNdXi( 1, 50 ) = a1*c2*db3;
            tdNdXi( 2, 50 ) = a1*b3*dc2;

            tdNdXi( 0, 51 ) = b3*c2*da2;
            tdNdXi( 1, 51 ) = a2*c2*db3;
            tdNdXi( 2, 51 ) = a2*b3*dc2;

            tdNdXi( 0, 52 ) = b1*c3*da1;
            tdNdXi( 1, 52 ) = a1*c3*db1;
            tdNdXi( 2, 52 ) = a1*b1*dc3;

            tdNdXi( 0, 53 ) = b1*c3*da2;
            tdNdXi( 1, 53 ) = a2*c3*db1;
            tdNdXi( 2, 53 ) = a2*b1*dc3;

            tdNdXi( 0, 54 ) = b2*c3*da2;
            tdNdXi( 1, 54 ) = a2*c3*db2;
            tdNdXi( 2, 54 ) = a2*b2*dc3;

            tdNdXi( 0, 55 ) = b2*c3*da1;
            tdNdXi( 1, 55 ) = a1*c3*db2;
            tdNdXi( 2, 55 ) = a1*b2*dc3;

            tdNdXi( 0, 56 ) = b1*c1*da1;
            tdNdXi( 1, 56 ) = a1*c1*db1;
            tdNdXi( 2, 56 ) = a1*b1*dc1;

            tdNdXi( 0, 57 ) = b1*c1*da2;
            tdNdXi( 1, 57 ) = a2*c1*db1;
            tdNdXi( 2, 57 ) = a2*b1*dc1;

            tdNdXi( 0, 58 ) = b2*c1*da2;
            tdNdXi( 1, 58 ) = a2*c1*db2;
            tdNdXi( 2, 58 ) = a2*b2*dc1;

            tdNdXi( 0, 59 ) = b2*c1*da1;
            tdNdXi( 1, 59 ) = a1*c1*db2;
            tdNdXi( 2, 59 ) = a1*b2*dc1;

            tdNdXi( 0, 60 ) = b1*c2*da1;
            tdNdXi( 1, 60 ) = a1*c2*db1;
            tdNdXi( 2, 60 ) = a1*b1*dc2;

            tdNdXi( 0, 61 ) = b1*c2*da2;
            tdNdXi( 1, 61 ) = a2*c2*db1;
            tdNdXi( 2, 61 ) = a2*b1*dc2;

            tdNdXi( 0, 62 ) = b2*c2*da2;
            tdNdXi( 1, 62 ) = a2*c2*db2;
            tdNdXi( 2, 62 ) = a2*b2*dc2;

            tdNdXi( 0, 63 ) = b2*c2*da1;
            tdNdXi( 1, 63 ) = a1*c2*db2;
            tdNdXi( 2, 63 ) = a1*b2*dc2;
            return tdNdXi;
        }
//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 64  >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // often used parameters
            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 ) * 0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) ) * 0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) ) * 0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 ) * 0.0625;

            real c0 =  ( zeta*( 1.0 + 9.0 * zeta * ( 1.0 - zeta ) ) - 1.0 )*0.0625;
            real c1 =  ( 9.0 - zeta * ( 27.0 + zeta*( 9.0 - 27.0*zeta ) ) )*0.0625;
            real c2 =  ( 9.0 + zeta * ( 27.0 - zeta*( 9.0 + 27.0*zeta ) ) )*0.0625;
            real c3 = ( -zeta*( 1.0 - 9.0 * zeta * ( 1.0 + zeta ) ) - 1.0 )*0.0625;

            real da0 = (   1.0 + xi*( 18.0 - 27.0*xi )) * 0.0625;
            real da1 = ( -27.0 - xi*( 18.0 - 81.0*xi )) * 0.0625;
            real da2 = (  27.0 - xi*( 18.0 + 81.0*xi )) * 0.0625;
            real da3 = (  -1.0 + xi*( 18.0 + 27.0*xi )) * 0.0625;

            real dda0 = ( 18.0 - 54.0*xi ) * 0.0625;
            real dda1 = ( 162.0*xi - 18.0 ) * 0.0625;
            real dda2 = ( - 162.0*xi - 18.0 ) * 0.0625;
            real dda3 = ( 54.0*xi + 18.0 ) * 0.0625;

            real db0 = (   1.0 + eta*( 18.0 - 27.0*eta )) * 0.0625;
            real db1 = ( -27.0 - eta*( 18.0 - 81.0*eta )) * 0.0625;
            real db2 = (  27.0 - eta*( 18.0 + 81.0*eta )) * 0.0625;
            real db3 = (  -1.0 + eta*( 18.0 + 27.0*eta )) * 0.0625;

            real ddb0 = ( 18.0 - 54.0*eta ) * 0.0625;
            real ddb1 = ( 162.0*eta - 18.0 ) * 0.0625;
            real ddb2 = ( - 162.0*eta - 18.0 ) * 0.0625;
            real ddb3 = ( 54.0*eta + 18.0 ) * 0.0625;

            real dc0 = (   1.0 + zeta*( 18.0 - 27.0*zeta )) * 0.0625;
            real dc1 = ( -27.0 - zeta*( 18.0 - 81.0*zeta )) * 0.0625;
            real dc2 = (  27.0 - zeta*( 18.0 + 81.0*zeta )) * 0.0625;
            real dc3 = (  -1.0 + zeta*( 18.0 + 27.0*zeta )) * 0.0625;

            real ddc0 = ( 18.0 - 54.0*zeta ) * 0.0625;
            real ddc1 = ( 162.0*zeta - 18.0 ) * 0.0625;
            real ddc2 = ( - 162.0*zeta - 18.0 ) * 0.0625;
            real ddc3 = ( 54.0*zeta + 18.0 ) * 0.0625;

            Matrix< DDRMat > td2NdXi2(6,64);
            td2NdXi2( 0,  0 ) = b0*c0*dda0;
            td2NdXi2( 1,  0 ) = a0*c0*ddb0;
            td2NdXi2( 2,  0 ) = a0*b0*ddc0;
            td2NdXi2( 3,  0 ) = a0*db0*dc0;
            td2NdXi2( 4,  0 ) = b0*da0*dc0;
            td2NdXi2( 5,  0 ) = c0*da0*db0;

            td2NdXi2( 0,  1 ) = b0*c0*dda3;
            td2NdXi2( 1,  1 ) = a3*c0*ddb0;
            td2NdXi2( 2,  1 ) = a3*b0*ddc0;
            td2NdXi2( 3,  1 ) = a3*db0*dc0;
            td2NdXi2( 4,  1 ) = b0*da3*dc0;
            td2NdXi2( 5,  1 ) = c0*da3*db0;

            td2NdXi2( 0,  2 ) = b3*c0*dda3;
            td2NdXi2( 1,  2 ) = a3*c0*ddb3;
            td2NdXi2( 2,  2 ) = a3*b3*ddc0;
            td2NdXi2( 3,  2 ) = a3*db3*dc0;
            td2NdXi2( 4,  2 ) = b3*da3*dc0;
            td2NdXi2( 5,  2 ) = c0*da3*db3;

            td2NdXi2( 0,  3 ) = b3*c0*dda0;
            td2NdXi2( 1,  3 ) = a0*c0*ddb3;
            td2NdXi2( 2,  3 ) = a0*b3*ddc0;
            td2NdXi2( 3,  3 ) = a0*db3*dc0;
            td2NdXi2( 4,  3 ) = b3*da0*dc0;
            td2NdXi2( 5,  3 ) = c0*da0*db3;

            td2NdXi2( 0,  4 ) = b0*c3*dda0;
            td2NdXi2( 1,  4 ) = a0*c3*ddb0;
            td2NdXi2( 2,  4 ) = a0*b0*ddc3;
            td2NdXi2( 3,  4 ) = a0*db0*dc3;
            td2NdXi2( 4,  4 ) = b0*da0*dc3;
            td2NdXi2( 5,  4 ) = c3*da0*db0;

            td2NdXi2( 0,  5 ) = b0*c3*dda3;
            td2NdXi2( 1,  5 ) = a3*c3*ddb0;
            td2NdXi2( 2,  5 ) = a3*b0*ddc3;
            td2NdXi2( 3,  5 ) = a3*db0*dc3;
            td2NdXi2( 4,  5 ) = b0*da3*dc3;
            td2NdXi2( 5,  5 ) = c3*da3*db0;

            td2NdXi2( 0,  6 ) = b3*c3*dda3;
            td2NdXi2( 1,  6 ) = a3*c3*ddb3;
            td2NdXi2( 2,  6 ) = a3*b3*ddc3;
            td2NdXi2( 3,  6 ) = a3*db3*dc3;
            td2NdXi2( 4,  6 ) = b3*da3*dc3;
            td2NdXi2( 5,  6 ) = c3*da3*db3;

            td2NdXi2( 0,  7 ) = b3*c3*dda0;
            td2NdXi2( 1,  7 ) = a0*c3*ddb3;
            td2NdXi2( 2,  7 ) = a0*b3*ddc3;
            td2NdXi2( 3,  7 ) = a0*db3*dc3;
            td2NdXi2( 4,  7 ) = b3*da0*dc3;
            td2NdXi2( 5,  7 ) = c3*da0*db3;

            td2NdXi2( 0,  8 ) = b0*c0*dda1;
            td2NdXi2( 1,  8 ) = a1*c0*ddb0;
            td2NdXi2( 2,  8 ) = a1*b0*ddc0;
            td2NdXi2( 3,  8 ) = a1*db0*dc0;
            td2NdXi2( 4,  8 ) = b0*da1*dc0;
            td2NdXi2( 5,  8 ) = c0*da1*db0;

            td2NdXi2( 0,  9 ) = b0*c0*dda2;
            td2NdXi2( 1,  9 ) = a2*c0*ddb0;
            td2NdXi2( 2,  9 ) = a2*b0*ddc0;
            td2NdXi2( 3,  9 ) = a2*db0*dc0;
            td2NdXi2( 4,  9 ) = b0*da2*dc0;
            td2NdXi2( 5,  9 ) = c0*da2*db0;

            td2NdXi2( 0, 10 ) = b1*c0*dda0;
            td2NdXi2( 1, 10 ) = a0*c0*ddb1;
            td2NdXi2( 2, 10 ) = a0*b1*ddc0;
            td2NdXi2( 3, 10 ) = a0*db1*dc0;
            td2NdXi2( 4, 10 ) = b1*da0*dc0;
            td2NdXi2( 5, 10 ) = c0*da0*db1;

            td2NdXi2( 0, 11 ) = b2*c0*dda0;
            td2NdXi2( 1, 11 ) = a0*c0*ddb2;
            td2NdXi2( 2, 11 ) = a0*b2*ddc0;
            td2NdXi2( 3, 11 ) = a0*db2*dc0;
            td2NdXi2( 4, 11 ) = b2*da0*dc0;
            td2NdXi2( 5, 11 ) = c0*da0*db2;

            td2NdXi2( 0, 12 ) = b0*c1*dda0;
            td2NdXi2( 1, 12 ) = a0*c1*ddb0;
            td2NdXi2( 2, 12 ) = a0*b0*ddc1;
            td2NdXi2( 3, 12 ) = a0*db0*dc1;
            td2NdXi2( 4, 12 ) = b0*da0*dc1;
            td2NdXi2( 5, 12 ) = c1*da0*db0;

            td2NdXi2( 0, 13 ) = b0*c2*dda0;
            td2NdXi2( 1, 13 ) = a0*c2*ddb0;
            td2NdXi2( 2, 13 ) = a0*b0*ddc2;
            td2NdXi2( 3, 13 ) = a0*db0*dc2;
            td2NdXi2( 4, 13 ) = b0*da0*dc2;
            td2NdXi2( 5, 13 ) = c2*da0*db0;

            td2NdXi2( 0, 14 ) = b1*c0*dda3;
            td2NdXi2( 1, 14 ) = a3*c0*ddb1;
            td2NdXi2( 2, 14 ) = a3*b1*ddc0;
            td2NdXi2( 3, 14 ) = a3*db1*dc0;
            td2NdXi2( 4, 14 ) = b1*da3*dc0;
            td2NdXi2( 5, 14 ) = c0*da3*db1;

            td2NdXi2( 0, 15 ) = b2*c0*dda3;
            td2NdXi2( 1, 15 ) = a3*c0*ddb2;
            td2NdXi2( 2, 15 ) = a3*b2*ddc0;
            td2NdXi2( 3, 15 ) = a3*db2*dc0;
            td2NdXi2( 4, 15 ) = b2*da3*dc0;
            td2NdXi2( 5, 15 ) = c0*da3*db2;

            td2NdXi2( 0, 16 ) = b0*c1*dda3;
            td2NdXi2( 1, 16 ) = a3*c1*ddb0;
            td2NdXi2( 2, 16 ) = a3*b0*ddc1;
            td2NdXi2( 3, 16 ) = a3*db0*dc1;
            td2NdXi2( 4, 16 ) = b0*da3*dc1;
            td2NdXi2( 5, 16 ) = c1*da3*db0;

            td2NdXi2( 0, 17 ) = b0*c2*dda3;
            td2NdXi2( 1, 17 ) = a3*c2*ddb0;
            td2NdXi2( 2, 17 ) = a3*b0*ddc2;
            td2NdXi2( 3, 17 ) = a3*db0*dc2;
            td2NdXi2( 4, 17 ) = b0*da3*dc2;
            td2NdXi2( 5, 17 ) = c2*da3*db0;

            td2NdXi2( 0, 18 ) = b3*c0*dda2;
            td2NdXi2( 1, 18 ) = a2*c0*ddb3;
            td2NdXi2( 2, 18 ) = a2*b3*ddc0;
            td2NdXi2( 3, 18 ) = a2*db3*dc0;
            td2NdXi2( 4, 18 ) = b3*da2*dc0;
            td2NdXi2( 5, 18 ) = c0*da2*db3;

            td2NdXi2( 0, 19 ) = b3*c0*dda1;
            td2NdXi2( 1, 19 ) = a1*c0*ddb3;
            td2NdXi2( 2, 19 ) = a1*b3*ddc0;
            td2NdXi2( 3, 19 ) = a1*db3*dc0;
            td2NdXi2( 4, 19 ) = b3*da1*dc0;
            td2NdXi2( 5, 19 ) = c0*da1*db3;

            td2NdXi2( 0, 20 ) = b3*c1*dda3;
            td2NdXi2( 1, 20 ) = a3*c1*ddb3;
            td2NdXi2( 2, 20 ) = a3*b3*ddc1;
            td2NdXi2( 3, 20 ) = a3*db3*dc1;
            td2NdXi2( 4, 20 ) = b3*da3*dc1;
            td2NdXi2( 5, 20 ) = c1*da3*db3;

            td2NdXi2( 0, 21 ) = b3*c2*dda3;
            td2NdXi2( 1, 21 ) = a3*c2*ddb3;
            td2NdXi2( 2, 21 ) = a3*b3*ddc2;
            td2NdXi2( 3, 21 ) = a3*db3*dc2;
            td2NdXi2( 4, 21 ) = b3*da3*dc2;
            td2NdXi2( 5, 21 ) = c2*da3*db3;

            td2NdXi2( 0, 22 ) = b3*c1*dda0;
            td2NdXi2( 1, 22 ) = a0*c1*ddb3;
            td2NdXi2( 2, 22 ) = a0*b3*ddc1;
            td2NdXi2( 3, 22 ) = a0*db3*dc1;
            td2NdXi2( 4, 22 ) = b3*da0*dc1;
            td2NdXi2( 5, 22 ) = c1*da0*db3;

            td2NdXi2( 0, 23 ) = b3*c2*dda0;
            td2NdXi2( 1, 23 ) = a0*c2*ddb3;
            td2NdXi2( 2, 23 ) = a0*b3*ddc2;
            td2NdXi2( 3, 23 ) = a0*db3*dc2;
            td2NdXi2( 4, 23 ) = b3*da0*dc2;
            td2NdXi2( 5, 23 ) = c2*da0*db3;

            td2NdXi2( 0, 24 ) = b0*c3*dda1;
            td2NdXi2( 1, 24 ) = a1*c3*ddb0;
            td2NdXi2( 2, 24 ) = a1*b0*ddc3;
            td2NdXi2( 3, 24 ) = a1*db0*dc3;
            td2NdXi2( 4, 24 ) = b0*da1*dc3;
            td2NdXi2( 5, 24 ) = c3*da1*db0;

            td2NdXi2( 0, 25 ) = b0*c3*dda2;
            td2NdXi2( 1, 25 ) = a2*c3*ddb0;
            td2NdXi2( 2, 25 ) = a2*b0*ddc3;
            td2NdXi2( 3, 25 ) = a2*db0*dc3;
            td2NdXi2( 4, 25 ) = b0*da2*dc3;
            td2NdXi2( 5, 25 ) = c3*da2*db0;

            td2NdXi2( 0, 26 ) = b1*c3*dda0;
            td2NdXi2( 1, 26 ) = a0*c3*ddb1;
            td2NdXi2( 2, 26 ) = a0*b1*ddc3;
            td2NdXi2( 3, 26 ) = a0*db1*dc3;
            td2NdXi2( 4, 26 ) = b1*da0*dc3;
            td2NdXi2( 5, 26 ) = c3*da0*db1;

            td2NdXi2( 0, 27 ) = b2*c3*dda0;
            td2NdXi2( 1, 27 ) = a0*c3*ddb2;
            td2NdXi2( 2, 27 ) = a0*b2*ddc3;
            td2NdXi2( 3, 27 ) = a0*db2*dc3;
            td2NdXi2( 4, 27 ) = b2*da0*dc3;
            td2NdXi2( 5, 27 ) = c3*da0*db2;

            td2NdXi2( 0, 28 ) = b1*c3*dda3;
            td2NdXi2( 1, 28 ) = a3*c3*ddb1;
            td2NdXi2( 2, 28 ) = a3*b1*ddc3;
            td2NdXi2( 3, 28 ) = a3*db1*dc3;
            td2NdXi2( 4, 28 ) = b1*da3*dc3;
            td2NdXi2( 5, 28 ) = c3*da3*db1;

            td2NdXi2( 0, 29 ) = b2*c3*dda3;
            td2NdXi2( 1, 29 ) = a3*c3*ddb2;
            td2NdXi2( 2, 29 ) = a3*b2*ddc3;
            td2NdXi2( 3, 29 ) = a3*db2*dc3;
            td2NdXi2( 4, 29 ) = b2*da3*dc3;
            td2NdXi2( 5, 29 ) = c3*da3*db2;

            td2NdXi2( 0, 30 ) = b3*c3*dda2;
            td2NdXi2( 1, 30 ) = a2*c3*ddb3;
            td2NdXi2( 2, 30 ) = a2*b3*ddc3;
            td2NdXi2( 3, 30 ) = a2*db3*dc3;
            td2NdXi2( 4, 30 ) = b3*da2*dc3;
            td2NdXi2( 5, 30 ) = c3*da2*db3;

            td2NdXi2( 0, 31 ) = b3*c3*dda1;
            td2NdXi2( 1, 31 ) = a1*c3*ddb3;
            td2NdXi2( 2, 31 ) = a1*b3*ddc3;
            td2NdXi2( 3, 31 ) = a1*db3*dc3;
            td2NdXi2( 4, 31 ) = b3*da1*dc3;
            td2NdXi2( 5, 31 ) = c3*da1*db3;

            td2NdXi2( 0, 32 ) = b1*c0*dda1;
            td2NdXi2( 1, 32 ) = a1*c0*ddb1;
            td2NdXi2( 2, 32 ) = a1*b1*ddc0;
            td2NdXi2( 3, 32 ) = a1*db1*dc0;
            td2NdXi2( 4, 32 ) = b1*da1*dc0;
            td2NdXi2( 5, 32 ) = c0*da1*db1;

            td2NdXi2( 0, 33 ) = b2*c0*dda1;
            td2NdXi2( 1, 33 ) = a1*c0*ddb2;
            td2NdXi2( 2, 33 ) = a1*b2*ddc0;
            td2NdXi2( 3, 33 ) = a1*db2*dc0;
            td2NdXi2( 4, 33 ) = b2*da1*dc0;
            td2NdXi2( 5, 33 ) = c0*da1*db2;

            td2NdXi2( 0, 34 ) = b2*c0*dda2;
            td2NdXi2( 1, 34 ) = a2*c0*ddb2;
            td2NdXi2( 2, 34 ) = a2*b2*ddc0;
            td2NdXi2( 3, 34 ) = a2*db2*dc0;
            td2NdXi2( 4, 34 ) = b2*da2*dc0;
            td2NdXi2( 5, 34 ) = c0*da2*db2;

            td2NdXi2( 0, 35 ) = b1*c0*dda2;
            td2NdXi2( 1, 35 ) = a2*c0*ddb1;
            td2NdXi2( 2, 35 ) = a2*b1*ddc0;
            td2NdXi2( 3, 35 ) = a2*db1*dc0;
            td2NdXi2( 4, 35 ) = b1*da2*dc0;
            td2NdXi2( 5, 35 ) = c0*da2*db1;

            td2NdXi2( 0, 36 ) = b0*c1*dda1;
            td2NdXi2( 1, 36 ) = a1*c1*ddb0;
            td2NdXi2( 2, 36 ) = a1*b0*ddc1;
            td2NdXi2( 3, 36 ) = a1*db0*dc1;
            td2NdXi2( 4, 36 ) = b0*da1*dc1;
            td2NdXi2( 5, 36 ) = c1*da1*db0;

            td2NdXi2( 0, 37 ) = b0*c1*dda2;
            td2NdXi2( 1, 37 ) = a2*c1*ddb0;
            td2NdXi2( 2, 37 ) = a2*b0*ddc1;
            td2NdXi2( 3, 37 ) = a2*db0*dc1;
            td2NdXi2( 4, 37 ) = b0*da2*dc1;
            td2NdXi2( 5, 37 ) = c1*da2*db0;

            td2NdXi2( 0, 38 ) = b0*c2*dda2;
            td2NdXi2( 1, 38 ) = a2*c2*ddb0;
            td2NdXi2( 2, 38 ) = a2*b0*ddc2;
            td2NdXi2( 3, 38 ) = a2*db0*dc2;
            td2NdXi2( 4, 38 ) = b0*da2*dc2;
            td2NdXi2( 5, 38 ) = c2*da2*db0;

            td2NdXi2( 0, 39 ) = b0*c2*dda1;
            td2NdXi2( 1, 39 ) = a1*c2*ddb0;
            td2NdXi2( 2, 39 ) = a1*b0*ddc2;
            td2NdXi2( 3, 39 ) = a1*db0*dc2;
            td2NdXi2( 4, 39 ) = b0*da1*dc2;
            td2NdXi2( 5, 39 ) = c2*da1*db0;

            td2NdXi2( 0, 40 ) = b1*c1*dda0;
            td2NdXi2( 1, 40 ) = a0*c1*ddb1;
            td2NdXi2( 2, 40 ) = a0*b1*ddc1;
            td2NdXi2( 3, 40 ) = a0*db1*dc1;
            td2NdXi2( 4, 40 ) = b1*da0*dc1;
            td2NdXi2( 5, 40 ) = c1*da0*db1;

            td2NdXi2( 0, 41 ) = b1*c2*dda0;
            td2NdXi2( 1, 41 ) = a0*c2*ddb1;
            td2NdXi2( 2, 41 ) = a0*b1*ddc2;
            td2NdXi2( 3, 41 ) = a0*db1*dc2;
            td2NdXi2( 4, 41 ) = b1*da0*dc2;
            td2NdXi2( 5, 41 ) = c2*da0*db1;

            td2NdXi2( 0, 42 ) = b2*c2*dda0;
            td2NdXi2( 1, 42 ) = a0*c2*ddb2;
            td2NdXi2( 2, 42 ) = a0*b2*ddc2;
            td2NdXi2( 3, 42 ) = a0*db2*dc2;
            td2NdXi2( 4, 42 ) = b2*da0*dc2;
            td2NdXi2( 5, 42 ) = c2*da0*db2;

            td2NdXi2( 0, 43 ) = b2*c1*dda0;
            td2NdXi2( 1, 43 ) = a0*c1*ddb2;
            td2NdXi2( 2, 43 ) = a0*b2*ddc1;
            td2NdXi2( 3, 43 ) = a0*db2*dc1;
            td2NdXi2( 4, 43 ) = b2*da0*dc1;
            td2NdXi2( 5, 43 ) = c1*da0*db2;

            td2NdXi2( 0, 44 ) = b1*c1*dda3;
            td2NdXi2( 1, 44 ) = a3*c1*ddb1;
            td2NdXi2( 2, 44 ) = a3*b1*ddc1;
            td2NdXi2( 3, 44 ) = a3*db1*dc1;
            td2NdXi2( 4, 44 ) = b1*da3*dc1;
            td2NdXi2( 5, 44 ) = c1*da3*db1;

            td2NdXi2( 0, 45 ) = b2*c1*dda3;
            td2NdXi2( 1, 45 ) = a3*c1*ddb2;
            td2NdXi2( 2, 45 ) = a3*b2*ddc1;
            td2NdXi2( 3, 45 ) = a3*db2*dc1;
            td2NdXi2( 4, 45 ) = b2*da3*dc1;
            td2NdXi2( 5, 45 ) = c1*da3*db2;

            td2NdXi2( 0, 46 ) = b2*c2*dda3;
            td2NdXi2( 1, 46 ) = a3*c2*ddb2;
            td2NdXi2( 2, 46 ) = a3*b2*ddc2;
            td2NdXi2( 3, 46 ) = a3*db2*dc2;
            td2NdXi2( 4, 46 ) = b2*da3*dc2;
            td2NdXi2( 5, 46 ) = c2*da3*db2;

            td2NdXi2( 0, 47 ) = b1*c2*dda3;
            td2NdXi2( 1, 47 ) = a3*c2*ddb1;
            td2NdXi2( 2, 47 ) = a3*b1*ddc2;
            td2NdXi2( 3, 47 ) = a3*db1*dc2;
            td2NdXi2( 4, 47 ) = b1*da3*dc2;
            td2NdXi2( 5, 47 ) = c2*da3*db1;

            td2NdXi2( 0, 48 ) = b3*c1*dda2;
            td2NdXi2( 1, 48 ) = a2*c1*ddb3;
            td2NdXi2( 2, 48 ) = a2*b3*ddc1;
            td2NdXi2( 3, 48 ) = a2*db3*dc1;
            td2NdXi2( 4, 48 ) = b3*da2*dc1;
            td2NdXi2( 5, 48 ) = c1*da2*db3;

            td2NdXi2( 0, 49 ) = b3*c1*dda1;
            td2NdXi2( 1, 49 ) = a1*c1*ddb3;
            td2NdXi2( 2, 49 ) = a1*b3*ddc1;
            td2NdXi2( 3, 49 ) = a1*db3*dc1;
            td2NdXi2( 4, 49 ) = b3*da1*dc1;
            td2NdXi2( 5, 49 ) = c1*da1*db3;

            td2NdXi2( 0, 50 ) = b3*c2*dda1;
            td2NdXi2( 1, 50 ) = a1*c2*ddb3;
            td2NdXi2( 2, 50 ) = a1*b3*ddc2;
            td2NdXi2( 3, 50 ) = a1*db3*dc2;
            td2NdXi2( 4, 50 ) = b3*da1*dc2;
            td2NdXi2( 5, 50 ) = c2*da1*db3;

            td2NdXi2( 0, 51 ) = b3*c2*dda2;
            td2NdXi2( 1, 51 ) = a2*c2*ddb3;
            td2NdXi2( 2, 51 ) = a2*b3*ddc2;
            td2NdXi2( 3, 51 ) = a2*db3*dc2;
            td2NdXi2( 4, 51 ) = b3*da2*dc2;
            td2NdXi2( 5, 51 ) = c2*da2*db3;

            td2NdXi2( 0, 52 ) = b1*c3*dda1;
            td2NdXi2( 1, 52 ) = a1*c3*ddb1;
            td2NdXi2( 2, 52 ) = a1*b1*ddc3;
            td2NdXi2( 3, 52 ) = a1*db1*dc3;
            td2NdXi2( 4, 52 ) = b1*da1*dc3;
            td2NdXi2( 5, 52 ) = c3*da1*db1;

            td2NdXi2( 0, 53 ) = b1*c3*dda2;
            td2NdXi2( 1, 53 ) = a2*c3*ddb1;
            td2NdXi2( 2, 53 ) = a2*b1*ddc3;
            td2NdXi2( 3, 53 ) = a2*db1*dc3;
            td2NdXi2( 4, 53 ) = b1*da2*dc3;
            td2NdXi2( 5, 53 ) = c3*da2*db1;

            td2NdXi2( 0, 54 ) = b2*c3*dda2;
            td2NdXi2( 1, 54 ) = a2*c3*ddb2;
            td2NdXi2( 2, 54 ) = a2*b2*ddc3;
            td2NdXi2( 3, 54 ) = a2*db2*dc3;
            td2NdXi2( 4, 54 ) = b2*da2*dc3;
            td2NdXi2( 5, 54 ) = c3*da2*db2;

            td2NdXi2( 0, 55 ) = b2*c3*dda1;
            td2NdXi2( 1, 55 ) = a1*c3*ddb2;
            td2NdXi2( 2, 55 ) = a1*b2*ddc3;
            td2NdXi2( 3, 55 ) = a1*db2*dc3;
            td2NdXi2( 4, 55 ) = b2*da1*dc3;
            td2NdXi2( 5, 55 ) = c3*da1*db2;

            td2NdXi2( 0, 56 ) = b1*c1*dda1;
            td2NdXi2( 1, 56 ) = a1*c1*ddb1;
            td2NdXi2( 2, 56 ) = a1*b1*ddc1;
            td2NdXi2( 3, 56 ) = a1*db1*dc1;
            td2NdXi2( 4, 56 ) = b1*da1*dc1;
            td2NdXi2( 5, 56 ) = c1*da1*db1;

            td2NdXi2( 0, 57 ) = b1*c1*dda2;
            td2NdXi2( 1, 57 ) = a2*c1*ddb1;
            td2NdXi2( 2, 57 ) = a2*b1*ddc1;
            td2NdXi2( 3, 57 ) = a2*db1*dc1;
            td2NdXi2( 4, 57 ) = b1*da2*dc1;
            td2NdXi2( 5, 57 ) = c1*da2*db1;

            td2NdXi2( 0, 58 ) = b2*c1*dda2;
            td2NdXi2( 1, 58 ) = a2*c1*ddb2;
            td2NdXi2( 2, 58 ) = a2*b2*ddc1;
            td2NdXi2( 3, 58 ) = a2*db2*dc1;
            td2NdXi2( 4, 58 ) = b2*da2*dc1;
            td2NdXi2( 5, 58 ) = c1*da2*db2;

            td2NdXi2( 0, 59 ) = b2*c1*dda1;
            td2NdXi2( 1, 59 ) = a1*c1*ddb2;
            td2NdXi2( 2, 59 ) = a1*b2*ddc1;
            td2NdXi2( 3, 59 ) = a1*db2*dc1;
            td2NdXi2( 4, 59 ) = b2*da1*dc1;
            td2NdXi2( 5, 59 ) = c1*da1*db2;

            td2NdXi2( 0, 60 ) = b1*c2*dda1;
            td2NdXi2( 1, 60 ) = a1*c2*ddb1;
            td2NdXi2( 2, 60 ) = a1*b1*ddc2;
            td2NdXi2( 3, 60 ) = a1*db1*dc2;
            td2NdXi2( 4, 60 ) = b1*da1*dc2;
            td2NdXi2( 5, 60 ) = c2*da1*db1;

            td2NdXi2( 0, 61 ) = b1*c2*dda2;
            td2NdXi2( 1, 61 ) = a2*c2*ddb1;
            td2NdXi2( 2, 61 ) = a2*b1*ddc2;
            td2NdXi2( 3, 61 ) = a2*db1*dc2;
            td2NdXi2( 4, 61 ) = b1*da2*dc2;
            td2NdXi2( 5, 61 ) = c2*da2*db1;

            td2NdXi2( 0, 62 ) = b2*c2*dda2;
            td2NdXi2( 1, 62 ) = a2*c2*ddb2;
            td2NdXi2( 2, 62 ) = a2*b2*ddc2;
            td2NdXi2( 3, 62 ) = a2*db2*dc2;
            td2NdXi2( 4, 62 ) = b2*da2*dc2;
            td2NdXi2( 5, 62 ) = c2*da2*db2;

            td2NdXi2( 0, 63 ) = b2*c2*dda1;
            td2NdXi2( 1, 63 ) = a1*c2*ddb2;
            td2NdXi2( 2, 63 ) = a1*b2*ddc2;
            td2NdXi2( 3, 63 ) = a1*db2*dc2;
            td2NdXi2( 4, 63 ) = b2*da1*dc2;
            td2NdXi2( 5, 63 ) = c2*da1*db2;
            return td2NdXi2;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX64_HPP_ */
