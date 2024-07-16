/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Quad16.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD16_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD16_HPP_

#include "assert.h"
#include "moris_typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_MTK_Interpolation_Function.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 16 >::get_number_of_param_dimensions() const
        {
            return 2;
        }

//------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 16 >::get_interpolation_order() const
        {
            return Interpolation_Order::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 16 >::get_param_coords( Matrix< DDRMat > & aXiHat ) const
        {
            aXiHat.set_size( 2, 16, 0.0 );
            real c = 1.0/3.0;
            aXiHat( 0,  0 ) = -1.000000;
            aXiHat( 1,  0 ) = -1.000000;
            aXiHat( 0,  1 ) =  1.000000;
            aXiHat( 1,  1 ) = -1.000000;
            aXiHat( 0,  2 ) =  1.000000;
            aXiHat( 1,  2 ) =  1.000000;
            aXiHat( 0,  3 ) = -1.000000;
            aXiHat( 1,  3 ) =  1.000000;
            aXiHat( 0,  4 ) = -c;
            aXiHat( 1,  4 ) = -1.000000;
            aXiHat( 0,  5 ) =  c;
            aXiHat( 1,  5 ) = -1.000000;
            aXiHat( 0,  6 ) =  1.000000;
            aXiHat( 1,  6 ) = -c;
            aXiHat( 0,  7 ) =  1.000000;
            aXiHat( 1,  7 ) =  c;
            aXiHat( 0,  8 ) =  c;
            aXiHat( 1,  8 ) =  1.000000;
            aXiHat( 0,  9 ) = -c;
            aXiHat( 1,  9 ) =  1.000000;
            aXiHat( 0, 10 ) = -1.000000;
            aXiHat( 1, 10 ) =  c;
            aXiHat( 0, 11 ) = -1.000000;
            aXiHat( 1, 11 ) = -c;
            aXiHat( 0, 12 ) = -c;
            aXiHat( 1, 12 ) = -c;
            aXiHat( 0, 13 ) =  c;
            aXiHat( 1, 13 ) = -c;
            aXiHat( 0, 14 ) =  c;
            aXiHat( 1, 14 ) =  c;
            aXiHat( 0, 15 ) = -c;
            aXiHat( 1, 15 ) =  c;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 16 >::eval_N(const Matrix< DDRMat > & aXi,
                                                                                                              Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD16 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

            // populate matrix with values
            aNXi.set_size(1,16);
            aNXi(  0 ) = a0*b0;
            aNXi(  1 ) = a3*b0;
            aNXi(  2 ) = a3*b3;
            aNXi(  3 ) = a0*b3;
            aNXi(  4 ) = a1*b0;
            aNXi(  5 ) = a2*b0;
            aNXi(  6 ) = a3*b1;
            aNXi(  7 ) = a3*b2;
            aNXi(  8 ) = a2*b3;
            aNXi(  9 ) = a1*b3;
            aNXi( 10 ) = a0*b2;
            aNXi( 11 ) = a0*b1;
            aNXi( 12 ) = a1*b1;
            aNXi( 13 ) = a2*b1;
            aNXi( 14 ) = a2*b2;
            aNXi( 15 ) = a1*b2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 16 >::eval_dNdXi( const Matrix< DDRMat > & aXi,
                                                                                                                   Matrix< DDRMat > & adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD16 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            // often used parameters
            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 ) * 0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) ) * 0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) ) * 0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 ) * 0.0625;

            real da0 = (   1.0 + xi*( 18.0 - 27.0*xi )) * 0.0625;
            real da1 = ( -27.0 - xi*( 18.0 - 81.0*xi )) * 0.0625;
            real da2 = (  27.0 - xi*( 18.0 + 81.0*xi )) * 0.0625;
            real da3 = (  -1.0 + xi*( 18.0 + 27.0*xi )) * 0.0625;

            real db0 = (   1.0 + eta*( 18.0 - 27.0*eta )) * 0.0625;
            real db1 = ( -27.0 - eta*( 18.0 - 81.0*eta )) * 0.0625;
            real db2 = (  27.0 - eta*( 18.0 + 81.0*eta )) * 0.0625;
            real db3 = (  -1.0 + eta*( 18.0 + 27.0*eta )) * 0.0625;

            // populate output matrix
            adNdXi.set_size( 2, 16 );
            adNdXi( 0,  0 ) = da0*b0;
            adNdXi( 1,  0 ) = a0*db0;

            adNdXi( 0,  1 ) = da3*b0;
            adNdXi( 1,  1 ) = a3*db0;

            adNdXi( 0,  2 ) = da3*b3;
            adNdXi( 1,  2 ) = a3*db3;

            adNdXi( 0,  3 ) = da0*b3;
            adNdXi( 1,  3 ) = a0*db3;

            adNdXi( 0,  4 ) = da1*b0;
            adNdXi( 1,  4 ) = a1*db0;

            adNdXi( 0,  5 ) = da2*b0;
            adNdXi( 1,  5 ) = a2*db0;

            adNdXi( 0,  6 ) = da3*b1;
            adNdXi( 1,  6 ) = a3*db1;

            adNdXi( 0,  7 ) = da3*b2;
            adNdXi( 1,  7 ) = a3*db2;

            adNdXi( 0,  8 ) = da2*b3;
            adNdXi( 1,  8 ) = a2*db3;

            adNdXi( 0,  9 ) = da1*b3;
            adNdXi( 1,  9 ) = a1*db3;

            adNdXi( 0, 10 ) = da0*b2;
            adNdXi( 1, 10 ) = a0*db2;

            adNdXi( 0, 11 ) = da0*b1;
            adNdXi( 1, 11 ) = a0*db1;

            adNdXi( 0, 12 ) = da1*b1;
            adNdXi( 1, 12 ) = a1*db1;

            adNdXi( 0, 13 ) = da2*b1;
            adNdXi( 1, 13 ) = a2*db1;

            adNdXi( 0, 14 ) = da2*b2;
            adNdXi( 1, 14 ) = a2*db2;

            adNdXi( 0, 15 ) = da1*b2;
            adNdXi( 1, 15 ) = a1*db2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 16 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi,
                                                                                                                     Matrix< DDRMat > & ad2NdXi2 ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD16 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            // often used parameters
            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 ) * 0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) ) * 0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) ) * 0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 ) * 0.0625;

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

            ad2NdXi2.set_size( 3, 16 );
            ad2NdXi2( 0,   0 ) = dda0*b0;
            ad2NdXi2( 1,   0 ) = a0*ddb0;
            ad2NdXi2( 2,   0 ) = da0*db0;

            ad2NdXi2( 0,   1 ) = dda3*b0;
            ad2NdXi2( 1,   1 ) = a3*ddb0;
            ad2NdXi2( 2,   1 ) = da3*db0;

            ad2NdXi2( 0,   2 ) = dda3*b3;
            ad2NdXi2( 1,   2 ) = a3*ddb3;
            ad2NdXi2( 2,   2 ) = da3*db3;

            ad2NdXi2( 0,   3 ) = dda0*b3;
            ad2NdXi2( 1,   3 ) = a0*ddb3;
            ad2NdXi2( 2,   3 ) = da0*db3;

            ad2NdXi2( 0,   4 ) = dda1*b0;
            ad2NdXi2( 1,   4 ) = a1*ddb0;
            ad2NdXi2( 2,   4 ) = da1*db0;

            ad2NdXi2( 0,   5 ) = dda2*b0;
            ad2NdXi2( 1,   5 ) = a2*ddb0;
            ad2NdXi2( 2,   5 ) = da2*db0;

            ad2NdXi2( 0,   6 ) = dda3*b1;
            ad2NdXi2( 1,   6 ) = a3*ddb1;
            ad2NdXi2( 2,   6 ) = da3*db1;

            ad2NdXi2( 0,   7 ) = dda3*b2;
            ad2NdXi2( 1,   7 ) = a3*ddb2;
            ad2NdXi2( 2,   7 ) = da3*db2;

            ad2NdXi2( 0,   8 ) = dda2*b3;
            ad2NdXi2( 1,   8 ) = a2*ddb3;
            ad2NdXi2( 2,   8 ) = da2*db3;

            ad2NdXi2( 0,   9 ) = dda1*b3;
            ad2NdXi2( 1,   9 ) = a1*ddb3;
            ad2NdXi2( 2,   9 ) = da1*db3;

            ad2NdXi2( 0,  10 ) = dda0*b2;
            ad2NdXi2( 1,  10 ) = a0*ddb2;
            ad2NdXi2( 2,  10 ) = da0*db2;

            ad2NdXi2( 0,  11 ) = dda0*b1;
            ad2NdXi2( 1,  11 ) = a0*ddb1;
            ad2NdXi2( 2,  11 ) = da0*db1;

            ad2NdXi2( 0,  12 ) = dda1*b1;
            ad2NdXi2( 1,  12 ) = a1*ddb1;
            ad2NdXi2( 2,  12 ) = da1*db1;

            ad2NdXi2( 0,  13 ) = dda2*b1;
            ad2NdXi2( 1,  13 ) = a2*ddb1;
            ad2NdXi2( 2,  13 ) = da2*db1;

            ad2NdXi2( 0,  14 ) = dda2*b2;
            ad2NdXi2( 1,  14 ) = a2*ddb2;
            ad2NdXi2( 2,  14 ) = da2*db2;

            ad2NdXi2( 0,  15 ) = dda1*b2;
            ad2NdXi2( 1,  15 ) = a1*ddb2;
            ad2NdXi2( 2,  15 ) = da1*db2;

        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 16 >::eval_d3NdXi3( const Matrix< DDRMat > & aXi,
                                                                                                                     Matrix< DDRMat > & ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD16 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            // often used parameters
            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 ) * 0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) ) * 0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) ) * 0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 ) * 0.0625;

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

            real ddda0 = -  54.0 * 0.0625;
            real ddda1 =   162.0 * 0.0625;
            real ddda2 = - 162.0 * 0.0625;
            real ddda3 =    54.0 * 0.0625;

            real dddb0 = -  54.0 * 0.0625;
            real dddb1 =   162.0 * 0.0625;
            real dddb2 = - 162.0 * 0.0625;
            real dddb3 =    54.0 * 0.0625;

            ad3NdXi3.set_size( 4, 16 );

            // 0th Node: (0,0)
            ad3NdXi3( 0,  0 ) = ddda0*   b0;
            ad3NdXi3( 1,  0 ) =    a0*dddb0;
            ad3NdXi3( 2,  0 ) =  dda0*  db0;
            ad3NdXi3( 3,  0 ) =   da0* ddb0;

            // 1th Node: (3,0)
            ad3NdXi3( 0,  1 ) = ddda3*   b0;
            ad3NdXi3( 1,  1 ) =    a3*dddb0;
            ad3NdXi3( 2,  1 ) =  dda3*  db0;
            ad3NdXi3( 3,  1 ) =   da3* ddb0;

            // 2th Node: (3,3)
            ad3NdXi3( 0,  2 ) = ddda3*   b3;
            ad3NdXi3( 1,  2 ) =    a3*dddb3;
            ad3NdXi3( 2,  2 ) =  dda3*  db3;
            ad3NdXi3( 3,  2 ) =   da3* ddb3;

            // 3th Node: (0,3)
            ad3NdXi3( 0,  3 ) = ddda0*   b3;
            ad3NdXi3( 1,  3 ) =    a0*dddb3;
            ad3NdXi3( 2,  3 ) =  dda0*  db3;
            ad3NdXi3( 3,  3 ) =   da0* ddb3;

            // 4th Node: (1,0)
            ad3NdXi3( 0,  4 ) = ddda1*   b0;
            ad3NdXi3( 1,  4 ) =    a1*dddb0;
            ad3NdXi3( 2,  4 ) =  dda1*  db0;
            ad3NdXi3( 3,  4 ) =   da1* ddb0;

            // 5th Node: (2,0)
            ad3NdXi3( 0,  5 ) = ddda2*   b0;
            ad3NdXi3( 1,  5 ) =    a2*dddb0;
            ad3NdXi3( 2,  5 ) =  dda2*  db0;
            ad3NdXi3( 3,  5 ) =   da2* ddb0;

            // 6th Node: (3,1)
            ad3NdXi3( 0,  6 ) = ddda3*   b1;
            ad3NdXi3( 1,  6 ) =    a3*dddb1;
            ad3NdXi3( 2,  6 ) =  dda3*  db1;
            ad3NdXi3( 3,  6 ) =   da3* ddb1;

            // 7th Node: (3,2)
            ad3NdXi3( 0,  7 ) = ddda3*   b2;
            ad3NdXi3( 1,  7 ) =    a3*dddb2;
            ad3NdXi3( 2,  7 ) =  dda3*  db2;
            ad3NdXi3( 3,  7 ) =   da3* ddb2;

            // 8th Node: (2,3)
            ad3NdXi3( 0,  8 ) = ddda2*   b3;
            ad3NdXi3( 1,  8 ) =    a2*dddb3;
            ad3NdXi3( 2,  8 ) =  dda2*  db3;
            ad3NdXi3( 3,  8 ) =   da2* ddb3;

            // 9th Node: (1,3)
            ad3NdXi3( 0,  9 ) = ddda1*   b3;
            ad3NdXi3( 1,  9 ) =    a1*dddb3;
            ad3NdXi3( 2,  9 ) =  dda1*  db3;
            ad3NdXi3( 3,  9 ) =   da1* ddb3;

            // 10th Node: (0,2)
            ad3NdXi3( 0, 10 ) = ddda0*   b2;
            ad3NdXi3( 1, 10 ) =    a0*dddb2;
            ad3NdXi3( 2, 10 ) =  dda0*  db2;
            ad3NdXi3( 3, 10 ) =   da0* ddb2;

            // 11th Node: (0,1)
            ad3NdXi3( 0, 11 ) = ddda0*   b1;
            ad3NdXi3( 1, 11 ) =    a0*dddb1;
            ad3NdXi3( 2, 11 ) =  dda0*  db1;
            ad3NdXi3( 3, 11 ) =   da0* ddb1;

            // 12th Node: (1,1)
            ad3NdXi3( 0, 12 ) = ddda1*   b1;
            ad3NdXi3( 1, 12 ) =    a1*dddb1;
            ad3NdXi3( 2, 12 ) =  dda1*  db1;
            ad3NdXi3( 3, 12 ) =   da1* ddb1;

            // 13th Node: (2,1)
            ad3NdXi3( 0, 13 ) = ddda2*   b1;
            ad3NdXi3( 1, 13 ) =    a2*dddb1;
            ad3NdXi3( 2, 13 ) =  dda2*  db1;
            ad3NdXi3( 3, 13 ) =   da2* ddb1;

            // 14th Node: (2,2)
            ad3NdXi3( 0, 14 ) = ddda2*   b2;
            ad3NdXi3( 1, 14 ) =    a2*dddb2;
            ad3NdXi3( 2, 14 ) =  dda2*  db2;
            ad3NdXi3( 3, 14 ) =   da2* ddb2;

            // 15th Node: (1,2)
            ad3NdXi3( 0, 15 ) = ddda1*   b2;
            ad3NdXi3( 1, 15 ) =    a1*dddb2;
            ad3NdXi3( 2, 15 ) =  dda1*  db2;
            ad3NdXi3( 3, 15 ) =   da1* ddb2;
        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD16_HPP_ */
