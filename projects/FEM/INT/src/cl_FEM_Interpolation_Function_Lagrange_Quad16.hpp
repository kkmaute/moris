/*
 * cl_FEM_Interpolation_Function_Quad16.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD16_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD16_HPP_

#include "assert.h"

#include "cl_FEM_Interpolation_Matrix.hpp"
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
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 16 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 16  >::get_param_coords(
                Mat<real> & aXihat ) const
        {
            aXihat.set_size( 2, 16 );

            real c = 1.0/3.0;

            aXihat( 0,  0 ) = -1.000000;
            aXihat( 1,  0 ) = -1.000000;
            aXihat( 0,  1 ) =  1.000000;
            aXihat( 1,  1 ) = -1.000000;
            aXihat( 0,  2 ) =  1.000000;
            aXihat( 1,  2 ) =  1.000000;
            aXihat( 0,  3 ) = -1.000000;
            aXihat( 1,  3 ) =  1.000000;
            aXihat( 0,  4 ) = -c;
            aXihat( 1,  4 ) = -1.000000;
            aXihat( 0,  5 ) =  c;
            aXihat( 1,  5 ) = -1.000000;
            aXihat( 0,  6 ) =  1.000000;
            aXihat( 1,  6 ) = -c;
            aXihat( 0,  7 ) =  1.000000;
            aXihat( 1,  7 ) =  c;
            aXihat( 0,  8 ) =  c;
            aXihat( 1,  8 ) =  1.000000;
            aXihat( 0,  9 ) = -c;
            aXihat( 1,  9 ) =  1.000000;
            aXihat( 0, 10 ) = -1.000000;
            aXihat( 1, 10 ) =  c;
            aXihat( 0, 11 ) = -1.000000;
            aXihat( 1, 11 ) = -c;
            aXihat( 0, 12 ) = -c;
            aXihat( 1, 12 ) = -c;
            aXihat( 0, 13 ) =  c;
            aXihat( 1, 13 ) = -c;
            aXihat( 0, 14 ) =  c;
            aXihat( 1, 14 ) =  c;
            aXihat( 0, 15 ) = -c;
            aXihat( 1, 15 ) =  c;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 16  >::eval_N(
                  Interpolation_Matrix  & aN,
            const Mat<real> & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( aN.n_cols() == 16,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( aN.n_rows() == 1,
                    "eval_shape: aN not allocated or hat wrong size." );


            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;


            // populate matrix with values
            aN(  0 ) = a0*b0;
            aN(  1 ) = a3*b0;
            aN(  2 ) = a3*b3;
            aN(  3 ) = a0*b3;
            aN(  4 ) = a1*b0;
            aN(  5 ) = a2*b0;
            aN(  6 ) = a3*b1;
            aN(  7 ) = a3*b2;
            aN(  8 ) = a2*b3;
            aN(  9 ) = a1*b3;
            aN( 10 ) = a0*b2;
            aN( 11 ) = a0*b1;
            aN( 12 ) = a1*b1;
            aN( 13 ) = a2*b1;
            aN( 14 ) = a2*b2;
            aN( 15 ) = a1*b2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 16  >::eval_dNdXi(
                       Interpolation_Matrix & adNdXi,
                const Mat<real> & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( adNdXi.n_cols() == 16,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( adNdXi.n_rows() == 2,
                    "eval_shape: aN not allocated or hat wrong size." );


            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

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
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 16  >::eval_d2NdXi2(
                      Interpolation_Matrix  & ad2NdXi2,
                const Mat<real> & aXi ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( ad2NdXi2.n_cols() == 16,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( ad2NdXi2.n_rows() == 3,
                    "eval_shape: aN not allocated or hat wrong size." );


            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

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
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD16_HPP_ */
