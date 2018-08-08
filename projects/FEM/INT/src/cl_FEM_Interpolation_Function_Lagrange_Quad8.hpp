/*
 * cl_FEM_Interpolation_Function_Quad4.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD8_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD8_HPP_

#include "assert.h"

#include "cl_FEM_Interpolation_Matrix.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src

// checked against femdoc:
// - N
// - adNdxi
// - d2Ndxi ( fixed bug in femdoc)
namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 8  >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::SERENDIPITY;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 8  >::get_param_coords(
                Mat<real> & aXihat ) const
        {
            aXihat.set_size( 2, 8 );

            aXihat( 0, 0 ) = -1.000000;
            aXihat( 1, 0 ) = -1.000000;
            aXihat( 0, 1 ) =  1.000000;
            aXihat( 1, 1 ) = -1.000000;
            aXihat( 0, 2 ) =  1.000000;
            aXihat( 1, 2 ) =  1.000000;
            aXihat( 0, 3 ) = -1.000000;
            aXihat( 1, 3 ) =  1.000000;
            aXihat( 0, 4 ) =  0.000000;
            aXihat( 1, 4 ) = -1.000000;
            aXihat( 0, 5 ) =  1.000000;
            aXihat( 1, 5 ) =  0.000000;
            aXihat( 0, 6 ) =  0.000000;
            aXihat( 1, 6 ) =  1.000000;
            aXihat( 0, 7 ) = -1.000000;
            aXihat( 1, 7 ) =  0.000000;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 8  >::eval_N(
                  Interpolation_Matrix  & aN,
            const Mat<real> & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( aN.n_cols() == 8,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( aN.n_rows() == 1,
                    "eval_shape: aN not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // often used constants
            auto  xi2 = std::pow(  xi, 2 );
            auto eta2 = std::pow( eta, 2 );

            // populate output matrix
            aN( 0 ) = - ( eta - 1.0 ) * ( xi - 1.0 ) * ( eta + xi + 1.0 ) * 0.25;
            aN( 1 ) =  ( eta - 1.0 ) * ( xi + 1.0 ) * ( eta - xi + 1.0 ) * 0.25;
            aN( 2 ) =  ( eta + 1.0 ) * ( xi + 1.0 ) * ( eta + xi - 1.0 ) * 0.25;
            aN( 3 ) =  ( eta + 1.0 ) * ( xi - 1.0 ) * (  - eta + xi + 1.0 ) * 0.25;
            aN( 4 ) =   ( xi2 - 1.0 ) * ( eta - 1.0 ) * 0.5;
            aN( 5 ) = - ( eta2 - 1.0 ) * ( xi + 1.0 ) * 0.5;
            aN( 6 ) = - ( xi2 - 1.0 ) * ( eta + 1.0 ) * 0.5;
            aN( 7 ) =   ( eta2 - 1.0 ) * ( xi - 1.0 ) * 0.5;


        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 8  >::eval_dNdXi(
                       Interpolation_Matrix & adNdXi,
                const Mat<real> & aXi) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( adNdXi.n_cols() == 8,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( adNdXi.n_rows() == 2,
                    "eval_shape: aN not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // often used constants
            auto  xi2 = std::pow(  xi, 2 );
            auto eta2 = std::pow( eta, 2 );

            // populate output matrix
            adNdXi( 0, 0 ) = -( eta+xi*2.0 )*( eta-1.0 )*0.25;
            adNdXi( 1, 0 ) = -( eta*2.0+xi )*( xi-1.0 )*0.25;

            adNdXi( 0, 1 ) = ( eta-xi*2.0 )*( eta-1.0 )*0.25;
            adNdXi( 1, 1 ) = ( xi+1.0 )*( eta*2.0-xi )*0.25;

            adNdXi( 0, 2 ) = ( eta+xi*2.0 )*( eta+1.0 )*0.25;
            adNdXi( 1, 2 ) = ( eta*2.0+xi )*( xi+1.0 )*0.25;

            adNdXi( 0, 3 ) = -( eta-xi*2.0 )*( eta+1.0 )*0.25;
            adNdXi( 1, 3 ) = -( xi-1.0 )*( eta*2.0-xi )*0.25;

            adNdXi( 0, 4 ) = xi*( eta-1.0 );
            adNdXi( 1, 4 ) = 0.5*xi2-0.5;

            adNdXi( 0, 5 ) = -0.5*eta2+0.5;
            adNdXi( 1, 5 ) = -eta*( xi+1.0 );

            adNdXi( 0, 6 ) = -xi*( eta+1.0 );
            adNdXi( 1, 6 ) = -0.5*xi2+0.5;

            adNdXi( 0, 7 ) = 0.5*eta2-0.5;
            adNdXi( 1, 7 ) = eta*( xi-1.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 8  >::eval_d2NdXi2(
                       Interpolation_Matrix & ad2NdXi2,
                const Mat<real> & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( ad2NdXi2.n_cols() == 8,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( ad2NdXi2.n_rows() == 3,
                    "eval_shape: aN not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // populate output matrix
            ad2NdXi2( 0, 0 ) = 0.5 - eta * 0.5;
            ad2NdXi2( 1, 0 ) = 0.5 - xi * 0.5;
            ad2NdXi2( 2, 0 ) = 0.25 - 0.5 * ( xi + eta );

            ad2NdXi2( 0, 1 ) = 0.5 - eta * 0.5;
            ad2NdXi2( 1, 1 ) = xi * 0.5 + 0.5;
            ad2NdXi2( 2, 1 ) = 0.5 * ( eta - xi ) - 0.25;

            ad2NdXi2( 0, 2 ) = eta * 0.5 + 0.5;
            ad2NdXi2( 1, 2 ) = xi * 0.5 + 0.5;
            ad2NdXi2( 2, 2 ) = 0.5 * ( xi + eta ) + 0.25;

            ad2NdXi2( 0, 3 ) = eta * 0.5 + 0.5;
            ad2NdXi2( 1, 3 ) = 0.5 - xi * 0.5;
            ad2NdXi2( 2, 3 ) = 0.5 * ( xi - eta ) - 0.25;

            ad2NdXi2( 0, 4 ) = eta - 1.0;
            ad2NdXi2( 1, 4 ) = 0;
            ad2NdXi2( 2, 4 ) = xi;

            ad2NdXi2( 0, 5 ) = 0.0;
            ad2NdXi2( 1, 5 ) =  - xi - 1.0;
            ad2NdXi2( 2, 5 ) =  - eta;

            ad2NdXi2( 0, 6 ) =  - eta - 1.0;
            ad2NdXi2( 1, 6 ) = 0.0;
            ad2NdXi2( 2, 6 ) =  - xi;

            ad2NdXi2( 0, 7 ) = 0.0;
            ad2NdXi2( 1, 7 ) = xi - 1.0;
            ad2NdXi2( 2, 7 ) = eta;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD8_HPP_ */
