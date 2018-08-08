/*
 * cl_FEM_Interpolation_Function_Lagrange_Hex8.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_

#include "assert.h"

#include "cl_FEM_Interpolation_Matrix.hpp"
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
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 8  >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 8  >::get_param_coords(
                Mat<real> & aXihat ) const
        {
            aXihat.set_size( 3, 8 );

            aXihat( 0, 0 ) = -1.000000;
            aXihat( 1, 0 ) = -1.000000;
            aXihat( 2, 0 ) = -1.000000;
            aXihat( 0, 1 ) =  1.000000;
            aXihat( 1, 1 ) = -1.000000;
            aXihat( 2, 1 ) = -1.000000;
            aXihat( 0, 2 ) =  1.000000;
            aXihat( 1, 2 ) =  1.000000;
            aXihat( 2, 2 ) = -1.000000;
            aXihat( 0, 3 ) = -1.000000;
            aXihat( 1, 3 ) =  1.000000;
            aXihat( 2, 3 ) = -1.000000;
            aXihat( 0, 4 ) = -1.000000;
            aXihat( 1, 4 ) = -1.000000;
            aXihat( 2, 4 ) =  1.000000;
            aXihat( 0, 5 ) =  1.000000;
            aXihat( 1, 5 ) = -1.000000;
            aXihat( 2, 5 ) =  1.000000;
            aXihat( 0, 6 ) =  1.000000;
            aXihat( 1, 6 ) =  1.000000;
            aXihat( 2, 6 ) =  1.000000;
            aXihat( 0, 7 ) = -1.000000;
            aXihat( 1, 7 ) =  1.000000;
            aXihat( 2, 7 ) =  1.000000;
       }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 8  >::eval_N(
                  Interpolation_Matrix  & aN,
            const Mat<real> & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( aN.n_cols() == 8,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( aN.n_rows() == 1,
                    "eval_shape: aN not allocated or hat wrong size." );


            // unpack xi and eta from input vector
            auto    xi = aXi( 0 );
            auto   eta = aXi( 1 );
            auto  zeta = aXi( 2 );

            // populate output matrix
            aN( 0 ) =  - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 1 ) =    ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 2 ) =  - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 3 ) =    ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aN( 4 ) =    ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aN( 5 ) =  - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aN( 6 ) =    ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aN( 7 ) =  - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 8  >::eval_dNdXi(
                      Interpolation_Matrix  & adNdXi,
                const Mat<real> & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( adNdXi.n_cols() == 8,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( adNdXi.n_rows() == 3,
                    "eval_shape: aN not allocated or hat wrong size." );


            // unpack xi and eta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // populate output matrix
            adNdXi( 0, 0 ) = -(  eta - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 0 ) = -(   xi - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 0 ) = -(  eta - 1 ) * (   xi - 1 ) * 0.125;

            adNdXi( 0, 1 ) =  (  eta - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 1 ) =  (   xi + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 1 ) =  (  eta - 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 2 ) = -(  eta + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 2 ) = -(   xi + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 2 ) = -(  eta + 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 3 ) =  (  eta + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 3 ) =  (   xi - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 3 ) =  (  eta + 1 ) * (   xi - 1 ) * 0.125;

            adNdXi( 0, 4 ) =  (  eta - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 4 ) =  (   xi - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 4 ) =  (  eta - 1 ) * (   xi - 1 ) * 0.125;

            adNdXi( 0, 5 ) = -(  eta - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 5 ) = -(   xi + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 5 ) = -(  eta - 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 6 ) =  (  eta + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 6 ) =  (   xi + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 6 ) =  (  eta + 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 7 ) = -(  eta + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 7 ) = -(   xi - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 7 ) = -(  eta + 1 ) * (   xi - 1 ) * 0.125;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 8  >::eval_d2NdXi2(
                      Interpolation_Matrix  & ad2NdXi2,
                const Mat<real> & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( ad2NdXi2.n_cols() == 8,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( ad2NdXi2.n_rows() == 6,
                    "eval_shape: aN not allocated or hat wrong size." );


            // unpack xi and eta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // populate output matrix

            ad2NdXi2( 0, 0 ) = 0.0;
            ad2NdXi2( 1, 0 ) = 0.0;
            ad2NdXi2( 2, 0 ) = 0.0;
            ad2NdXi2( 3, 0 ) = 0.125 * ( - xi + 1.0 );
            ad2NdXi2( 4, 0 ) = 0.125 * ( - eta + 1.0 );
            ad2NdXi2( 5, 0 ) = 0.125 * ( - zeta + 1.0 );

            ad2NdXi2( 0, 1 ) = 0.0;
            ad2NdXi2( 1, 1 ) = 0.0;
            ad2NdXi2( 2, 1 ) = 0.0;
            ad2NdXi2( 3, 1 ) = 0.125 * ( xi + 1.0 );
            ad2NdXi2( 4, 1 ) = 0.125 * ( eta - 1.0 );
            ad2NdXi2( 5, 1 ) = 0.125 * ( zeta - 1.0 );

            ad2NdXi2( 0, 2 ) = 0.0;
            ad2NdXi2( 1, 2 ) = 0.0;
            ad2NdXi2( 2, 2 ) = 0.0;
            ad2NdXi2( 3, 2 ) = 0.125 * ( - xi - 1.0 );
            ad2NdXi2( 4, 2 ) = 0.125 * ( - eta - 1.0 );
            ad2NdXi2( 5, 2 ) = 0.125 * ( - zeta + 1.0 );

            ad2NdXi2( 0, 3 ) = 0.0;
            ad2NdXi2( 1, 3 ) = 0.0;
            ad2NdXi2( 2, 3 ) = 0.0;
            ad2NdXi2( 3, 3 ) = 0.125 * ( xi - 1.0 );
            ad2NdXi2( 4, 3 ) = 0.125 * ( eta + 1.0 );
            ad2NdXi2( 5, 3 ) = 0.125 * ( zeta - 1.0 );

            ad2NdXi2( 0, 4 ) = 0.0;
            ad2NdXi2( 1, 4 ) = 0.0;
            ad2NdXi2( 2, 4 ) = 0.0;
            ad2NdXi2( 3, 4 ) = 0.125 * ( xi - 1.0 );
            ad2NdXi2( 4, 4 ) = 0.125 * ( eta - 1.0 );
            ad2NdXi2( 5, 4 ) = 0.125 * ( zeta + 1.0 );

            ad2NdXi2( 0, 5 ) = 0.0;
            ad2NdXi2( 1, 5 ) = 0.0;
            ad2NdXi2( 2, 5 ) = 0.0;
            ad2NdXi2( 3, 5 ) = 0.125 * ( - xi - 1.0 );
            ad2NdXi2( 4, 5 ) = 0.125 * ( - eta + 1.0 );
            ad2NdXi2( 5, 5 ) = 0.125 * ( - zeta - 1.0 );

            ad2NdXi2( 0, 6 ) = 0.0;
            ad2NdXi2( 1, 6 ) = 0.0;
            ad2NdXi2( 2, 6 ) = 0.0;
            ad2NdXi2( 3, 6 ) = 0.125 * ( xi + 1.0 );
            ad2NdXi2( 4, 6 ) = 0.125 * ( eta + 1.0 );
            ad2NdXi2( 5, 6 ) = 0.125 * ( zeta + 1.0 );

            ad2NdXi2( 0, 7 ) = 0.0;
            ad2NdXi2( 1, 7 ) = 0.0;
            ad2NdXi2( 2, 7 ) = 0.0;
            ad2NdXi2( 3, 7 ) = 0.125 * ( - xi + 1.0 );
            ad2NdXi2( 4, 7 ) = 0.125 * ( - eta - 1.0 );
            ad2NdXi2( 5, 7 ) = 0.125 * ( - zeta - 1.0 );
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_ */
