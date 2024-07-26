/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Hex8.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_

#include "assert.h"
#include "moris_typedefs.hpp"                   //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src

// checked against femdoc:
// - N
// - dNdxi
// - d2Ndxi2
// - d3Ndxi3 by FD

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::get_number_of_param_dimensions() const
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::get_interpolation_order() const
        {
            return Interpolation_Order::LINEAR;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::get_param_coords( Matrix< DDRMat >& aXiHat ) const
        {
            aXiHat.set_size( 3, 8 );

            aXiHat( 0, 0 ) = -1.0;
            aXiHat( 1, 0 ) = -1.0;
            aXiHat( 2, 0 ) = -1.0;
            aXiHat( 0, 1 ) = 1.0;
            aXiHat( 1, 1 ) = -1.0;
            aXiHat( 2, 1 ) = -1.0;
            aXiHat( 0, 2 ) = 1.0;
            aXiHat( 1, 2 ) = 1.0;
            aXiHat( 2, 2 ) = -1.0;
            aXiHat( 0, 3 ) = -1.0;
            aXiHat( 1, 3 ) = 1.0;
            aXiHat( 2, 3 ) = -1.0;
            aXiHat( 0, 4 ) = -1.0;
            aXiHat( 1, 4 ) = -1.0;
            aXiHat( 2, 4 ) = 1.0;
            aXiHat( 0, 5 ) = 1.0;
            aXiHat( 1, 5 ) = -1.0;
            aXiHat( 2, 5 ) = 1.0;
            aXiHat( 0, 6 ) = 1.0;
            aXiHat( 1, 6 ) = 1.0;
            aXiHat( 2, 6 ) = 1.0;
            aXiHat( 0, 7 ) = -1.0;
            aXiHat( 1, 7 ) = 1.0;
            aXiHat( 2, 7 ) = 1.0;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_N(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            const real xi   = aXi( 0 );
            const real eta  = aXi( 1 );
            const real zeta = aXi( 2 );

            // populate output matrix
            aNXi.set_size( 1, 8 );

            // get raw pointer to memory
            real* tData = aNXi.data();

            // evaluate shape functions
            tData[ 0 ] = -( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            tData[ 1 ] = ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            tData[ 2 ] = -( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            tData[ 3 ] = ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            tData[ 4 ] = ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
            tData[ 5 ] = -( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            tData[ 6 ] = ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            tData[ 7 ] = -( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_dNdXi(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            const real xi   = aXi( 0 );
            const real eta  = aXi( 1 );
            const real zeta = aXi( 2 );

            // populate output matrix
            adNdXi.set_size( 3, 8 );

            // get raw pointer to memory
            real* tData = adNdXi.data();

            // evaluate derivative of shape functions
            tData[ 0 + 0 * 3 ] = -( eta - 1 ) * ( zeta - 1 ) * 0.125;
            tData[ 1 + 0 * 3 ] = -( xi - 1 ) * ( zeta - 1 ) * 0.125;
            tData[ 2 + 0 * 3 ] = -( eta - 1 ) * ( xi - 1 ) * 0.125;

            tData[ 0 + 1 * 3 ] = ( eta - 1 ) * ( zeta - 1 ) * 0.125;
            tData[ 1 + 1 * 3 ] = ( xi + 1 ) * ( zeta - 1 ) * 0.125;
            tData[ 2 + 1 * 3 ] = ( eta - 1 ) * ( xi + 1 ) * 0.125;

            tData[ 0 + 2 * 3 ] = -( eta + 1 ) * ( zeta - 1 ) * 0.125;
            tData[ 1 + 2 * 3 ] = -( xi + 1 ) * ( zeta - 1 ) * 0.125;
            tData[ 2 + 2 * 3 ] = -( eta + 1 ) * ( xi + 1 ) * 0.125;

            tData[ 0 + 3 * 3 ] = ( eta + 1 ) * ( zeta - 1 ) * 0.125;
            tData[ 1 + 3 * 3 ] = ( xi - 1 ) * ( zeta - 1 ) * 0.125;
            tData[ 2 + 3 * 3 ] = ( eta + 1 ) * ( xi - 1 ) * 0.125;

            tData[ 0 + 4 * 3 ] = ( eta - 1 ) * ( zeta + 1 ) * 0.125;
            tData[ 1 + 4 * 3 ] = ( xi - 1 ) * ( zeta + 1 ) * 0.125;
            tData[ 2 + 4 * 3 ] = ( eta - 1 ) * ( xi - 1 ) * 0.125;

            tData[ 0 + 5 * 3 ] = -( eta - 1 ) * ( zeta + 1 ) * 0.125;
            tData[ 1 + 5 * 3 ] = -( xi + 1 ) * ( zeta + 1 ) * 0.125;
            tData[ 2 + 5 * 3 ] = -( eta - 1 ) * ( xi + 1 ) * 0.125;

            tData[ 0 + 6 * 3 ] = ( eta + 1 ) * ( zeta + 1 ) * 0.125;
            tData[ 1 + 6 * 3 ] = ( xi + 1 ) * ( zeta + 1 ) * 0.125;
            tData[ 2 + 6 * 3 ] = ( eta + 1 ) * ( xi + 1 ) * 0.125;

            tData[ 0 + 7 * 3 ] = -( eta + 1 ) * ( zeta + 1 ) * 0.125;
            tData[ 1 + 7 * 3 ] = -( xi - 1 ) * ( zeta + 1 ) * 0.125;
            tData[ 2 + 7 * 3 ] = -( eta + 1 ) * ( xi - 1 ) * 0.125;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_d2NdXi2(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            const real xi   = aXi( 0 );
            const real eta  = aXi( 1 );
            const real zeta = aXi( 2 );

            // populate 6 x 8 matrix that is ordered as follows:
            // 1. row: d2N over dXi_1 dXi_1
            // 2. row: d2N over dXi_2 dXi_2
            // 3. row: d2N over dXi_3 dXi_3
            // 4. row: d2N over dXi_2 dXi_3
            // 5. row: d2N over dXi_1 dXi_3
            // 6. row: d2N over dXi_1 dXi_2

            ad2NdXi2.set_size( 6, 8, 0.0 );
            ad2NdXi2( 3, 0 ) = 0.125 * ( -xi + 1.0 );
            ad2NdXi2( 4, 0 ) = 0.125 * ( -eta + 1.0 );
            ad2NdXi2( 5, 0 ) = 0.125 * ( -zeta + 1.0 );

            ad2NdXi2( 3, 1 ) = 0.125 * ( xi + 1.0 );
            ad2NdXi2( 4, 1 ) = 0.125 * ( eta - 1.0 );
            ad2NdXi2( 5, 1 ) = 0.125 * ( zeta - 1.0 );

            ad2NdXi2( 3, 2 ) = 0.125 * ( -xi - 1.0 );
            ad2NdXi2( 4, 2 ) = 0.125 * ( -eta - 1.0 );
            ad2NdXi2( 5, 2 ) = 0.125 * ( -zeta + 1.0 );

            ad2NdXi2( 3, 3 ) = 0.125 * ( xi - 1.0 );
            ad2NdXi2( 4, 3 ) = 0.125 * ( eta + 1.0 );
            ad2NdXi2( 5, 3 ) = 0.125 * ( zeta - 1.0 );

            ad2NdXi2( 3, 4 ) = 0.125 * ( xi - 1.0 );
            ad2NdXi2( 4, 4 ) = 0.125 * ( eta - 1.0 );
            ad2NdXi2( 5, 4 ) = 0.125 * ( zeta + 1.0 );

            ad2NdXi2( 3, 5 ) = 0.125 * ( -xi - 1.0 );
            ad2NdXi2( 4, 5 ) = 0.125 * ( -eta + 1.0 );
            ad2NdXi2( 5, 5 ) = 0.125 * ( -zeta - 1.0 );

            ad2NdXi2( 3, 6 ) = 0.125 * ( xi + 1.0 );
            ad2NdXi2( 4, 6 ) = 0.125 * ( eta + 1.0 );
            ad2NdXi2( 5, 6 ) = 0.125 * ( zeta + 1.0 );

            ad2NdXi2( 3, 7 ) = 0.125 * ( -xi + 1.0 );
            ad2NdXi2( 4, 7 ) = 0.125 * ( -eta - 1.0 );
            ad2NdXi2( 5, 7 ) = 0.125 * ( -zeta - 1.0 );
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_d3NdXi3(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            // populate output matrix
            ad3NdXi3.set_size( 10, 8, 0.0 );

            ad3NdXi3( 9, 0 ) = -0.125;
            ad3NdXi3( 9, 1 ) = 0.125;
            ad3NdXi3( 9, 2 ) = -0.125;
            ad3NdXi3( 9, 3 ) = 0.125;
            ad3NdXi3( 9, 4 ) = 0.125;
            ad3NdXi3( 9, 5 ) = -0.125;
            ad3NdXi3( 9, 6 ) = 0.125;
            ad3NdXi3( 9, 7 ) = -0.125;
        }

        //------------------------------------------------------------------------------

    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_ */
