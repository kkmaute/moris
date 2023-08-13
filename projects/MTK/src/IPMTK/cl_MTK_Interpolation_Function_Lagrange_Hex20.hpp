/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Hex20.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_HEX20_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_HEX20_HPP_

#include "assert.h"
#include "typedefs.hpp"                         //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src

// checked against femdoc:
// - N
// - dNdxi
// - d2Ndxi2

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 20 >::get_number_of_param_dimensions() const
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 20 >::get_interpolation_order() const
        {
            return Interpolation_Order::SERENDIPITY;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 20 >::get_param_coords( Matrix< DDRMat >& aXiHat ) const
        {
            aXiHat.set_size( 3, 20, 0.0 );
            aXiHat( 0, 0 )  = -1.000000;
            aXiHat( 1, 0 )  = -1.000000;
            aXiHat( 2, 0 )  = -1.000000;
            aXiHat( 0, 1 )  = 1.000000;
            aXiHat( 1, 1 )  = -1.000000;
            aXiHat( 2, 1 )  = -1.000000;
            aXiHat( 0, 2 )  = 1.000000;
            aXiHat( 1, 2 )  = 1.000000;
            aXiHat( 2, 2 )  = -1.000000;
            aXiHat( 0, 3 )  = -1.000000;
            aXiHat( 1, 3 )  = 1.000000;
            aXiHat( 2, 3 )  = -1.000000;
            aXiHat( 0, 4 )  = -1.000000;
            aXiHat( 1, 4 )  = -1.000000;
            aXiHat( 2, 4 )  = 1.000000;
            aXiHat( 0, 5 )  = 1.000000;
            aXiHat( 1, 5 )  = -1.000000;
            aXiHat( 2, 5 )  = 1.000000;
            aXiHat( 0, 6 )  = 1.000000;
            aXiHat( 1, 6 )  = 1.000000;
            aXiHat( 2, 6 )  = 1.000000;
            aXiHat( 0, 7 )  = -1.000000;
            aXiHat( 1, 7 )  = 1.000000;
            aXiHat( 2, 7 )  = 1.000000;
            aXiHat( 0, 8 )  = 0.000000;
            aXiHat( 1, 8 )  = -1.000000;
            aXiHat( 2, 8 )  = -1.000000;
            aXiHat( 0, 9 )  = 1.000000;
            aXiHat( 1, 9 )  = 0.000000;
            aXiHat( 2, 9 )  = -1.000000;
            aXiHat( 0, 10 ) = 0.000000;
            aXiHat( 1, 10 ) = 1.000000;
            aXiHat( 2, 10 ) = -1.000000;
            aXiHat( 0, 11 ) = -1.000000;
            aXiHat( 1, 11 ) = 0.000000;
            aXiHat( 2, 11 ) = -1.000000;
            aXiHat( 0, 12 ) = -1.000000;
            aXiHat( 1, 12 ) = -1.000000;
            aXiHat( 2, 12 ) = 0.000000;
            aXiHat( 0, 13 ) = 1.000000;
            aXiHat( 1, 13 ) = -1.000000;
            aXiHat( 2, 13 ) = 0.000000;
            aXiHat( 0, 14 ) = 1.000000;
            aXiHat( 1, 14 ) = 1.000000;
            aXiHat( 2, 14 ) = 0.000000;
            aXiHat( 0, 15 ) = -1.000000;
            aXiHat( 1, 15 ) = 1.000000;
            aXiHat( 2, 15 ) = 0.000000;
            aXiHat( 0, 16 ) = 0.000000;
            aXiHat( 1, 16 ) = -1.000000;
            aXiHat( 2, 16 ) = 1.000000;
            aXiHat( 0, 17 ) = 1.000000;
            aXiHat( 1, 17 ) = 0.000000;
            aXiHat( 2, 17 ) = 1.000000;
            aXiHat( 0, 18 ) = 0.000000;
            aXiHat( 1, 18 ) = 1.000000;
            aXiHat( 2, 18 ) = 1.000000;
            aXiHat( 0, 19 ) = -1.000000;
            aXiHat( 1, 19 ) = 0.000000;
            aXiHat( 2, 19 ) = 1.000000;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 20 >::eval_N(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX20 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real xi   = aXi( 0 );
            real eta  = aXi( 1 );
            real zeta = aXi( 2 );

            // often used constants
            real xi2   = std::pow( xi, 2 );
            real eta2  = std::pow( eta, 2 );
            real zeta2 = std::pow( zeta, 2 );

            // populate output matrix
            aNXi.set_size( 1, 20 );
            aNXi( 0 )  = ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * ( eta + xi + zeta + 2.0 ) * 0.125;
            aNXi( 1 )  = -( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * ( eta - xi + zeta + 2.0 ) * 0.125;
            aNXi( 2 )  = -( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * ( eta + xi - zeta - 2.0 ) * 0.125;
            aNXi( 3 )  = -( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * ( -eta + xi + zeta + 2.0 ) * 0.125;
            aNXi( 4 )  = -( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * ( eta + xi - zeta + 2.0 ) * 0.125;
            aNXi( 5 )  = ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * ( eta - xi - zeta + 2.0 ) * 0.125;
            aNXi( 6 )  = ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * ( eta + xi + zeta - 2.0 ) * 0.125;
            aNXi( 7 )  = -( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * ( eta - xi + zeta - 2.0 ) * 0.125;
            aNXi( 8 )  = -( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 ) * 0.25;
            aNXi( 9 )  = ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.25;
            aNXi( 10 ) = ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 ) * 0.25;
            aNXi( 11 ) = -( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.25;
            aNXi( 12 ) = -( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 ) * 0.25;
            aNXi( 13 ) = ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 ) * 0.25;
            aNXi( 14 ) = -( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 ) * 0.25;
            aNXi( 15 ) = ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 ) * 0.25;
            aNXi( 16 ) = ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 ) * 0.25;
            aNXi( 17 ) = -( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.25;
            aNXi( 18 ) = -( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 ) * 0.25;
            aNXi( 19 ) = ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.25;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 20 >::eval_dNdXi( const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&                                                                                              adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX20 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real xi   = aXi( 0 );
            real eta  = aXi( 1 );
            real zeta = aXi( 2 );

            // often used constants
            real xi2   = std::pow( xi, 2 );
            real eta2  = std::pow( eta, 2 );
            real zeta2 = std::pow( zeta, 2 );

            // populate adNdXi
            adNdXi.set_size( 3, 20 );
            adNdXi( 0, 0 ) = ( ( eta - 1.0 ) * ( zeta - 1.0 ) * ( eta + 2.0 * xi + zeta + 1.0 ) ) * 0.125;
            adNdXi( 1, 0 ) = ( ( xi - 1.0 ) * ( zeta - 1.0 ) * ( 2.0 * eta + xi + zeta + 1.0 ) ) * 0.125;
            adNdXi( 2, 0 ) = ( ( eta - 1.0 ) * ( xi - 1.0 ) * ( eta + xi + 2.0 * zeta + 1.0 ) ) * 0.125;

            adNdXi( 0, 1 ) = -( ( eta - 1.0 ) * ( zeta - 1.0 ) * ( eta - 2.0 * xi + zeta + 1.0 ) ) * 0.125;
            adNdXi( 1, 1 ) = -( ( xi + 1.0 ) * ( zeta - 1.0 ) * ( 2.0 * eta - xi + zeta + 1.0 ) ) * 0.125;
            adNdXi( 2, 1 ) = -( ( eta - 1.0 ) * ( xi + 1.0 ) * ( eta - xi + 2.0 * zeta + 1.0 ) ) * 0.125;

            adNdXi( 0, 2 ) = -( ( eta + 1.0 ) * ( zeta - 1.0 ) * ( eta + 2.0 * xi - zeta - 1.0 ) ) * 0.125;
            adNdXi( 1, 2 ) = -( ( xi + 1.0 ) * ( zeta - 1.0 ) * ( 2.0 * eta + xi - zeta - 1.0 ) ) * 0.125;
            adNdXi( 2, 2 ) = -( ( eta + 1.0 ) * ( xi + 1.0 ) * ( eta + xi - 2.0 * zeta - 1.0 ) ) * 0.125;

            adNdXi( 0, 3 ) = -( ( eta + 1.0 ) * ( zeta - 1.0 ) * ( 2.0 * xi - eta + zeta + 1.0 ) ) * 0.125;
            adNdXi( 1, 3 ) = -( ( xi - 1.0 ) * ( zeta - 1.0 ) * ( xi - 2.0 * eta + zeta + 1.0 ) ) * 0.125;
            adNdXi( 2, 3 ) = -( ( eta + 1.0 ) * ( xi - 1.0 ) * ( xi - eta + 2.0 * zeta + 1.0 ) ) * 0.125;

            adNdXi( 0, 4 ) = -( ( eta - 1.0 ) * ( zeta + 1.0 ) * ( eta + 2.0 * xi - zeta + 1.0 ) ) * 0.125;
            adNdXi( 1, 4 ) = -( ( xi - 1.0 ) * ( zeta + 1.0 ) * ( 2.0 * eta + xi - zeta + 1.0 ) ) * 0.125;
            adNdXi( 2, 4 ) = -( ( eta - 1.0 ) * ( xi - 1.0 ) * ( eta + xi - 2.0 * zeta + 1.0 ) ) * 0.125;

            adNdXi( 0, 5 ) = ( ( eta - 1.0 ) * ( zeta + 1.0 ) * ( eta - 2.0 * xi - zeta + 1.0 ) ) * 0.125;
            adNdXi( 1, 5 ) = ( ( xi + 1.0 ) * ( zeta + 1.0 ) * ( 2.0 * eta - xi - zeta + 1.0 ) ) * 0.125;
            adNdXi( 2, 5 ) = ( ( eta - 1.0 ) * ( xi + 1.0 ) * ( eta - xi - 2.0 * zeta + 1.0 ) ) * 0.125;

            adNdXi( 0, 6 ) = ( ( eta + 1.0 ) * ( zeta + 1.0 ) * ( eta + 2.0 * xi + zeta - 1.0 ) ) * 0.125;
            adNdXi( 1, 6 ) = ( ( xi + 1.0 ) * ( zeta + 1.0 ) * ( 2.0 * eta + xi + zeta - 1.0 ) ) * 0.125;
            adNdXi( 2, 6 ) = ( ( eta + 1.0 ) * ( xi + 1.0 ) * ( eta + xi + 2.0 * zeta - 1.0 ) ) * 0.125;

            adNdXi( 0, 7 ) = -( ( eta + 1.0 ) * ( zeta + 1.0 ) * ( eta - 2.0 * xi + zeta - 1.0 ) ) * 0.125;
            adNdXi( 1, 7 ) = -( ( xi - 1.0 ) * ( zeta + 1.0 ) * ( 2.0 * eta - xi + zeta - 1.0 ) ) * 0.125;
            adNdXi( 2, 7 ) = -( ( eta + 1.0 ) * ( xi - 1.0 ) * ( eta - xi + 2.0 * zeta - 1.0 ) ) * 0.125;

            adNdXi( 0, 8 ) = -( xi * ( eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            adNdXi( 1, 8 ) = -( ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            adNdXi( 2, 8 ) = -( ( xi2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;

            adNdXi( 0, 9 ) = ( ( eta2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            adNdXi( 1, 9 ) = ( eta * ( xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            adNdXi( 2, 9 ) = ( ( eta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;

            adNdXi( 0, 10 ) = ( xi * ( eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            adNdXi( 1, 10 ) = ( ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            adNdXi( 2, 10 ) = ( ( xi2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;

            adNdXi( 0, 11 ) = -( ( eta2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            adNdXi( 1, 11 ) = -( eta * ( xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            adNdXi( 2, 11 ) = -( ( eta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;

            adNdXi( 0, 12 ) = -( ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            adNdXi( 1, 12 ) = -( ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
            adNdXi( 2, 12 ) = -( zeta * ( eta - 1.0 ) * ( xi - 1.0 ) ) * 0.5;

            adNdXi( 0, 13 ) = ( ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            adNdXi( 1, 13 ) = ( ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
            adNdXi( 2, 13 ) = ( zeta * ( eta - 1.0 ) * ( xi + 1.0 ) ) * 0.5;

            adNdXi( 0, 14 ) = -( ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            adNdXi( 1, 14 ) = -( ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
            adNdXi( 2, 14 ) = -( zeta * ( eta + 1.0 ) * ( xi + 1.0 ) ) * 0.5;

            adNdXi( 0, 15 ) = ( ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            adNdXi( 1, 15 ) = ( ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
            adNdXi( 2, 15 ) = ( zeta * ( eta + 1.0 ) * ( xi - 1.0 ) ) * 0.5;

            adNdXi( 0, 16 ) = ( xi * ( eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            adNdXi( 1, 16 ) = ( ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            adNdXi( 2, 16 ) = ( ( xi2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;

            adNdXi( 0, 17 ) = -( ( eta2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            adNdXi( 1, 17 ) = -( eta * ( xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            adNdXi( 2, 17 ) = -( ( eta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;

            adNdXi( 0, 18 ) = -( xi * ( eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            adNdXi( 1, 18 ) = -( ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            adNdXi( 2, 18 ) = -( ( xi2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;

            adNdXi( 0, 19 ) = ( ( eta2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            adNdXi( 1, 19 ) = ( eta * ( xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            adNdXi( 2, 19 ) = ( ( eta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
        }
        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 20 >::eval_d2NdXi2( const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&                                                                                                ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX20 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real xi   = aXi( 0 );
            real eta  = aXi( 1 );
            real zeta = aXi( 2 );

            // often used constants
            real xi2   = std::pow( xi, 2 );
            real eta2  = std::pow( eta, 2 );
            real zeta2 = std::pow( zeta, 2 );

            ad2NdXi2.set_size( 6, 20 );
            ad2NdXi2( 0, 0 ) = ( ( eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 1, 0 ) = ( ( xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 2, 0 ) = ( ( eta - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
            ad2NdXi2( 3, 0 ) = ( ( xi - 1.0 ) * ( 2.0 * ( eta + zeta ) + xi ) ) * 0.125;
            ad2NdXi2( 4, 0 ) = ( ( eta - 1.0 ) * ( 2.0 * ( xi + zeta ) + eta ) ) * 0.125;
            ad2NdXi2( 5, 0 ) = ( ( zeta - 1.0 ) * ( 2.0 * ( xi + eta ) + zeta ) ) * 0.125;

            ad2NdXi2( 0, 1 ) = ( ( eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 1, 1 ) = -( ( xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 2, 1 ) = -( ( eta - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
            ad2NdXi2( 3, 1 ) = -( ( xi + 1.0 ) * ( 2.0 * ( eta + zeta ) - xi ) ) * 0.125;
            ad2NdXi2( 4, 1 ) = -( ( eta - 1.0 ) * ( eta - 2.0 * ( xi - zeta ) ) ) * 0.125;
            ad2NdXi2( 5, 1 ) = -( ( zeta - 1.0 ) * ( 2.0 * ( eta - xi ) + zeta ) ) * 0.125;

            ad2NdXi2( 0, 2 ) = -( ( eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 1, 2 ) = -( ( xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 2, 2 ) = ( ( eta + 1.0 ) * ( xi + 1.0 ) ) * 0.25;
            ad2NdXi2( 3, 2 ) = -( ( xi + 1.0 ) * ( 2.0 * eta + xi - 2.0 * zeta ) ) * 0.125;
            ad2NdXi2( 4, 2 ) = -( ( eta + 1.0 ) * ( eta + 2.0 * xi - 2.0 * zeta ) ) * 0.125;
            ad2NdXi2( 5, 2 ) = -( ( zeta - 1.0 ) * ( 2.0 * eta + 2.0 * xi - zeta ) ) * 0.125;

            ad2NdXi2( 0, 3 ) = -( ( eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 1, 3 ) = ( ( xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 2, 3 ) = -( ( eta + 1.0 ) * ( xi - 1.0 ) ) * 0.25;
            ad2NdXi2( 3, 3 ) = -( ( xi - 1.0 ) * ( xi - 2.0 * ( eta - zeta ) ) ) * 0.125;
            ad2NdXi2( 4, 3 ) = -( ( eta + 1.0 ) * ( 2.0 * ( xi + zeta ) - eta ) ) * 0.125;
            ad2NdXi2( 5, 3 ) = -( ( zeta - 1.0 ) * ( zeta - 2.0 * ( eta - xi ) ) ) * 0.125;

            ad2NdXi2( 0, 4 ) = -( ( eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 1, 4 ) = -( ( xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 2, 4 ) = ( ( eta - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
            ad2NdXi2( 3, 4 ) = -( ( xi - 1.0 ) * ( 2.0 * eta + xi - 2.0 * zeta ) ) * 0.125;
            ad2NdXi2( 4, 4 ) = -( ( eta - 1.0 ) * ( eta + 2.0 * xi - 2.0 * zeta ) ) * 0.125;
            ad2NdXi2( 5, 4 ) = -( ( zeta + 1.0 ) * ( 2.0 * eta + 2.0 * xi - zeta ) ) * 0.125;

            ad2NdXi2( 0, 5 ) = -( ( eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 1, 5 ) = ( ( xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 2, 5 ) = -( ( eta - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
            ad2NdXi2( 3, 5 ) = -( ( xi + 1.0 ) * ( xi - 2.0 * ( eta - zeta ) ) ) * 0.125;
            ad2NdXi2( 4, 5 ) = -( ( eta - 1.0 ) * ( 2.0 * ( xi + zeta ) - eta ) ) * 0.125;
            ad2NdXi2( 5, 5 ) = -( ( zeta + 1.0 ) * ( zeta - 2.0 * ( eta - xi ) ) ) * 0.125;

            ad2NdXi2( 0, 6 ) = ( ( eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 1, 6 ) = ( ( xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 2, 6 ) = ( ( eta + 1.0 ) * ( xi + 1.0 ) ) * 0.25;
            ad2NdXi2( 3, 6 ) = ( ( xi + 1.0 ) * ( 2.0 * ( eta + zeta ) + xi ) ) * 0.125;
            ad2NdXi2( 4, 6 ) = ( ( eta + 1.0 ) * ( 2.0 * ( xi + zeta ) + eta ) ) * 0.125;
            ad2NdXi2( 5, 6 ) = ( ( zeta + 1.0 ) * ( 2.0 * ( xi + eta ) + zeta ) ) * 0.125;

            ad2NdXi2( 0, 7 ) = ( ( eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 1, 7 ) = -( ( xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 2, 7 ) = -( ( eta + 1.0 ) * ( xi - 1.0 ) ) * 0.25;
            ad2NdXi2( 3, 7 ) = -( ( xi - 1.0 ) * ( 2.0 * ( eta + zeta ) - xi ) ) * 0.125;
            ad2NdXi2( 4, 7 ) = -( ( eta + 1.0 ) * ( eta - 2.0 * ( xi - zeta ) ) ) * 0.125;
            ad2NdXi2( 5, 7 ) = -( ( zeta + 1.0 ) * ( 2.0 * ( eta - xi ) + zeta ) ) * 0.125;

            ad2NdXi2( 0, 8 ) = -( ( eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 8 ) = 0.0;
            ad2NdXi2( 2, 8 ) = 0.0;
            ad2NdXi2( 3, 8 ) = 0.25 - xi2 * 0.25;
            ad2NdXi2( 4, 8 ) = -( xi * ( eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 8 ) = -( xi * ( zeta - 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 9 ) = 0.0;
            ad2NdXi2( 1, 9 ) = ( ( xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 9 ) = 0.0;
            ad2NdXi2( 3, 9 ) = ( eta * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 9 ) = eta2 * 0.25 - 0.25;
            ad2NdXi2( 5, 9 ) = ( eta * ( zeta - 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 10 ) = ( ( eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 10 ) = 0.0;
            ad2NdXi2( 2, 10 ) = 0.0;
            ad2NdXi2( 3, 10 ) = xi2 * 0.25 - 0.25;
            ad2NdXi2( 4, 10 ) = ( xi * ( eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 10 ) = ( xi * ( zeta - 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 11 ) = 0.0;
            ad2NdXi2( 1, 11 ) = -( ( xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 11 ) = 0.0;
            ad2NdXi2( 3, 11 ) = -( eta * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 11 ) = 0.25 - eta2 * 0.25;
            ad2NdXi2( 5, 11 ) = -( eta * ( zeta - 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 12 ) = 0.0;
            ad2NdXi2( 1, 12 ) = 0.0;
            ad2NdXi2( 2, 12 ) = -( ( eta - 1.0 ) * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 12 ) = -( zeta * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 12 ) = -( zeta * ( eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 12 ) = 0.25 - zeta2 * 0.25;

            ad2NdXi2( 0, 13 ) = 0.0;
            ad2NdXi2( 1, 13 ) = 0.0;
            ad2NdXi2( 2, 13 ) = ( ( eta - 1.0 ) * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 13 ) = ( zeta * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 13 ) = ( zeta * ( eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 13 ) = zeta2 * 0.25 - 0.25;

            ad2NdXi2( 0, 14 ) = 0.0;
            ad2NdXi2( 1, 14 ) = 0.0;
            ad2NdXi2( 2, 14 ) = -( ( eta + 1.0 ) * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 14 ) = -( zeta * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 14 ) = -( zeta * ( eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 14 ) = 0.25 - zeta2 * 0.25;

            ad2NdXi2( 0, 15 ) = 0.0;
            ad2NdXi2( 1, 15 ) = 0.0;
            ad2NdXi2( 2, 15 ) = ( ( eta + 1.0 ) * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 15 ) = ( zeta * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 15 ) = ( zeta * ( eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 15 ) = zeta2 * 0.25 - 0.25;

            ad2NdXi2( 0, 16 ) = ( ( eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 16 ) = 0.0;
            ad2NdXi2( 2, 16 ) = 0.0;
            ad2NdXi2( 3, 16 ) = xi2 * 0.25 - 0.25;
            ad2NdXi2( 4, 16 ) = ( xi * ( eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 16 ) = ( xi * ( zeta + 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 17 ) = 0.0;
            ad2NdXi2( 1, 17 ) = -( ( xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 17 ) = 0.0;
            ad2NdXi2( 3, 17 ) = -( eta * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 17 ) = 0.25 - eta2 * 0.25;
            ad2NdXi2( 5, 17 ) = -( eta * ( zeta + 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 18 ) = -( ( eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 18 ) = 0.0;
            ad2NdXi2( 2, 18 ) = 0.0;
            ad2NdXi2( 3, 18 ) = 0.25 - xi2 * 0.25;
            ad2NdXi2( 4, 18 ) = -( xi * ( eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 18 ) = -( xi * ( zeta + 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 19 ) = 0.0;
            ad2NdXi2( 1, 19 ) = ( ( xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 19 ) = 0.0;
            ad2NdXi2( 3, 19 ) = ( eta * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 19 ) = eta2 * 0.25 - 0.25;
            ad2NdXi2( 5, 19 ) = ( eta * ( zeta + 1.0 ) ) * 0.5;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 20 >::eval_d3NdXi3( const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&                                                                                                ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( false, "HEX20 - eval_d3NdXi3: 3rd order derivatives not implemented for this element" );

            ad3NdXi3.set_size( 10, 20, 0.0 );
        }

        //------------------------------------------------------------------------------

    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_HEX20_HPP_ */
