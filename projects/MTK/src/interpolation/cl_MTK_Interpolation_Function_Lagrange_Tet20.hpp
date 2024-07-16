/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Tet20.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET20_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET20_HPP_

#include "assert.h"

#include "moris_typedefs.hpp"                         //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::get_number_of_param_dimensions() const
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::get_interpolation_order() const
        {
            return Interpolation_Order::CUBIC;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::get_param_coords( Matrix< DDRMat >& aXiHat ) const
        {
            real t13 = 1.0 / 3.0;
            real t23 = 2.0 / 3.0;

            aXiHat = {
                { 1.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t23, t13, t23, t13, 0.0, 0.0, 0.0, 0.0, t13, t13, 0.0, t13 },
                { 0.0, 1.0, 0.0, 0.0, t13, t23, t23, t13, 0.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t13, t13, t13, 0.0 },
                { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, t13, t23, t13, t23, t13, t23, 0.0, t13, t13, t13 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_N(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TET20 - eval_N: aXi not allocated or hat wrong size." );

            // unpack the tet coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = 1.0 - aXi( 0 ) - aXi( 1 ) - aXi( 2 );
            real zeta4 = aXi( 2 );

            // populate matrix with values
            aNXi.set_size( 1, 20 );
            aNXi( 0 ) = 0.5 * zeta1 * ( 3.0 * zeta1 - 1.0 ) * ( 3.0 * zeta1 - 2.0 );
            aNXi( 1 ) = 0.5 * zeta2 * ( 3.0 * zeta2 - 1.0 ) * ( 3.0 * zeta2 - 2.0 );
            aNXi( 2 ) = 0.5 * zeta3 * ( 3.0 * zeta3 - 1.0 ) * ( 3.0 * zeta3 - 2.0 );
            aNXi( 3 ) = 0.5 * zeta4 * ( 3.0 * zeta4 - 1.0 ) * ( 3.0 * zeta4 - 2.0 );

            aNXi( 4 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta1 - 1.0 );
            aNXi( 5 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta2 - 1.0 );

            aNXi( 6 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta2 - 1.0 );
            aNXi( 7 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta3 - 1.0 );

            aNXi( 8 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta1 - 1.0 );
            aNXi( 9 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta3 - 1.0 );

            aNXi( 10 ) = 4.5 * zeta1 * zeta4 * ( 3.0 * zeta1 - 1.0 );
            aNXi( 11 ) = 4.5 * zeta1 * zeta4 * ( 3.0 * zeta4 - 1.0 );

            aNXi( 12 ) = 4.5 * zeta2 * zeta4 * ( 3.0 * zeta2 - 1.0 );
            aNXi( 13 ) = 4.5 * zeta2 * zeta4 * ( 3.0 * zeta4 - 1.0 );

            aNXi( 14 ) = 4.5 * zeta3 * zeta4 * ( 3.0 * zeta3 - 1.0 );
            aNXi( 15 ) = 4.5 * zeta3 * zeta4 * ( 3.0 * zeta4 - 1.0 );

            aNXi( 16 ) = 27.0 * zeta1 * zeta2 * zeta3;
            aNXi( 17 ) = 27.0 * zeta1 * zeta2 * zeta4;
            aNXi( 18 ) = 27.0 * zeta2 * zeta3 * zeta4;
            aNXi( 19 ) = 27.0 * zeta1 * zeta3 * zeta4;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_dNdXi(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TET20 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack the tet coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta4 = aXi( 2 );

            // populate output matrix
            adNdXi.set_size( 3, 20 );

            adNdXi( 0, 0 ) = ( 27.0 * zeta1 * zeta1 ) / 2.0 - 9.0 * zeta1 + 1.0;
            adNdXi( 0, 1 ) = 0.0;
            adNdXi( 0, 2 ) = 18.0 * zeta1 + 18.0 * zeta2 + 18.0 * zeta4 - 27.0 * zeta1 * zeta2 - 27.0 * zeta1 * zeta4 - 27.0 * zeta2 * zeta4    //
                           - ( 27.0 * zeta1 * zeta1 ) / 2.0 - ( 27.0 * zeta2 * zeta2 ) / 2.0 - ( 27.0 * zeta4 * zeta4 ) / 2.0 - 11.0 / 2.0;
            adNdXi( 0, 3 ) = 0;
            adNdXi( 0, 4 ) = ( 9.0 * zeta2 * ( 6.0 * zeta1 - 1.0 ) ) / 2.0;
            adNdXi( 0, 5 ) = ( 9.0 * zeta2 * ( 3.0 * zeta2 - 1.0 ) ) / 2.0;
            adNdXi( 0, 6 ) = -( 9.0 * zeta2 * ( 3.0 * zeta2 - 1.0 ) ) / 2.0;
            adNdXi( 0, 7 ) = ( 9.0 * zeta2 * ( 6.0 * zeta1 + 6.0 * zeta2 + 6.0 * zeta4 - 5.0 ) ) / 2.0;
            adNdXi( 0, 8 ) = 36.0 * zeta1 + ( 9.0 * zeta2 ) / 2.0 + ( 9.0 * zeta4 ) / 2.0 - 27.0 * zeta1 * zeta2 - 27.0 * zeta1 * zeta4    //
                           - ( 81.0 * zeta1 * zeta1 ) / 2.0 - 9.0 / 2.0;
            adNdXi( 0, 9 ) = ( 9.0 * ( zeta1 + zeta2 + zeta4 - 1.0 ) * ( 3.0 * zeta1 + 3.0 * zeta2 + 3.0 * zeta4 - 2.0 ) ) / 2.0    //
                           + ( 27.0 * zeta1 * ( zeta1 + zeta2 + zeta4 - 1.0 ) ) / 2.0                                               //
                           + ( 9.0 * zeta1 * ( 3.0 * zeta1 + 3.0 * zeta2 + 3.0 * zeta4 - 2.0 ) ) / 2.0;
            adNdXi( 0, 10 ) = ( 9.0 * zeta4 * ( 6.0 * zeta1 - 1.0 ) ) / 2.0;
            adNdXi( 0, 11 ) = ( 9.0 * zeta4 * ( 3.0 * zeta4 - 1.0 ) ) / 2.0;
            adNdXi( 0, 12 ) = 0.0;
            adNdXi( 0, 13 ) = 0.0;
            adNdXi( 0, 14 ) = ( 9.0 * zeta4 * ( 6.0 * zeta1 + 6.0 * zeta2 + 6.0 * zeta4 - 5.0 ) ) / 2.0;
            adNdXi( 0, 15 ) = -( 9.0 * zeta4 * ( 3.0 * zeta4 - 1.0 ) ) / 2.0;
            adNdXi( 0, 16 ) = -27.0 * zeta2 * ( 2.0 * zeta1 + zeta2 + zeta4 - 1.0 );
            adNdXi( 0, 17 ) = 27.0 * zeta2 * zeta4;
            adNdXi( 0, 18 ) = -27.0 * zeta2 * zeta4;
            adNdXi( 0, 19 ) = -27.0 * zeta4 * ( 2.0 * zeta1 + zeta2 + zeta4 - 1.0 );

            adNdXi( 1, 0 ) = 0.0;
            adNdXi( 1, 1 ) = ( 27.0 * zeta2 * zeta2 ) / 2.0 - 9.0 * zeta2 + 1.0;
            adNdXi( 1, 2 ) = 18.0 * zeta1 + 18.0 * zeta2 + 18.0 * zeta4 - 27.0 * zeta1 * zeta2 - 27.0 * zeta1 * zeta4 - 27.0 * zeta2 * zeta4    //
                           - ( 27.0 * zeta1 * zeta1 ) / 2.0 - ( 27.0 * zeta2 * zeta2 ) / 2.0 - ( 27.0 * zeta4 * zeta4 ) / 2.0 - 11.0 / 2.0;
            adNdXi( 1, 3 ) = 0.0;
            adNdXi( 1, 4 ) = ( 9.0 * zeta1 * ( 3.0 * zeta1 - 1.0 ) ) / 2.0;
            adNdXi( 1, 5 ) = ( 9.0 * zeta1 * ( 6.0 * zeta2 - 1.0 ) ) / 2.0;
            adNdXi( 1, 6 ) = ( 9.0 * zeta1 ) / 2.0 + 36.0 * zeta2 + ( 9.0 * zeta4 ) / 2.0 - 27.0 * zeta1 * zeta2 - 27.0 * zeta2 * zeta4    //
                           - ( 81.0 * zeta2 * zeta2 ) / 2.0 - 9.0 / 2.0;
            adNdXi( 1, 7 ) = ( 9.0 * ( zeta1 + zeta2 + zeta4 - 1.0 ) * ( 3.0 * zeta1 + 3.0 * zeta2 + 3.0 * zeta4 - 2.0 ) ) / 2.0    //
                           + ( 27.0 * zeta2 * ( zeta1 + zeta2 + zeta4 - 1.0 ) ) / 2.0                                               //
                           + ( 9.0 * zeta2 * ( 3.0 * zeta1 + 3.0 * zeta2 + 3.0 * zeta4 - 2.0 ) ) / 2.0;
            adNdXi( 1, 8 )  = -( 9.0 * zeta1 * ( 3.0 * zeta1 - 1.0 ) ) / 2.0;
            adNdXi( 1, 9 )  = ( 9.0 * zeta1 * ( 6.0 * zeta1 + 6.0 * zeta2 + 6.0 * zeta4 - 5.0 ) ) / 2.0;
            adNdXi( 1, 10 ) = 0.0;
            adNdXi( 1, 11 ) = 0.0;
            adNdXi( 1, 12 ) = ( 9.0 * zeta4 * ( 6.0 * zeta2 - 1.0 ) ) / 2.0;
            adNdXi( 1, 13 ) = ( 9.0 * zeta4 * ( 3.0 * zeta4 - 1.0 ) ) / 2.0;
            adNdXi( 1, 14 ) = ( 9.0 * zeta4 * ( 6.0 * zeta1 + 6.0 * zeta2 + 6.0 * zeta4 - 5.0 ) ) / 2.0;
            adNdXi( 1, 15 ) = -( 9.0 * zeta4 * ( 3.0 * zeta4 - 1.0 ) ) / 2.0;
            adNdXi( 1, 16 ) = -27.0 * zeta1 * ( zeta1 + 2.0 * zeta2 + zeta4 - 1.0 );
            adNdXi( 1, 17 ) = 27.0 * zeta1 * zeta4;
            adNdXi( 1, 18 ) = -27.0 * zeta4 * ( zeta1 + 2.0 * zeta2 + zeta4 - 1.0 );
            adNdXi( 1, 19 ) = -27.0 * zeta1 * zeta4;

            adNdXi( 2, 0 ) = 0.0;
            adNdXi( 2, 1 ) = 0.0;
            adNdXi( 2, 2 ) = 18.0 * zeta1 + 18.0 * zeta2 + 18.0 * zeta4 - 27.0 * zeta1 * zeta2 - 27.0 * zeta1 * zeta4 - 27.0 * zeta2 * zeta4    //
                           - ( 27.0 * zeta1 * zeta1 ) / 2.0 - ( 27.0 * zeta2 * zeta2 ) / 2.0 - ( 27.0 * zeta4 * zeta4 ) / 2.0 - 11.0 / 2.0;
            adNdXi( 2, 3 )  = ( 27.0 * zeta4 * zeta4 ) / 2.0 - 9.0 * zeta4 + 1.0;
            adNdXi( 2, 4 )  = 0.0;
            adNdXi( 2, 5 )  = 0.0;
            adNdXi( 2, 6 )  = -( 9.0 * zeta2 * ( 3.0 * zeta2 - 1.0 ) ) / 2.0;
            adNdXi( 2, 7 )  = ( 9.0 * zeta2 * ( 6.0 * zeta1 + 6.0 * zeta2 + 6.0 * zeta4 - 5.0 ) ) / 2.0;
            adNdXi( 2, 8 )  = -( 9.0 * zeta1 * ( 3.0 * zeta1 - 1.0 ) ) / 2.0;
            adNdXi( 2, 9 )  = ( 9.0 * zeta1 * ( 6.0 * zeta1 + 6.0 * zeta2 + 6.0 * zeta4 - 5.0 ) ) / 2.0;
            adNdXi( 2, 10 ) = ( 9.0 * zeta1 * ( 3.0 * zeta1 - 1.0 ) ) / 2.0;
            adNdXi( 2, 11 ) = ( 9.0 * zeta1 * ( 6.0 * zeta4 - 1.0 ) ) / 2.0;
            adNdXi( 2, 12 ) = ( 9.0 * zeta2 * ( 3.0 * zeta2 - 1.0 ) ) / 2.0;
            adNdXi( 2, 13 ) = ( 9.0 * zeta2 * ( 6.0 * zeta4 - 1.0 ) ) / 2.0;
            adNdXi( 2, 14 ) = ( 27.0 * zeta1 * zeta1 ) / 2.0 + 27.0 * zeta1 * zeta2 + 54.0 * zeta1 * zeta4 - ( 45.0 * zeta1 ) / 2.0    //
                            + ( 27.0 * zeta2 * zeta2 ) / 2.0                                                                           //
                            + 54.0 * zeta2 * zeta4 - ( 45.0 * zeta2 ) / 2.0 + ( 81.0 * zeta4 * zeta4 ) / 2.0 - 45.0 * zeta4 + 9.0;
            adNdXi( 2, 15 ) = ( 9.0 * zeta1 ) / 2.0 + ( 9.0 * zeta2 ) / 2.0 + 36.0 * zeta4 - 27.0 * zeta1 * zeta4 - 27.0 * zeta2 * zeta4    //
                            - ( 81.0 * zeta4 * zeta4 ) / 2.0 - 9.0 / 2.0;
            adNdXi( 2, 16 ) = -27.0 * zeta1 * zeta2;
            adNdXi( 2, 17 ) = 27.0 * zeta1 * zeta2;
            adNdXi( 2, 18 ) = -27.0 * zeta2 * ( zeta1 + zeta2 + 2.0 * zeta4 - 1.0 );
            adNdXi( 2, 19 ) = -27.0 * zeta1 * ( zeta1 + zeta2 + 2.0 * zeta4 - 1.0 );
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_d2NdXi2(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad2NdXi2 ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TET20 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack the tet coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta4 = aXi( 2 );

            // populate output matrix
            ad2NdXi2.set_size( 6, 20, 0.0 );

            ad2NdXi2( 0, 0 )  = 27.0 * zeta1 - 9.0;
            ad2NdXi2( 0, 2 )  = 18.0 - 27.0 * zeta2 - 27.0 * zeta4 - 27.0 * zeta1;
            ad2NdXi2( 0, 4 )  = 27.0 * zeta2;
            ad2NdXi2( 0, 7 )  = 27.0 * zeta2;
            ad2NdXi2( 0, 8 )  = 36.0 - 27.0 * zeta2 - 27.0 * zeta4 - 81.0 * zeta1;
            ad2NdXi2( 0, 9 )  = 81.0 * zeta1 + 54.0 * zeta2 + 54.0 * zeta4 - 45.0;
            ad2NdXi2( 0, 10 ) = 27.0 * zeta4;
            ad2NdXi2( 0, 14 ) = 27.0 * zeta4;
            ad2NdXi2( 0, 16 ) = -54.0 * zeta2;
            ad2NdXi2( 0, 19 ) = -54.0 * zeta4;

            ad2NdXi2( 1, 1 )  = 27.0 * zeta2 - 9.0;
            ad2NdXi2( 1, 2 )  = 18.0 - 27.0 * zeta2 - 27.0 * zeta4 - 27.0 * zeta1;
            ad2NdXi2( 1, 5 )  = 27.0 * zeta1;
            ad2NdXi2( 1, 6 )  = 36.0 - 81.0 * zeta2 - 27.0 * zeta4 - 27.0 * zeta1;
            ad2NdXi2( 1, 7 )  = 54.0 * zeta1 + 81.0 * zeta2 + 54.0 * zeta4 - 45.0;
            ad2NdXi2( 1, 9 )  = 27.0 * zeta1;
            ad2NdXi2( 1, 12 ) = 27.0 * zeta4;
            ad2NdXi2( 1, 14 ) = 27.0 * zeta4;
            ad2NdXi2( 1, 16 ) = -54.0 * zeta1;
            ad2NdXi2( 1, 18 ) = -54.0 * zeta4;

            ad2NdXi2( 2, 2 )  = 18.0 - 27.0 * zeta2 - 27.0 * zeta4 - 27.0 * zeta1;
            ad2NdXi2( 2, 3 )  = 27.0 * zeta4 - 9.0;
            ad2NdXi2( 2, 7 )  = 27.0 * zeta2;
            ad2NdXi2( 2, 9 )  = 27.0 * zeta1;
            ad2NdXi2( 2, 11 ) = 27.0 * zeta1;
            ad2NdXi2( 2, 13 ) = 27.0 * zeta2;
            ad2NdXi2( 2, 14 ) = 54.0 * zeta1 + 54.0 * zeta2 + 81.0 * zeta4 - 45.0;
            ad2NdXi2( 2, 15 ) = 36.0 - 27.0 * zeta2 - 81.0 * zeta4 - 27.0 * zeta1;
            ad2NdXi2( 2, 18 ) = -54.0 * zeta2;
            ad2NdXi2( 2, 19 ) = -54.0 * zeta1;

            ad2NdXi2( 3, 2 )  = 18.0 - 27.0 * zeta2 - 27.0 * zeta4 - 27.0 * zeta1;
            ad2NdXi2( 3, 6 )  = 9.0 / 2.0 - 27.0 * zeta2;
            ad2NdXi2( 3, 7 )  = 27.0 * zeta1 + 54.0 * zeta2 + 27.0 * zeta4 - 45.0 / 2.0;
            ad2NdXi2( 3, 9 )  = 27.0 * zeta1;
            ad2NdXi2( 3, 12 ) = 27.0 * zeta2 - 9.0 / 2.0;
            ad2NdXi2( 3, 13 ) = 27.0 * zeta4 - 9.0 / 2.0;
            ad2NdXi2( 3, 14 ) = 27.0 * zeta1 + 27.0 * zeta2 + 54.0 * zeta4 - 45.0 / 2.0;
            ad2NdXi2( 3, 15 ) = 9.0 / 2.0 - 27.0 * zeta4;
            ad2NdXi2( 3, 16 ) = -27.0 * zeta1;
            ad2NdXi2( 3, 17 ) = 27.0 * zeta1;
            ad2NdXi2( 3, 18 ) = 27.0 - 54.0 * zeta2 - 54.0 * zeta4 - 27.0 * zeta1;
            ad2NdXi2( 3, 19 ) = -27.0 * zeta1;

            ad2NdXi2( 4, 2 )  = 18.0 - 27.0 * zeta2 - 27.0 * zeta4 - 27.0 * zeta1;
            ad2NdXi2( 4, 7 )  = 27.0 * zeta2;
            ad2NdXi2( 4, 8 )  = 9.0 / 2.0 - 27.0 * zeta1;
            ad2NdXi2( 4, 9 )  = 54.0 * zeta1 + 27.0 * zeta2 + 27.0 * zeta4 - 45.0 / 2.0;
            ad2NdXi2( 4, 10 ) = 27.0 * zeta1 - 9.0 / 2.0;
            ad2NdXi2( 4, 11 ) = 27.0 * zeta4 - 9.0 / 2.0;
            ad2NdXi2( 4, 14 ) = 27.0 * zeta1 + 27.0 * zeta2 + 54.0 * zeta4 - 45.0 / 2.0;
            ad2NdXi2( 4, 15 ) = 9.0 / 2.0 - 27.0 * zeta4;
            ad2NdXi2( 4, 16 ) = -27.0 * zeta2;
            ad2NdXi2( 4, 17 ) = 27.0 * zeta2;
            ad2NdXi2( 4, 18 ) = -27.0 * zeta2;
            ad2NdXi2( 4, 19 ) = 27.0 - 27.0 * zeta2 - 54.0 * zeta4 - 54.0 * zeta1;

            ad2NdXi2( 5, 2 )  = 18.0 - 27.0 * zeta2 - 27.0 * zeta4 - 27.0 * zeta1;
            ad2NdXi2( 5, 4 )  = 27.0 * zeta1 - 9.0 / 2.0;
            ad2NdXi2( 5, 5 )  = 27.0 * zeta2 - 9.0 / 2.0;
            ad2NdXi2( 5, 6 )  = 9.0 / 2.0 - 27.0 * zeta2;
            ad2NdXi2( 5, 7 )  = 27.0 * zeta1 + 54.0 * zeta2 + 27.0 * zeta4 - 45.0 / 2.0;
            ad2NdXi2( 5, 8 )  = 9.0 / 2.0 - 27.0 * zeta1;
            ad2NdXi2( 5, 9 )  = 54.0 * zeta1 + 27.0 * zeta2 + 27.0 * zeta4 - 45.0 / 2.0;
            ad2NdXi2( 5, 14 ) = 27.0 * zeta4;
            ad2NdXi2( 5, 16 ) = 27.0 - 54.0 * zeta2 - 27.0 * zeta4 - 54.0 * zeta1;
            ad2NdXi2( 5, 17 ) = 27.0 * zeta4;
            ad2NdXi2( 5, 18 ) = -27.0 * zeta4;
            ad2NdXi2( 5, 19 ) = -27.0 * zeta4;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_d3NdXi3(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad3NdXi3 ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( false, "TET20 - eval_d3NdXi3: 3rd order derivatives not implemented for this element." );

            ad3NdXi3.set_size( 1, 20, 0.0 );
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET20_HPP_ */
