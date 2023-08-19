/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Tri3.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TRI3_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TRI3_HPP_

#include "assert.h"
#include "typedefs.hpp"                         //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 3 >::get_number_of_param_dimensions() const
        {
            return 2;
        }

        //------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 3 >::get_interpolation_order() const
        {
            return Interpolation_Order::LINEAR;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 3 >::get_param_coords( Matrix< DDRMat >& aXiHat ) const
        {
            aXiHat = {
                { 1.0, 0.0, 0.0 },
                { 0.0, 1.0, 0.0 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 3 >::eval_N(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.numel() >= 2, "TRI3 - eval_N: aXi not allocated or hat wrong size." );

            // get the triangular coordinates
            const real zeta1 = aXi( 0 );
            const real zeta2 = aXi( 1 );
            const real zeta3 = 1.0 - aXi( 0 ) - aXi( 1 );

            // populate matrix with values
            aNXi.set_size( 1, 3 );
            aNXi( 0 ) = zeta1;
            aNXi( 1 ) = zeta2;
            aNXi( 2 ) = zeta3;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 3 >::eval_dNdXi(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.numel() >= 2, "TRI3 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // populate output matrix
            adNdXi.set_size( 2, 3 );
            adNdXi( 0, 0 ) = 1.0;
            adNdXi( 0, 1 ) = 0.0;
            adNdXi( 0, 2 ) = -1.0;

            adNdXi( 1, 0 ) = 0.0;
            adNdXi( 1, 1 ) = 1.0;
            adNdXi( 1, 2 ) = -1.0;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 3 >::eval_d2NdXi2(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "TRI3 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // populate output matrix
            ad2NdXi2.set_size( 3, 3, 0.0 );
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 3 >::eval_d3NdXi3(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "TRI3 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            ad3NdXi3.set_size( 4, 3, 0.0 );
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TRI3_HPP_ */

