/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Bar4.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_BAR4_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_BAR4_HPP_

#include "assert.h"
#include "moris_typedefs.hpp"                   //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src
namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::get_number_of_param_dimensions() const
        {
            return 1;
        }

        //------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::get_interpolation_order() const
        {
            return Interpolation_Order::CUBIC;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::get_param_coords( Matrix< DDRMat > &aXiHat ) const
        {
            aXiHat.set_size( 1, 4, 0.0 );
            aXiHat( 0 ) = -1.000000;
            aXiHat( 1 ) = 1.000000;
            aXiHat( 2 ) = -1.0 / 3.0;
            aXiHat( 3 ) = 1.0 / 3.0;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::eval_N( const Matrix< DDRMat > &aXi,
                Matrix< DDRMat >                                                                                          &aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "LINE4 - eval_N: aXi not allocated or hat wrong size." );

            real xi   = aXi( 0 );
            real t116 = 1.0 / 16.0;

            aNXi.set_size( 1, 4, 0.0 );
            aNXi( 0 ) = t116 * ( 3.0 * xi + 1.0 ) * ( 3.0 * xi - 1.0 ) * ( 1.0 - xi );
            aNXi( 1 ) = t116 * ( 3.0 * xi + 1.0 ) * ( 3.0 * xi - 1.0 ) * ( 1.0 + xi );
            aNXi( 2 ) = 9.0 * t116 * ( xi + 1.0 ) * ( 3.0 * xi - 1.0 ) * ( xi - 1.0 );
            aNXi( 3 ) = 9.0 * t116 * ( xi + 1.0 ) * ( 3.0 * xi + 1.0 ) * ( 1.0 - xi );
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::eval_dNdXi( const Matrix< DDRMat > &aXi,
                Matrix< DDRMat >                                                                                              &adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "LINE4 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack param point
            real xi   = aXi( 0 );
            real xi2  = std::pow( xi, 2 );
            real t116 = 1.0 / 16.0;

            // set adNdXi
            adNdXi.set_size( 1, 4 );
            adNdXi( 0 ) = t116 * ( -27.0 * xi2 + 18.0 * xi + 1.0 );
            adNdXi( 1 ) = t116 * ( 27.0 * xi2 + 18.0 * xi - 1.0 );
            adNdXi( 2 ) = 9.0 * t116 * ( 9.0 * xi2 - 2.0 * xi - 3.0 );
            adNdXi( 3 ) = 9.0 * t116 * ( -9.0 * xi2 - 2.0 * xi + 3.0 );
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::eval_d2NdXi2( const Matrix< DDRMat > &aXi,
                Matrix< DDRMat >                                                                                                &ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "LINE4 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            real xi  = aXi( 0 );
            real t98 = 9.0 / 8.0;

            ad2NdXi2.set_size( 1, 4 );
            ad2NdXi2( 0 ) = t98 * ( 1.0 - 3.0 * xi );
            ad2NdXi2( 1 ) = t98 * ( 1.0 + 3.0 * xi );
            ad2NdXi2( 2 ) = t98 * ( 9.0 * xi - 1.0 );
            ad2NdXi2( 3 ) = -t98 * ( 9.0 * xi + 1.0 );
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::eval_d3NdXi3( const Matrix< DDRMat > &aXi,
                Matrix< DDRMat >                                                                                                &ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "LINE4 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            real t98 = 9.0 / 8.0;

            ad3NdXi3.set_size( 1, 4 );
            ad3NdXi3( 0 ) = -t98 * 3.0;
            ad3NdXi3( 1 ) = t98 * 3.0;
            ad3NdXi3( 2 ) = t98 * 9.0;
            ad3NdXi3( 3 ) = -t98 * 9.0;
        }
        //------------------------------------------------------------------------------

    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_BAR4_HPP_ */
