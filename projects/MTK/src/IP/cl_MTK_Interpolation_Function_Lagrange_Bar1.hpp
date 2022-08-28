/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Bar1.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_

#include "assert.h"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_MTK_Interpolation_Function.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1 >::get_number_of_param_dimensions() const
        {
            return 1;
        }

//------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1 >::get_interpolation_order() const
        {
            return Interpolation_Order::CONSTANT;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1 >::get_param_coords( Matrix< DDRMat > & aXiHat ) const
        {
            aXiHat.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1  >::eval_N( const Matrix< DDRMat > & aXi,
                                                                                                               Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "LINE1 - eval_N: aXi not allocated or hat wrong size." );

            aNXi.set_size( 1, 1, 1.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1 >::eval_dNdXi( const Matrix< DDRMat > & aXi,
                                                                                                                  Matrix< DDRMat > & adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "LINE1 - eval_dNdXi: aXi not allocated or hat wrong size." );

            adNdXi.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1  >::eval_d2NdXi2( const Matrix< DDRMat > & aXi,
                                                                                                                     Matrix< DDRMat > & ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "LINE1 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            ad2NdXi2.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1  >::eval_d3NdXi3( const Matrix< DDRMat > & aXi,
                                                                                                                     Matrix< DDRMat > & ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "LINE1 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            ad3NdXi3.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

//------------------------------------------------------------------------------
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_ */
