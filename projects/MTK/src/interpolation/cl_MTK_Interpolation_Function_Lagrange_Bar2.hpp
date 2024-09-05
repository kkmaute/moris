/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Bar2.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_BAR2_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_BAR2_HPP_

#include "assert.h"
#include "moris_typedefs.hpp"                   //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    template<>
    uint
    Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 2 >::get_number_of_param_dimensions() const
    {
        return 1;
    }

    //------------------------------------------------------------------------------

    template<>
    Interpolation_Order
    Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 2 >::get_interpolation_order() const
    {
        return Interpolation_Order::LINEAR;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 2 >::get_param_coords( Matrix< DDRMat > &aXiHat ) const
    {
        aXiHat.set_size( 1, 2, 0.0 );
        aXiHat( 0 ) = -1.000000;
        aXiHat( 1 ) = 1.000000;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 2 >::eval_N( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                          &aNXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 1, "LINE2 - eval_N: aXi not allocated or hat wrong size." );

        // set aNXi
        real xi = aXi( 0 );
        aNXi.set_size( 1, 2 );
        aNXi( 0 ) = 0.5 * ( 1.0 - xi );
        aNXi( 1 ) = 0.5 * ( 1.0 + xi );
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 2 >::eval_dNdXi( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                              &adNdXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 1, "LINE2 - eval_dNdXi: aXi not allocated or hat wrong size." );

        // set adNdXi
        adNdXi.set_size( 1, 2 );
        adNdXi( 0 ) = -0.5;
        adNdXi( 1 ) = 0.5;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 2 >::eval_d2NdXi2( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                                &ad2NdXi2 ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 1, "LINE2 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

        ad2NdXi2.set_size( 1, 2, 0.0 );
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 2 >::eval_d3NdXi3( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                                &ad3NdXi3 ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 1, "LINE2 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

        ad3NdXi3.set_size( 1, 2, 0.0 );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::mtk

//------------------------------------------------------------------------------
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_BAR2_HPP_ */
