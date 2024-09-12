/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Tet4.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET4_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET4_HPP_

#include "assert.h"
#include "moris_typedefs.hpp"                   //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    template<>
    uint
    Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 4 >::get_number_of_param_dimensions() const
    {
        return 3;
    }

    //------------------------------------------------------------------------------

    template<>
    Interpolation_Order
    Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 4 >::get_interpolation_order() const
    {
        return Interpolation_Order::LINEAR;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 4 >::get_param_coords( Matrix< DDRMat >& aXiHat ) const
    {
        aXiHat = {
            { 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000 }
        };
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 4 >::eval_N(
            const Matrix< DDRMat >& aXi,
            Matrix< DDRMat >&       aNXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 3, "TET4 - eval_N: aXi not allocated or hat wrong size." );

        // unpack zeta1, zeta2, and  zeta3 from input vector
        real zeta1 = aXi( 0 );
        real zeta2 = aXi( 1 );
        real zeta3 = 1.0 - aXi( 0 ) - aXi( 1 ) - aXi( 2 );
        real zeta4 = aXi( 2 );

        // populate matrix with values
        aNXi.set_size( 1, 4 );
        aNXi( 0 ) = zeta1;
        aNXi( 1 ) = zeta2;
        aNXi( 2 ) = zeta3;
        aNXi( 3 ) = zeta4;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 4 >::eval_dNdXi(
            const Matrix< DDRMat >& aXi,
            Matrix< DDRMat >&       adNdXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 3, "TET4 - eval_dNdXi: aXi not allocated or hat wrong size." );

        // populate output matrix
        adNdXi.set_size( 3, 4 );
        adNdXi( 0, 0 ) = 1.0;
        adNdXi( 0, 1 ) = 0.0;
        adNdXi( 0, 2 ) = -1.0;
        adNdXi( 0, 3 ) = 0.0;

        adNdXi( 1, 0 ) = 0.0;
        adNdXi( 1, 1 ) = 1.0;
        adNdXi( 1, 2 ) = -1.0;
        adNdXi( 1, 3 ) = 0.0;

        adNdXi( 2, 0 ) = 0.0;
        adNdXi( 2, 1 ) = 0.0;
        adNdXi( 2, 2 ) = -1.0;
        adNdXi( 2, 3 ) = 1.0;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 4 >::eval_d2NdXi2(
            const Matrix< DDRMat >& aXi,
            Matrix< DDRMat >&       ad2NdXi2 ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 3, "TET4 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

        // populate output matrix
        ad2NdXi2.set_size( 6, 4, 0.0 );
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 4 >::eval_d3NdXi3(
            const Matrix< DDRMat >& aXi,
            Matrix< DDRMat >&       ad3NdXi3 ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( false, "TET4 - eval_d3NdXi3: 3rd order derivatives not implemented for this element." );

        ad3NdXi3.set_size( 1, 4, 0.0 );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET4_HPP_ */
