/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Quad8.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD8_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD8_HPP_

#include "assert.h"
#include "moris_typedefs.hpp"                   //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src

// checked against femdoc:
// - N
// - adNdxi
// - d2Ndxi ( fixed bug in femdoc)

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    template<>
    uint
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 8 >::get_number_of_param_dimensions() const
    {
        return 2;
    }

    //------------------------------------------------------------------------------

    template<>
    Interpolation_Order
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 8 >::get_interpolation_order() const
    {
        return Interpolation_Order::SERENDIPITY;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 8 >::get_param_coords( Matrix< DDRMat > &aXiHat ) const
    {
        aXiHat.set_size( 2, 8, 0.0 );
        aXiHat( 0, 0 ) = -1.000000;
        aXiHat( 1, 0 ) = -1.000000;
        aXiHat( 0, 1 ) = 1.000000;
        aXiHat( 1, 1 ) = -1.000000;
        aXiHat( 0, 2 ) = 1.000000;
        aXiHat( 1, 2 ) = 1.000000;
        aXiHat( 0, 3 ) = -1.000000;
        aXiHat( 1, 3 ) = 1.000000;
        aXiHat( 0, 4 ) = 0.000000;
        aXiHat( 1, 4 ) = -1.000000;
        aXiHat( 0, 5 ) = 1.000000;
        aXiHat( 1, 5 ) = 0.000000;
        aXiHat( 0, 6 ) = 0.000000;
        aXiHat( 1, 6 ) = 1.000000;
        aXiHat( 0, 7 ) = -1.000000;
        aXiHat( 1, 7 ) = 0.000000;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 8 >::eval_N( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                          &aNXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 2,
                "QUAD8 - eval_N: aXi not allocated or hat wrong size." );

        // unpack xi and eta from input vector
        real xi  = aXi( 0 );
        real eta = aXi( 1 );

        // often used constants
        real xi2  = std::pow( xi, 2 );
        real eta2 = std::pow( eta, 2 );

        // populate output matrix
        aNXi.set_size( 1, 8 );
        aNXi( 0 ) = -( eta - 1.0 ) * ( xi - 1.0 ) * ( eta + xi + 1.0 ) * 0.25;
        aNXi( 1 ) = ( eta - 1.0 ) * ( xi + 1.0 ) * ( eta - xi + 1.0 ) * 0.25;
        aNXi( 2 ) = ( eta + 1.0 ) * ( xi + 1.0 ) * ( eta + xi - 1.0 ) * 0.25;
        aNXi( 3 ) = ( eta + 1.0 ) * ( xi - 1.0 ) * ( -eta + xi + 1.0 ) * 0.25;
        aNXi( 4 ) = ( xi2 - 1.0 ) * ( eta - 1.0 ) * 0.5;
        aNXi( 5 ) = -( eta2 - 1.0 ) * ( xi + 1.0 ) * 0.5;
        aNXi( 6 ) = -( xi2 - 1.0 ) * ( eta + 1.0 ) * 0.5;
        aNXi( 7 ) = ( eta2 - 1.0 ) * ( xi - 1.0 ) * 0.5;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 8 >::eval_dNdXi( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                              &adNdXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 2, "QUAD8 - eval_dNdXi: aXi not allocated or hat wrong size." );

        // unpack xi and eta from input vector
        real xi  = aXi( 0 );
        real eta = aXi( 1 );

        // often used constants
        real xi2  = std::pow( xi, 2 );
        real eta2 = std::pow( eta, 2 );

        // populate output matrix
        adNdXi.set_size( 2, 8 );
        adNdXi( 0, 0 ) = -( eta + xi * 2.0 ) * ( eta - 1.0 ) * 0.25;
        adNdXi( 1, 0 ) = -( eta * 2.0 + xi ) * ( xi - 1.0 ) * 0.25;

        adNdXi( 0, 1 ) = ( eta - xi * 2.0 ) * ( eta - 1.0 ) * 0.25;
        adNdXi( 1, 1 ) = ( xi + 1.0 ) * ( eta * 2.0 - xi ) * 0.25;

        adNdXi( 0, 2 ) = ( eta + xi * 2.0 ) * ( eta + 1.0 ) * 0.25;
        adNdXi( 1, 2 ) = ( eta * 2.0 + xi ) * ( xi + 1.0 ) * 0.25;

        adNdXi( 0, 3 ) = -( eta - xi * 2.0 ) * ( eta + 1.0 ) * 0.25;
        adNdXi( 1, 3 ) = -( xi - 1.0 ) * ( eta * 2.0 - xi ) * 0.25;

        adNdXi( 0, 4 ) = xi * ( eta - 1.0 );
        adNdXi( 1, 4 ) = 0.5 * xi2 - 0.5;

        adNdXi( 0, 5 ) = -0.5 * eta2 + 0.5;
        adNdXi( 1, 5 ) = -eta * ( xi + 1.0 );

        adNdXi( 0, 6 ) = -xi * ( eta + 1.0 );
        adNdXi( 1, 6 ) = -0.5 * xi2 + 0.5;

        adNdXi( 0, 7 ) = 0.5 * eta2 - 0.5;
        adNdXi( 1, 7 ) = eta * ( xi - 1.0 );
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 8 >::eval_d2NdXi2( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                                &ad2NdXi2 ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 2, "QUAD8 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

        // unpack xi and eta from input vector
        real xi  = aXi( 0 );
        real eta = aXi( 1 );

        // populate output matrix
        ad2NdXi2.set_size( 3, 8 );
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
        ad2NdXi2( 1, 5 ) = -xi - 1.0;
        ad2NdXi2( 2, 5 ) = -eta;

        ad2NdXi2( 0, 6 ) = -eta - 1.0;
        ad2NdXi2( 1, 6 ) = 0.0;
        ad2NdXi2( 2, 6 ) = -xi;

        ad2NdXi2( 0, 7 ) = 0.0;
        ad2NdXi2( 1, 7 ) = xi - 1.0;
        ad2NdXi2( 2, 7 ) = eta;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 8 >::eval_d3NdXi3( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                                &ad3NdXi3 ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( false, "QUAD8 - eval_d3NdXi3: 3rd order derivatives not implemented for this element." );

        ad3NdXi3.set_size( 4, 8, 0.0 );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD8_HPP_ */
