/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Lagrange_Quad25.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD25_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD25_HPP_

#include "assert.h"
#include "moris_typedefs.hpp"                   //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    template<>
    uint
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 25 >::get_number_of_param_dimensions() const
    {
        return 2;
    }

    //------------------------------------------------------------------------------

    template<>
    Interpolation_Order
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 25 >::get_interpolation_order() const
    {
        return Interpolation_Order::QUARTIC;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 25 >::get_param_coords( Matrix< DDRMat > &aXiHat ) const
    {
        aXiHat.set_size( 2, 25, 0.0 );
        real c = 0.5;

        // Node layout (xi, eta), row-major ordering
        aXiHat( 0, 0 ) = -1.0;
        aXiHat( 1, 0 ) = -1.0;
        aXiHat( 0, 1 ) = -c;
        aXiHat( 1, 1 ) = -1.0;
        aXiHat( 0, 2 ) = 0.0;
        aXiHat( 1, 2 ) = -1.0;
        aXiHat( 0, 3 ) = c;
        aXiHat( 1, 3 ) = -1.0;
        aXiHat( 0, 4 ) = 1.0;
        aXiHat( 1, 4 ) = -1.0;

        aXiHat( 0, 5 ) = -1.0;
        aXiHat( 1, 5 ) = -c;
        aXiHat( 0, 6 ) = -c;
        aXiHat( 1, 6 ) = -c;
        aXiHat( 0, 7 ) = 0.0;
        aXiHat( 1, 7 ) = -c;
        aXiHat( 0, 8 ) = c;
        aXiHat( 1, 8 ) = -c;
        aXiHat( 0, 9 ) = 1.0;
        aXiHat( 1, 9 ) = -c;

        aXiHat( 0, 10 ) = -1.0;
        aXiHat( 1, 10 ) = 0.0;
        aXiHat( 0, 11 ) = -c;
        aXiHat( 1, 11 ) = 0.0;
        aXiHat( 0, 12 ) = 0.0;
        aXiHat( 1, 12 ) = 0.0;
        aXiHat( 0, 13 ) = c;
        aXiHat( 1, 13 ) = 0.0;
        aXiHat( 0, 14 ) = 1.0;
        aXiHat( 1, 14 ) = 0.0;

        aXiHat( 0, 15 ) = -1.0;
        aXiHat( 1, 15 ) = c;
        aXiHat( 0, 16 ) = -c;
        aXiHat( 1, 16 ) = c;
        aXiHat( 0, 17 ) = 0.0;
        aXiHat( 1, 17 ) = c;
        aXiHat( 0, 18 ) = c;
        aXiHat( 1, 18 ) = c;
        aXiHat( 0, 19 ) = 1.0;
        aXiHat( 1, 19 ) = c;

        aXiHat( 0, 20 ) = -1.0;
        aXiHat( 1, 20 ) = 1.0;
        aXiHat( 0, 21 ) = -c;
        aXiHat( 1, 21 ) = 1.0;
        aXiHat( 0, 22 ) = 0.0;
        aXiHat( 1, 22 ) = 1.0;
        aXiHat( 0, 23 ) = c;
        aXiHat( 1, 23 ) = 1.0;
        aXiHat( 0, 24 ) = 1.0;
        aXiHat( 1, 24 ) = 1.0;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 25 >::eval_N( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                           &aNXi ) const
    {
        MORIS_ASSERT( aXi.length() >= 2, "QUAD25 - eval_N: aXi not allocated or has wrong size." );

        real xi  = aXi( 0 );
        real eta = aXi( 1 );

        real a0 = ( 2.0 / 3.0 ) * ( xi + 0.5 ) * xi * ( xi - 0.5 ) * ( xi - 1.0 );
        real a1 = -( 8.0 / 3.0 ) * ( xi + 1.0 ) * xi * ( xi - 0.5 ) * ( xi - 1.0 );
        real a2 = 4.0 * ( xi + 1.0 ) * ( xi + 0.5 ) * ( xi - 0.5 ) * ( xi - 1.0 );
        real a3 = -( 8.0 / 3.0 ) * ( xi + 1.0 ) * ( xi + 0.5 ) * xi * ( xi - 1.0 );
        real a4 = ( 2.0 / 3.0 ) * ( xi + 1.0 ) * ( xi + 0.5 ) * xi * ( xi - 0.5 );

        real b0 = ( 2.0 / 3.0 ) * ( eta + 0.5 ) * eta * ( eta - 0.5 ) * ( eta - 1.0 );
        real b1 = -( 8.0 / 3.0 ) * ( eta + 1.0 ) * eta * ( eta - 0.5 ) * ( eta - 1.0 );
        real b2 = 4.0 * ( eta + 1.0 ) * ( eta + 0.5 ) * ( eta - 0.5 ) * ( eta - 1.0 );
        real b3 = -( 8.0 / 3.0 ) * ( eta + 1.0 ) * ( eta + 0.5 ) * eta * ( eta - 1.0 );
        real b4 = ( 2.0 / 3.0 ) * ( eta + 1.0 ) * ( eta + 0.5 ) * eta * ( eta - 0.5 );

        aNXi.set_size( 1, 25 );
        aNXi( 0 )  = a0 * b0;
        aNXi( 1 )  = a1 * b0;
        aNXi( 2 )  = a2 * b0;
        aNXi( 3 )  = a3 * b0;
        aNXi( 4 )  = a4 * b0;
        aNXi( 5 )  = a0 * b1;
        aNXi( 6 )  = a1 * b1;
        aNXi( 7 )  = a2 * b1;
        aNXi( 8 )  = a3 * b1;
        aNXi( 9 )  = a4 * b1;
        aNXi( 10 ) = a0 * b2;
        aNXi( 11 ) = a1 * b2;
        aNXi( 12 ) = a2 * b2;
        aNXi( 13 ) = a3 * b2;
        aNXi( 14 ) = a4 * b2;
        aNXi( 15 ) = a0 * b3;
        aNXi( 16 ) = a1 * b3;
        aNXi( 17 ) = a2 * b3;
        aNXi( 18 ) = a3 * b3;
        aNXi( 19 ) = a4 * b3;
        aNXi( 20 ) = a0 * b4;
        aNXi( 21 ) = a1 * b4;
        aNXi( 22 ) = a2 * b4;
        aNXi( 23 ) = a3 * b4;
        aNXi( 24 ) = a4 * b4;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 25 >::eval_dNdXi( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                               &adNdXi ) const
    {
        MORIS_ASSERT( aXi.length() >= 2, "QUAD25 - eval_dNdXi: aXi not allocated or has wrong size." );

        real xi  = aXi( 0 );
        real eta = aXi( 1 );

        real a0 = ( 2.0 / 3.0 ) * ( xi + 0.5 ) * xi * ( xi - 0.5 ) * ( xi - 1.0 );
        real a1 = -( 8.0 / 3.0 ) * ( xi + 1.0 ) * xi * ( xi - 0.5 ) * ( xi - 1.0 );
        real a2 = 4.0 * ( xi + 1.0 ) * ( xi + 0.5 ) * ( xi - 0.5 ) * ( xi - 1.0 );
        real a3 = -( 8.0 / 3.0 ) * ( xi + 1.0 ) * ( xi + 0.5 ) * xi * ( xi - 1.0 );
        real a4 = ( 2.0 / 3.0 ) * ( xi + 1.0 ) * ( xi + 0.5 ) * xi * ( xi - 0.5 );

        real b0 = ( 2.0 / 3.0 ) * ( eta + 0.5 ) * eta * ( eta - 0.5 ) * ( eta - 1.0 );
        real b1 = -( 8.0 / 3.0 ) * ( eta + 1.0 ) * eta * ( eta - 0.5 ) * ( eta - 1.0 );
        real b2 = 4.0 * ( eta + 1.0 ) * ( eta + 0.5 ) * ( eta - 0.5 ) * ( eta - 1.0 );
        real b3 = -( 8.0 / 3.0 ) * ( eta + 1.0 ) * ( eta + 0.5 ) * eta * ( eta - 1.0 );
        real b4 = ( 2.0 / 3.0 ) * ( eta + 1.0 ) * ( eta + 0.5 ) * eta * ( eta - 0.5 );

        real da0 = ( 8.0 / 3.0 ) * xi * xi * xi - 2.0 * xi * xi - ( 1.0 / 3.0 ) * xi + 1.0 / 6.0;
        real da1 = -( 32.0 / 3.0 ) * xi * xi * xi + 12.0 * xi * xi - 4.0 / 3.0;
        real da2 = 16.0 * xi * xi * xi - 10.0 * xi;
        real da3 = -( 32.0 / 3.0 ) * xi * xi * xi - 4.0 * xi * xi + ( 16.0 / 3.0 ) * xi + 4.0 / 3.0;
        real da4 = ( 8.0 / 3.0 ) * xi * xi * xi + 2.0 * xi * xi - ( 1.0 / 3.0 ) * xi - 1.0 / 6.0;

        real db0 = ( 8.0 / 3.0 ) * eta * eta * eta - 2.0 * eta * eta - ( 1.0 / 3.0 ) * eta + 1.0 / 6.0;
        real db1 = -( 32.0 / 3.0 ) * eta * eta * eta + 12.0 * eta * eta - 4.0 / 3.0;
        real db2 = 16.0 * eta * eta * eta - 10.0 * eta;
        real db3 = -( 32.0 / 3.0 ) * eta * eta * eta - 4.0 * eta * eta + ( 16.0 / 3.0 ) * eta + 4.0 / 3.0;
        real db4 = ( 8.0 / 3.0 ) * eta * eta * eta + 2.0 * eta * eta - ( 1.0 / 3.0 ) * eta - 1.0 / 6.0;

        adNdXi.set_size( 2, 25 );
        adNdXi( 0, 0 ) = da0 * b0;
        adNdXi( 1, 0 ) = a0 * db0;
        adNdXi( 0, 1 ) = da1 * b0;
        adNdXi( 1, 1 ) = a1 * db0;
        adNdXi( 0, 2 ) = da2 * b0;
        adNdXi( 1, 2 ) = a2 * db0;
        adNdXi( 0, 3 ) = da3 * b0;
        adNdXi( 1, 3 ) = a3 * db0;
        adNdXi( 0, 4 ) = da4 * b0;
        adNdXi( 1, 4 ) = a4 * db0;

        adNdXi( 0, 5 ) = da0 * b1;
        adNdXi( 1, 5 ) = a0 * db1;
        adNdXi( 0, 6 ) = da1 * b1;
        adNdXi( 1, 6 ) = a1 * db1;
        adNdXi( 0, 7 ) = da2 * b1;
        adNdXi( 1, 7 ) = a2 * db1;
        adNdXi( 0, 8 ) = da3 * b1;
        adNdXi( 1, 8 ) = a3 * db1;
        adNdXi( 0, 9 ) = da4 * b1;
        adNdXi( 1, 9 ) = a4 * db1;

        adNdXi( 0, 10 ) = da0 * b2;
        adNdXi( 1, 10 ) = a0 * db2;
        adNdXi( 0, 11 ) = da1 * b2;
        adNdXi( 1, 11 ) = a1 * db2;
        adNdXi( 0, 12 ) = da2 * b2;
        adNdXi( 1, 12 ) = a2 * db2;
        adNdXi( 0, 13 ) = da3 * b2;
        adNdXi( 1, 13 ) = a3 * db2;
        adNdXi( 0, 14 ) = da4 * b2;
        adNdXi( 1, 14 ) = a4 * db2;

        adNdXi( 0, 15 ) = da0 * b3;
        adNdXi( 1, 15 ) = a0 * db3;
        adNdXi( 0, 16 ) = da1 * b3;
        adNdXi( 1, 16 ) = a1 * db3;
        adNdXi( 0, 17 ) = da2 * b3;
        adNdXi( 1, 17 ) = a2 * db3;
        adNdXi( 0, 18 ) = da3 * b3;
        adNdXi( 1, 18 ) = a3 * db3;
        adNdXi( 0, 19 ) = da4 * b3;
        adNdXi( 1, 19 ) = a4 * db3;

        adNdXi( 0, 20 ) = da0 * b4;
        adNdXi( 1, 20 ) = a0 * db4;
        adNdXi( 0, 21 ) = da1 * b4;
        adNdXi( 1, 21 ) = a1 * db4;
        adNdXi( 0, 22 ) = da2 * b4;
        adNdXi( 1, 22 ) = a2 * db4;
        adNdXi( 0, 23 ) = da3 * b4;
        adNdXi( 1, 23 ) = a3 * db4;
        adNdXi( 0, 24 ) = da4 * b4;
        adNdXi( 1, 24 ) = a4 * db4;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 25 >::eval_d2NdXi2( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                                 &ad2NdXi2 ) const
    {
        MORIS_ASSERT( aXi.length() >= 2, "QUAD25 - eval_d2NdXi2: aXi not allocated or has wrong size." );

        real xi  = aXi( 0 );
        real eta = aXi( 1 );

        real a0 = ( 2.0 / 3.0 ) * ( xi + 0.5 ) * xi * ( xi - 0.5 ) * ( xi - 1.0 );
        real a1 = -( 8.0 / 3.0 ) * ( xi + 1.0 ) * xi * ( xi - 0.5 ) * ( xi - 1.0 );
        real a2 = 4.0 * ( xi + 1.0 ) * ( xi + 0.5 ) * ( xi - 0.5 ) * ( xi - 1.0 );
        real a3 = -( 8.0 / 3.0 ) * ( xi + 1.0 ) * ( xi + 0.5 ) * xi * ( xi - 1.0 );
        real a4 = ( 2.0 / 3.0 ) * ( xi + 1.0 ) * ( xi + 0.5 ) * xi * ( xi - 0.5 );

        real b0 = ( 2.0 / 3.0 ) * ( eta + 0.5 ) * eta * ( eta - 0.5 ) * ( eta - 1.0 );
        real b1 = -( 8.0 / 3.0 ) * ( eta + 1.0 ) * eta * ( eta - 0.5 ) * ( eta - 1.0 );
        real b2 = 4.0 * ( eta + 1.0 ) * ( eta + 0.5 ) * ( eta - 0.5 ) * ( eta - 1.0 );
        real b3 = -( 8.0 / 3.0 ) * ( eta + 1.0 ) * ( eta + 0.5 ) * eta * ( eta - 1.0 );
        real b4 = ( 2.0 / 3.0 ) * ( eta + 1.0 ) * ( eta + 0.5 ) * eta * ( eta - 0.5 );

        real da0 = ( 8.0 / 3.0 ) * xi * xi * xi - 2.0 * xi * xi - ( 1.0 / 3.0 ) * xi + 1.0 / 6.0;
        real da1 = -( 32.0 / 3.0 ) * xi * xi * xi + 12.0 * xi * xi - 4.0 / 3.0;
        real da2 = 16.0 * xi * xi * xi - 10.0 * xi;
        real da3 = -( 32.0 / 3.0 ) * xi * xi * xi - 4.0 * xi * xi + ( 16.0 / 3.0 ) * xi + 4.0 / 3.0;
        real da4 = ( 8.0 / 3.0 ) * xi * xi * xi + 2.0 * xi * xi - ( 1.0 / 3.0 ) * xi - 1.0 / 6.0;

        real db0 = ( 8.0 / 3.0 ) * eta * eta * eta - 2.0 * eta * eta - ( 1.0 / 3.0 ) * eta + 1.0 / 6.0;
        real db1 = -( 32.0 / 3.0 ) * eta * eta * eta + 12.0 * eta * eta - 4.0 / 3.0;
        real db2 = 16.0 * eta * eta * eta - 10.0 * eta;
        real db3 = -( 32.0 / 3.0 ) * eta * eta * eta - 4.0 * eta * eta + ( 16.0 / 3.0 ) * eta + 4.0 / 3.0;
        real db4 = ( 8.0 / 3.0 ) * eta * eta * eta + 2.0 * eta * eta - ( 1.0 / 3.0 ) * eta - 1.0 / 6.0;

        real dda0 = 8.0 * xi * xi - 4.0 * xi - 1.0 / 3.0;
        real dda1 = -32.0 * xi * xi + 24.0 * xi;
        real dda2 = 48.0 * xi * xi - 10.0;
        real dda3 = -32.0 * xi * xi - 8.0 * xi + 16.0 / 3.0;
        real dda4 = 8.0 * xi * xi + 4.0 * xi - 1.0 / 3.0;

        real ddb0 = 8.0 * eta * eta - 4.0 * eta - 1.0 / 3.0;
        real ddb1 = -32.0 * eta * eta + 24.0 * eta;
        real ddb2 = 48.0 * eta * eta - 10.0;
        real ddb3 = -32.0 * eta * eta - 8.0 * eta + 16.0 / 3.0;
        real ddb4 = 8.0 * eta * eta + 4.0 * eta - 1.0 / 3.0;

        ad2NdXi2.set_size( 3, 25 );
        ad2NdXi2( 0, 0 ) = dda0 * b0;
        ad2NdXi2( 1, 0 ) = a0 * ddb0;
        ad2NdXi2( 2, 0 ) = da0 * db0;
        ad2NdXi2( 0, 1 ) = dda1 * b0;
        ad2NdXi2( 1, 1 ) = a1 * ddb0;
        ad2NdXi2( 2, 1 ) = da1 * db0;
        ad2NdXi2( 0, 2 ) = dda2 * b0;
        ad2NdXi2( 1, 2 ) = a2 * ddb0;
        ad2NdXi2( 2, 2 ) = da2 * db0;
        ad2NdXi2( 0, 3 ) = dda3 * b0;
        ad2NdXi2( 1, 3 ) = a3 * ddb0;
        ad2NdXi2( 2, 3 ) = da3 * db0;
        ad2NdXi2( 0, 4 ) = dda4 * b0;
        ad2NdXi2( 1, 4 ) = a4 * ddb0;
        ad2NdXi2( 2, 4 ) = da4 * db0;

        ad2NdXi2( 0, 5 ) = dda0 * b1;
        ad2NdXi2( 1, 5 ) = a0 * ddb1;
        ad2NdXi2( 2, 5 ) = da0 * db1;
        ad2NdXi2( 0, 6 ) = dda1 * b1;
        ad2NdXi2( 1, 6 ) = a1 * ddb1;
        ad2NdXi2( 2, 6 ) = da1 * db1;
        ad2NdXi2( 0, 7 ) = dda2 * b1;
        ad2NdXi2( 1, 7 ) = a2 * ddb1;
        ad2NdXi2( 2, 7 ) = da2 * db1;
        ad2NdXi2( 0, 8 ) = dda3 * b1;
        ad2NdXi2( 1, 8 ) = a3 * ddb1;
        ad2NdXi2( 2, 8 ) = da3 * db1;
        ad2NdXi2( 0, 9 ) = dda4 * b1;
        ad2NdXi2( 1, 9 ) = a4 * ddb1;
        ad2NdXi2( 2, 9 ) = da4 * db1;

        ad2NdXi2( 0, 10 ) = dda0 * b2;
        ad2NdXi2( 1, 10 ) = a0 * ddb2;
        ad2NdXi2( 2, 10 ) = da0 * db2;
        ad2NdXi2( 0, 11 ) = dda1 * b2;
        ad2NdXi2( 1, 11 ) = a1 * ddb2;
        ad2NdXi2( 2, 11 ) = da1 * db2;
        ad2NdXi2( 0, 12 ) = dda2 * b2;
        ad2NdXi2( 1, 12 ) = a2 * ddb2;
        ad2NdXi2( 2, 12 ) = da2 * db2;
        ad2NdXi2( 0, 13 ) = dda3 * b2;
        ad2NdXi2( 1, 13 ) = a3 * ddb2;
        ad2NdXi2( 2, 13 ) = da3 * db2;
        ad2NdXi2( 0, 14 ) = dda4 * b2;
        ad2NdXi2( 1, 14 ) = a4 * ddb2;
        ad2NdXi2( 2, 14 ) = da4 * db2;

        ad2NdXi2( 0, 15 ) = dda0 * b3;
        ad2NdXi2( 1, 15 ) = a0 * ddb3;
        ad2NdXi2( 2, 15 ) = da0 * db3;
        ad2NdXi2( 0, 16 ) = dda1 * b3;
        ad2NdXi2( 1, 16 ) = a1 * ddb3;
        ad2NdXi2( 2, 16 ) = da1 * db3;
        ad2NdXi2( 0, 17 ) = dda2 * b3;
        ad2NdXi2( 1, 17 ) = a2 * ddb3;
        ad2NdXi2( 2, 17 ) = da2 * db3;
        ad2NdXi2( 0, 18 ) = dda3 * b3;
        ad2NdXi2( 1, 18 ) = a3 * ddb3;
        ad2NdXi2( 2, 18 ) = da3 * db3;
        ad2NdXi2( 0, 19 ) = dda4 * b3;
        ad2NdXi2( 1, 19 ) = a4 * ddb3;
        ad2NdXi2( 2, 19 ) = da4 * db3;

        ad2NdXi2( 0, 20 ) = dda0 * b4;
        ad2NdXi2( 1, 20 ) = a0 * ddb4;
        ad2NdXi2( 2, 20 ) = da0 * db4;
        ad2NdXi2( 0, 21 ) = dda1 * b4;
        ad2NdXi2( 1, 21 ) = a1 * ddb4;
        ad2NdXi2( 2, 21 ) = da1 * db4;
        ad2NdXi2( 0, 22 ) = dda2 * b4;
        ad2NdXi2( 1, 22 ) = a2 * ddb4;
        ad2NdXi2( 2, 22 ) = da2 * db4;
        ad2NdXi2( 0, 23 ) = dda3 * b4;
        ad2NdXi2( 1, 23 ) = a3 * ddb4;
        ad2NdXi2( 2, 23 ) = da3 * db4;
        ad2NdXi2( 0, 24 ) = dda4 * b4;
        ad2NdXi2( 1, 24 ) = a4 * ddb4;
        ad2NdXi2( 2, 24 ) = da4 * db4;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 25 >::eval_d3NdXi3( const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >                                                                                                 &ad3NdXi3 ) const
    {
        MORIS_ASSERT( aXi.length() >= 2, "QUAD25 - eval_d3NdXi3: aXi not allocated or has wrong size." );

        real xi  = aXi( 0 );
        real eta = aXi( 1 );

        real a0 = ( 2.0 / 3.0 ) * ( xi + 0.5 ) * xi * ( xi - 0.5 ) * ( xi - 1.0 );
        real a1 = -( 8.0 / 3.0 ) * ( xi + 1.0 ) * xi * ( xi - 0.5 ) * ( xi - 1.0 );
        real a2 = 4.0 * ( xi + 1.0 ) * ( xi + 0.5 ) * ( xi - 0.5 ) * ( xi - 1.0 );
        real a3 = -( 8.0 / 3.0 ) * ( xi + 1.0 ) * ( xi + 0.5 ) * xi * ( xi - 1.0 );
        real a4 = ( 2.0 / 3.0 ) * ( xi + 1.0 ) * ( xi + 0.5 ) * xi * ( xi - 0.5 );

        real b0 = ( 2.0 / 3.0 ) * ( eta + 0.5 ) * eta * ( eta - 0.5 ) * ( eta - 1.0 );
        real b1 = -( 8.0 / 3.0 ) * ( eta + 1.0 ) * eta * ( eta - 0.5 ) * ( eta - 1.0 );
        real b2 = 4.0 * ( eta + 1.0 ) * ( eta + 0.5 ) * ( eta - 0.5 ) * ( eta - 1.0 );
        real b3 = -( 8.0 / 3.0 ) * ( eta + 1.0 ) * ( eta + 0.5 ) * eta * ( eta - 1.0 );
        real b4 = ( 2.0 / 3.0 ) * ( eta + 1.0 ) * ( eta + 0.5 ) * eta * ( eta - 0.5 );

        real da0 = ( 8.0 / 3.0 ) * xi * xi * xi - 2.0 * xi * xi - ( 1.0 / 3.0 ) * xi + 1.0 / 6.0;
        real da1 = -( 32.0 / 3.0 ) * xi * xi * xi + 12.0 * xi * xi - 4.0 / 3.0;
        real da2 = 16.0 * xi * xi * xi - 10.0 * xi;
        real da3 = -( 32.0 / 3.0 ) * xi * xi * xi - 4.0 * xi * xi + ( 16.0 / 3.0 ) * xi + 4.0 / 3.0;
        real da4 = ( 8.0 / 3.0 ) * xi * xi * xi + 2.0 * xi * xi - ( 1.0 / 3.0 ) * xi - 1.0 / 6.0;

        real db0 = ( 8.0 / 3.0 ) * eta * eta * eta - 2.0 * eta * eta - ( 1.0 / 3.0 ) * eta + 1.0 / 6.0;
        real db1 = -( 32.0 / 3.0 ) * eta * eta * eta + 12.0 * eta * eta - 4.0 / 3.0;
        real db2 = 16.0 * eta * eta * eta - 10.0 * eta;
        real db3 = -( 32.0 / 3.0 ) * eta * eta * eta - 4.0 * eta * eta + ( 16.0 / 3.0 ) * eta + 4.0 / 3.0;
        real db4 = ( 8.0 / 3.0 ) * eta * eta * eta + 2.0 * eta * eta - ( 1.0 / 3.0 ) * eta - 1.0 / 6.0;

        real dda0 = 8.0 * xi * xi - 4.0 * xi - 1.0 / 3.0;
        real dda1 = -32.0 * xi * xi + 24.0 * xi;
        real dda2 = 48.0 * xi * xi - 10.0;
        real dda3 = -32.0 * xi * xi - 8.0 * xi + 16.0 / 3.0;
        real dda4 = 8.0 * xi * xi + 4.0 * xi - 1.0 / 3.0;

        real ddb0 = 8.0 * eta * eta - 4.0 * eta - 1.0 / 3.0;
        real ddb1 = -32.0 * eta * eta + 24.0 * eta;
        real ddb2 = 48.0 * eta * eta - 10.0;
        real ddb3 = -32.0 * eta * eta - 8.0 * eta + 16.0 / 3.0;
        real ddb4 = 8.0 * eta * eta + 4.0 * eta - 1.0 / 3.0;

        real ddda0 = 16.0 * xi - 4.0;
        real ddda1 = -64.0 * xi + 24.0;
        real ddda2 = 96.0 * xi;
        real ddda3 = -64.0 * xi - 8.0;
        real ddda4 = 16.0 * xi + 4.0;

        real dddb0 = 16.0 * eta - 4.0;
        real dddb1 = -64.0 * eta + 24.0;
        real dddb2 = 96.0 * eta;
        real dddb3 = -64.0 * eta - 8.0;
        real dddb4 = 16.0 * eta + 4.0;

        ad3NdXi3.set_size( 4, 25 );
        ad3NdXi3( 0, 0 ) = ddda0 * b0;
        ad3NdXi3( 1, 0 ) = a0 * dddb0;
        ad3NdXi3( 2, 0 ) = dda0 * db0;
        ad3NdXi3( 3, 0 ) = da0 * ddb0;

        ad3NdXi3( 0, 1 ) = ddda1 * b0;
        ad3NdXi3( 1, 1 ) = a1 * dddb0;
        ad3NdXi3( 2, 1 ) = dda1 * db0;
        ad3NdXi3( 3, 1 ) = da1 * ddb0;

        ad3NdXi3( 0, 2 ) = ddda2 * b0;
        ad3NdXi3( 1, 2 ) = a2 * dddb0;
        ad3NdXi3( 2, 2 ) = dda2 * db0;
        ad3NdXi3( 3, 2 ) = da2 * ddb0;

        ad3NdXi3( 0, 3 ) = ddda3 * b0;
        ad3NdXi3( 1, 3 ) = a3 * dddb0;
        ad3NdXi3( 2, 3 ) = dda3 * db0;
        ad3NdXi3( 3, 3 ) = da3 * ddb0;

        ad3NdXi3( 0, 4 ) = ddda4 * b0;
        ad3NdXi3( 1, 4 ) = a4 * dddb0;
        ad3NdXi3( 2, 4 ) = dda4 * db0;
        ad3NdXi3( 3, 4 ) = da4 * ddb0;

        ad3NdXi3( 0, 5 ) = ddda0 * b1;
        ad3NdXi3( 1, 5 ) = a0 * dddb1;
        ad3NdXi3( 2, 5 ) = dda0 * db1;
        ad3NdXi3( 3, 5 ) = da0 * ddb1;

        ad3NdXi3( 0, 6 ) = ddda1 * b1;
        ad3NdXi3( 1, 6 ) = a1 * dddb1;
        ad3NdXi3( 2, 6 ) = dda1 * db1;
        ad3NdXi3( 3, 6 ) = da1 * ddb1;

        ad3NdXi3( 0, 7 ) = ddda2 * b1;
        ad3NdXi3( 1, 7 ) = a2 * dddb1;
        ad3NdXi3( 2, 7 ) = dda2 * db1;
        ad3NdXi3( 3, 7 ) = da2 * ddb1;

        ad3NdXi3( 0, 8 ) = ddda3 * b1;
        ad3NdXi3( 1, 8 ) = a3 * dddb1;
        ad3NdXi3( 2, 8 ) = dda3 * db1;
        ad3NdXi3( 3, 8 ) = da3 * ddb1;

        ad3NdXi3( 0, 9 ) = ddda4 * b1;
        ad3NdXi3( 1, 9 ) = a4 * dddb1;
        ad3NdXi3( 2, 9 ) = dda4 * db1;
        ad3NdXi3( 3, 9 ) = da4 * ddb1;

        ad3NdXi3( 0, 10 ) = ddda0 * b2;
        ad3NdXi3( 1, 10 ) = a0 * dddb2;
        ad3NdXi3( 2, 10 ) = dda0 * db2;
        ad3NdXi3( 3, 10 ) = da0 * ddb2;

        ad3NdXi3( 0, 11 ) = ddda1 * b2;
        ad3NdXi3( 1, 11 ) = a1 * dddb2;
        ad3NdXi3( 2, 11 ) = dda1 * db2;
        ad3NdXi3( 3, 11 ) = da1 * ddb2;

        ad3NdXi3( 0, 12 ) = ddda2 * b2;
        ad3NdXi3( 1, 12 ) = a2 * dddb2;
        ad3NdXi3( 2, 12 ) = dda2 * db2;
        ad3NdXi3( 3, 12 ) = da2 * ddb2;

        ad3NdXi3( 0, 13 ) = ddda3 * b2;
        ad3NdXi3( 1, 13 ) = a3 * dddb2;
        ad3NdXi3( 2, 13 ) = dda3 * db2;
        ad3NdXi3( 3, 13 ) = da3 * ddb2;

        ad3NdXi3( 0, 14 ) = ddda4 * b2;
        ad3NdXi3( 1, 14 ) = a4 * dddb2;
        ad3NdXi3( 2, 14 ) = dda4 * db2;
        ad3NdXi3( 3, 14 ) = da4 * ddb2;

        ad3NdXi3( 0, 15 ) = ddda0 * b3;
        ad3NdXi3( 1, 15 ) = a0 * dddb3;
        ad3NdXi3( 2, 15 ) = dda0 * db3;
        ad3NdXi3( 3, 15 ) = da0 * ddb3;

        ad3NdXi3( 0, 16 ) = ddda1 * b3;
        ad3NdXi3( 1, 16 ) = a1 * dddb3;
        ad3NdXi3( 2, 16 ) = dda1 * db3;
        ad3NdXi3( 3, 16 ) = da1 * ddb3;

        ad3NdXi3( 0, 17 ) = ddda2 * b3;
        ad3NdXi3( 1, 17 ) = a2 * dddb3;
        ad3NdXi3( 2, 17 ) = dda2 * db3;
        ad3NdXi3( 3, 17 ) = da2 * ddb3;

        ad3NdXi3( 0, 18 ) = ddda3 * b3;
        ad3NdXi3( 1, 18 ) = a3 * dddb3;
        ad3NdXi3( 2, 18 ) = dda3 * db3;
        ad3NdXi3( 3, 18 ) = da3 * ddb3;

        ad3NdXi3( 0, 19 ) = ddda4 * b3;
        ad3NdXi3( 1, 19 ) = a4 * dddb3;
        ad3NdXi3( 2, 19 ) = dda4 * db3;
        ad3NdXi3( 3, 19 ) = da4 * ddb3;

        ad3NdXi3( 0, 20 ) = ddda0 * b4;
        ad3NdXi3( 1, 20 ) = a0 * dddb4;
        ad3NdXi3( 2, 20 ) = dda0 * db4;
        ad3NdXi3( 3, 20 ) = da0 * ddb4;

        ad3NdXi3( 0, 21 ) = ddda1 * b4;
        ad3NdXi3( 1, 21 ) = a1 * dddb4;
        ad3NdXi3( 2, 21 ) = dda1 * db4;
        ad3NdXi3( 3, 21 ) = da1 * ddb4;

        ad3NdXi3( 0, 22 ) = ddda2 * b4;
        ad3NdXi3( 1, 22 ) = a2 * dddb4;
        ad3NdXi3( 2, 22 ) = dda2 * db4;
        ad3NdXi3( 3, 22 ) = da2 * ddb4;

        ad3NdXi3( 0, 23 ) = ddda3 * b4;
        ad3NdXi3( 1, 23 ) = a3 * dddb4;
        ad3NdXi3( 2, 23 ) = dda3 * db4;
        ad3NdXi3( 3, 23 ) = da3 * ddb4;

        ad3NdXi3( 0, 24 ) = ddda4 * b4;
        ad3NdXi3( 1, 24 ) = a4 * dddb4;
        ad3NdXi3( 2, 24 ) = dda4 * db4;
        ad3NdXi3( 3, 24 ) = da4 * ddb4;
    }

    //------------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD25_HPP_ */