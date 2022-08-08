/*
 * cl_MTK_Interpolation_Function_Lagrange_Quad9.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD9_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD9_HPP_

#include "assert.h"

//#include "cl_MTK_Interpolation_Matrix.hpp"
#include "typedefs.hpp"                         //MRS/COR/src
#include "cl_MTK_Enums.hpp"                     //MTK/src
#include "cl_MTK_Interpolation_Function.hpp"    //MTK/src

// does not exist in femdoc

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::get_number_of_param_dimensions() const
        {
            return 2;
        }

        //------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::get_interpolation_order() const
        {
            return Interpolation_Order::QUADRATIC;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::get_param_coords( Matrix< DDRMat >& aXiHat ) const
        {
            aXiHat.set_size( 2, 9, 0.0 );
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
            aXiHat( 0, 8 ) = 0.000000;
            aXiHat( 1, 8 ) = 0.000000;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::eval_N(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "QUAD9 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real xi  = aXi( 0 );
            real eta = aXi( 1 );

            // often used constants
            real c    = xi * eta * 0.25;
            real xi2  = std::pow( xi, 2 );
            real eta2 = std::pow( eta, 2 );

            // populate output matrix
            aNXi.set_size( 1, 9 );
            aNXi( 0 ) = ( c * ( eta - 1.0 ) * ( xi - 1.0 ) );
            aNXi( 1 ) = ( c * ( eta - 1.0 ) * ( xi + 1.0 ) );
            aNXi( 2 ) = ( c * ( eta + 1.0 ) * ( xi + 1.0 ) );
            aNXi( 3 ) = ( c * ( eta + 1.0 ) * ( xi - 1.0 ) );
            aNXi( 4 ) = ( eta * ( 1.0 - xi2 ) * ( eta - 1.0 ) ) * 0.5;
            aNXi( 5 ) = ( xi * ( 1.0 - eta2 ) * ( xi + 1.0 ) ) * 0.5;
            aNXi( 6 ) = ( eta * ( 1.0 - xi2 ) * ( eta + 1.0 ) ) * 0.5;
            aNXi( 7 ) = ( xi * ( 1.0 - eta2 ) * ( xi - 1.0 ) ) * 0.5;
            aNXi( 8 ) = ( eta2 - 1.0 ) * ( xi2 - 1.0 );
        }
        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::eval_dNdXi(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD9 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real xi  = aXi( 0 );
            real eta = aXi( 1 );

            // often used constants
            real c    = xi * eta;
            real xi2  = std::pow( xi, 2 );
            real eta2 = std::pow( eta, 2 );

            // populate output matrix
            adNdXi.set_size( 2, 9 );
            adNdXi( 0, 0 ) = ( eta * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            adNdXi( 1, 0 ) = ( xi * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) ) * 0.25;

            adNdXi( 0, 1 ) = ( eta * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            adNdXi( 1, 1 ) = ( xi * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) ) * 0.25;

            adNdXi( 0, 2 ) = ( eta * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            adNdXi( 1, 2 ) = ( xi * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) ) * 0.25;

            adNdXi( 0, 3 ) = ( eta * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            adNdXi( 1, 3 ) = ( xi * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) ) * 0.25;

            adNdXi( 0, 4 ) = -c * ( eta - 1.0 );
            adNdXi( 1, 4 ) = -( ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) ) * 0.5;

            adNdXi( 0, 5 ) = -( ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.5;
            adNdXi( 1, 5 ) = -c * ( xi + 1.0 );

            adNdXi( 0, 6 ) = -c * ( eta + 1.0 );
            adNdXi( 1, 6 ) = -( ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) ) * 0.5;

            adNdXi( 0, 7 ) = -( ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.5;
            adNdXi( 1, 7 ) = -c * ( xi - 1.0 );

            adNdXi( 0, 8 ) = 2.0 * xi * ( eta2 - 1.0 );
            adNdXi( 1, 8 ) = 2.0 * eta * ( xi2 - 1.0 );
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::eval_d2NdXi2(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD9 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real xi  = aXi( 0 );
            real eta = aXi( 1 );

            // often used constants
            real xi2  = std::pow( xi, 2 );
            real eta2 = std::pow( eta, 2 );

            // populate output matrix
            ad2NdXi2.set_size( 3, 9 );
            ad2NdXi2( 0, 0 ) = ( eta * ( eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 0 ) = ( xi * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 0 ) = ( ( 2.0 * eta - 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 1 ) = ( eta * ( eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 1 ) = ( xi * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 1 ) = ( ( 2.0 * eta - 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 2 ) = ( eta * ( eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 2 ) = ( xi * ( xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 2 ) = ( ( 2.0 * eta + 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 3 ) = ( eta * ( eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 3 ) = ( xi * ( xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 3 ) = ( ( 2.0 * eta + 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 4 ) = -eta * ( eta - 1.0 );
            ad2NdXi2( 1, 4 ) = 1.0 - xi2;
            ad2NdXi2( 2, 4 ) = -xi * ( 2.0 * eta - 1.0 );

            ad2NdXi2( 0, 5 ) = 1.0 - eta2;
            ad2NdXi2( 1, 5 ) = -xi * ( xi + 1.0 );
            ad2NdXi2( 2, 5 ) = -eta * ( 2.0 * xi + 1.0 );

            ad2NdXi2( 0, 6 ) = -eta * ( eta + 1.0 );
            ad2NdXi2( 1, 6 ) = 1.0 - xi2;
            ad2NdXi2( 2, 6 ) = -xi * ( 2.0 * eta + 1.0 );

            ad2NdXi2( 0, 7 ) = 1.0 - eta2;
            ad2NdXi2( 1, 7 ) = -xi * ( xi - 1.0 );
            ad2NdXi2( 2, 7 ) = -eta * ( 2.0 * xi - 1.0 );

            ad2NdXi2( 0, 8 ) = 2.0 * eta2 - 2.0;
            ad2NdXi2( 1, 8 ) = 2.0 * xi2 - 2.0;
            ad2NdXi2( 2, 8 ) = 4.0 * eta * xi;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::eval_d3NdXi3(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD9 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real xi  = aXi( 0 );
            real eta = aXi( 1 );

            // often used parameters
            // 1st dimension
            real da0 = xi - 0.5;
            real da1 = -xi * 2.0;
            real da2 = xi + 0.5;

            real dda0 = 1.0;
            real dda1 = -2.0;
            real dda2 = 1.0;

            // 2nd dimension
            real db0 = eta - 0.5;
            real db1 = -eta * 2.0;
            real db2 = eta + 0.5;

            real ddb0 = 1.0;
            real ddb1 = -2.0;
            real ddb2 = 1.0;

            // 3rd derivatives = 0 for all dimensions
            // row 0: dxi3; row 1: deta3; row 2: dxi2deta; row 3: dxideta2
            ad3NdXi3.set_size( 4, 9, 0.0 );

            // 0th node: (0,0)
            ad3NdXi3( 2, 0 ) = dda0 * db0;
            ad3NdXi3( 3, 0 ) = da0 * ddb0;

            // 1th node: (2,0)
            ad3NdXi3( 2, 1 ) = dda2 * db0;
            ad3NdXi3( 3, 1 ) = da2 * ddb0;

            // 2th node: (2,2)
            ad3NdXi3( 2, 2 ) = dda2 * db2;
            ad3NdXi3( 3, 2 ) = da2 * ddb2;

            // 3th node: (0,2)
            ad3NdXi3( 2, 3 ) = dda0 * db2;
            ad3NdXi3( 3, 3 ) = da0 * ddb2;

            // 4th node: (1,0)
            ad3NdXi3( 2, 4 ) = dda1 * db0;
            ad3NdXi3( 3, 4 ) = da1 * ddb0;

            // 5th node: (2,1)
            ad3NdXi3( 2, 5 ) = dda2 * db1;
            ad3NdXi3( 3, 5 ) = da2 * ddb1;

            // 6th node: (1,2)
            ad3NdXi3( 2, 6 ) = dda1 * db2;
            ad3NdXi3( 3, 6 ) = da1 * ddb2;

            // 7th node: (0,1)
            ad3NdXi3( 2, 7 ) = dda0 * db1;
            ad3NdXi3( 3, 7 ) = da0 * ddb1;

            // 8th node: (1,1)
            ad3NdXi3( 2, 8 ) = dda1 * db1;
            ad3NdXi3( 3, 8 ) = da1 * ddb1;
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_QUAD9_HPP_ */
