/*
 * cl_MTK_Interpolation_Function_Tri6.hpp
 *
 *  Created on: Apr 03, 2019
 *      Author: noel
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TRI6_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TRI6_HPP_

#include "assert.h"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_MTK_Interpolation_Function.hpp" //MTK/src
#include "op_times.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::get_number_of_param_dimensions() const
        {
            return 3;
        }

//------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::get_interpolation_order() const
        {
            return Interpolation_Order::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::get_param_coords( Matrix< DDRMat > & aXiHat ) const
        {
            aXiHat =
            {
                { 1.000000000000000, 0.000000000000000, 0.000000000000000,
                  0.500000000000000, 0.000000000000000, 0.500000000000000 },
                { 0.000000000000000, 1.000000000000000, 0.000000000000000,
                  0.500000000000000, 0.500000000000000, 0.000000000000000 },
                { 0.000000000000000, 0.000000000000000, 1.000000000000000,
                  0.000000000000000, 0.500000000000000, 0.500000000000000 }
            };
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::eval_N( const Matrix< DDRMat > & aXi,
                                                                                                             Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TRI6 - eval_N: aXi not allocated or hat wrong size." );

            // unpack  the triangular coordinates input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

            // populate matrix with values
            aNXi.set_size( 1, 6 );
            aNXi( 0 ) = zeta1 * ( 2.0 * zeta1 - 1.0 );
            aNXi( 1 ) = zeta2 * ( 2.0 * zeta2 - 1.0 );
            aNXi( 2 ) = zeta3 * ( 2.0 * zeta3 - 1.0 );
            aNXi( 3 ) = 4.0 * zeta1 * zeta2;
            aNXi( 4 ) = 4.0 * zeta2 * zeta3;
            aNXi( 5 ) = 4.0 * zeta3 * zeta1;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::eval_dNdXi( const Matrix< DDRMat > & aXi,
                                                                                                                 Matrix< DDRMat > & adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "TRI6 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack  the triangular coordinates input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

            // populate output matrix
            adNdXi.set_size( 3, 6 );
            adNdXi( 0, 0 ) = 4 * zeta1 - 1.0;
            adNdXi( 0, 1 ) = 0.0;
            adNdXi( 0, 2 ) = 0.0;
            adNdXi( 0, 3 ) = 4 * zeta2;
            adNdXi( 0, 4 ) = 0.0;
            adNdXi( 0, 5 ) = 4 * zeta3;

            adNdXi( 1, 0 ) = 0.0;
            adNdXi( 1, 1 ) = 4 * zeta2 - 1.0;
            adNdXi( 1, 2 ) = 0.0;
            adNdXi( 1, 3 ) = 4 * zeta1;
            adNdXi( 1, 4 ) = 4 * zeta3;
            adNdXi( 1, 5 ) = 0.0;

            adNdXi( 2, 0 ) = 0.0;
            adNdXi( 2, 1 ) = 0.0;
            adNdXi( 2, 2 ) = 4 * zeta3 - 1.0;
            adNdXi( 2, 3 ) = 0.0;
            adNdXi( 2, 4 ) = 4 * zeta2;
            adNdXi( 2, 5 ) = 4 * zeta1;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi,
                                                                                                                   Matrix< DDRMat > & ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TRI6 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // populate output matrix
            ad2NdXi2.set_size( 6, 6, 0.0 );
            ad2NdXi2( 0, 0 ) = 4.0;
            ad2NdXi2( 1, 1 ) = 4.0;
            ad2NdXi2( 2, 2 ) = 4.0;
            ad2NdXi2( 5, 3 ) = 4.0;
            ad2NdXi2( 3, 4 ) = 4.0;
            ad2NdXi2( 4, 5 ) = 4.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::eval_d3NdXi3( const Matrix< DDRMat > & aXi,
                                                                                                                   Matrix< DDRMat > & ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( false , "TRI6 - eval_d3NdXi3: 3rd order derivatives not implemented for this element." );

            ad3NdXi3.set_size(1,6,0.0);
        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TRI6_HPP_ */