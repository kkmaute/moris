/*
 * cl_MTK_Interpolation_Function_Tri10.hpp
 *
 *  Created on: Apr 04, 2019
 *      Author: noel
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TRI10_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TRI10_HPP_

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
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::get_number_of_param_dimensions() const
        {
            return 3;
        }

//------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::get_interpolation_order() const
        {
            return Interpolation_Order::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::get_param_coords( Matrix< DDRMat > & aXiHat ) const
        {
            real t13 = 1.0/3.0;
            real t23 = 2.0/3.0;

           aXiHat =
            {
                { 1.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t13, t23, t13 },
                { 0.0, 1.0, 0.0, t13, t23, t23, t13, 0.0, 0.0, t13 },
                { 0.0, 0.0, 1.0, 0.0, 0.0, t13, t23, t23, t13, t13 }
            };
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::eval_N( const Matrix< DDRMat > & aXi,
                                                                                                              Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TRI10 - eval_N: aXi not allocated or hat wrong size." );

            // unpack the triangular coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

            // populate matrix with values
            aNXi.set_size( 1, 10 );
            aNXi( 0 ) = 0.5 * ( 3.0 * zeta1 - 1.0 ) * ( 3.0 * zeta1 - 2.0 ) * zeta1;
            aNXi( 1 ) = 0.5 * ( 3.0 * zeta2 - 1.0 ) * ( 3.0 * zeta2 - 2.0 ) * zeta2;
            aNXi( 2 ) = 0.5 * ( 3.0 * zeta3 - 1.0 ) * ( 3.0 * zeta3 - 2.0 ) * zeta3;
            aNXi( 3 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta1 - 1.0 );
            aNXi( 4 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            aNXi( 5 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta2 - 1.0 );
            aNXi( 6 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            aNXi( 7 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            aNXi( 8 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta1 - 1.0 );
            aNXi( 9 ) = 27.0 * zeta1 * zeta2 * zeta3;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::eval_dNdXi( const Matrix< DDRMat > & aXi,
                                                                                                                  Matrix< DDRMat > & adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "TRI10 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack the triangular coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

            real zeta12 = std::pow( zeta1, 2 );
            real zeta22 = std::pow( zeta2, 2 );
            real zeta32 = std::pow( zeta3, 2 );

            // populate output matrix
            adNdXi.set_size( 3, 10 );
            adNdXi( 0, 0 ) = 0.5 * ( 27.0 * zeta12 - 18.0 * zeta1 + 2.0 );
            adNdXi( 0, 1 ) = 0.0;
            adNdXi( 0, 2 ) = 0.0;
            adNdXi( 0, 3 ) = 4.5 * zeta2 * ( 6.0 * zeta1 - 1.0 );
            adNdXi( 0, 4 ) = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            adNdXi( 0, 5 ) = 0.0;
            adNdXi( 0, 6 ) = 0.0;
            adNdXi( 0, 7 ) = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            adNdXi( 0, 8 ) = 4.5 * zeta3 * ( 6.0 * zeta1 - 1.0 );
            adNdXi( 0, 9 ) = 27.0 * zeta2 * zeta3;

            adNdXi( 1, 0 ) = 0.0;
            adNdXi( 1, 1 ) = 0.5 * ( 27.0 * zeta22 - 18.0 * zeta2 + 2.0 );
            adNdXi( 1, 2 ) = 0.0;
            adNdXi( 1, 3 ) = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );
            adNdXi( 1, 4 ) = 4.5 * zeta1 * ( 6.0 * zeta2 - 1.0 );
            adNdXi( 1, 5 ) = 4.5 * zeta3 * ( 6.0 * zeta2 - 1.0 );
            adNdXi( 1, 6 ) = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            adNdXi( 1, 7 ) = 0.0;
            adNdXi( 1, 8 ) = 0.0;
            adNdXi( 1, 9 ) = 27.0 * zeta1 * zeta3;

            adNdXi( 2, 0 ) = 0.0;
            adNdXi( 2, 1 ) = 0.0;
            adNdXi( 2, 2 ) = 0.5 * ( 27.0 * zeta32 - 18.0 * zeta3 + 2.0 );
            adNdXi( 2, 3 ) = 0.0;
            adNdXi( 2, 4 ) = 0.0;
            adNdXi( 2, 5 ) = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            adNdXi( 2, 6 ) = 4.5 * zeta2 * ( 6.0 * zeta3 - 1.0 );
            adNdXi( 2, 7 ) = 4.5 * zeta1 * ( 6.0 * zeta3 - 1.0 );
            adNdXi( 2, 8 ) = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );
            adNdXi( 2, 9 ) = 27.0 * zeta1 * zeta2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi,
                                                                                                                    Matrix< DDRMat > & ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TRI10 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack the triangular coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

            // populate output matrix
            ad2NdXi2.set_size( 6, 10, 0.0 );
            ad2NdXi2( 0, 0 ) = 9.0 * ( 3.0 * zeta1 - 1.0 );
            ad2NdXi2( 1, 1 ) = 9.0 * ( 3.0 * zeta2 - 1.0 );
            ad2NdXi2( 2, 2 ) = 9.0 * ( 3.0 * zeta3 - 1.0 );

            ad2NdXi2( 0, 3 ) = 27.0 * zeta2;
            ad2NdXi2( 5, 3 ) = 4.5 * ( 6.0 * zeta1 - 1 );

            ad2NdXi2( 1, 4 ) = 27.0 * zeta1;
            ad2NdXi2( 5, 4 ) = 4.5 * ( 6.0 * zeta2 - 1 );

            ad2NdXi2( 1, 5 ) = 27.0 * zeta3;
            ad2NdXi2( 3, 5 ) = 4.5 * ( 6.0 * zeta2 - 1 );

            ad2NdXi2( 2, 6 ) = 27.0 * zeta2;
            ad2NdXi2( 3, 6 ) = 4.5 * ( 6.0 * zeta3 - 1 );

            ad2NdXi2( 2, 7 ) = 27.0 * zeta1;
            ad2NdXi2( 4, 7 ) = 4.5 * ( 6.0 * zeta3 - 1 );

            ad2NdXi2( 0, 8 ) = 27.0 * zeta3;
            ad2NdXi2( 4, 8 ) = 4.5 * ( 6.0 * zeta1 - 1 );

            ad2NdXi2( 3, 9 ) = 27.0 * zeta1;
            ad2NdXi2( 4, 9 ) = 27.0 * zeta2;
            ad2NdXi2( 5, 9 ) = 27.0 * zeta3;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::eval_d3NdXi3( const Matrix< DDRMat > & aXi,
                                                                                                                    Matrix< DDRMat > & ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( false, "TRI10 - eval_d3NdXi3: 3rd order derivatives not implemented for this element." );

            ad3NdXi3.set_size( 1, 10, 0.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TRI10_HPP_ */