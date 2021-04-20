/*
 * cl_MTK_Interpolation_Function_Tet20.hpp
 *
 *  Created on: Apr 04, 2019
 *      Author: noel
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET20_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET20_HPP_

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
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::get_number_of_param_dimensions() const
        {
            return 4;
        }

//------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::get_interpolation_order() const
        {
            return Interpolation_Order::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::get_param_coords( Matrix< DDRMat > & aXiHat ) const
        {
            real t13 = 1.0/3.0;
            real t23 = 2.0/3.0;

            aXiHat =
            {
                { 1.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t23, t13, t23, t13, 0.0, 0.0, 0.0, 0.0, t13, t13, 0.0, t13 },
                { 0.0, 1.0, 0.0, 0.0, t13, t23, t23, t13, 0.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t13, t13, t13, 0.0 },
                { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, t13, t23, t13, t23, 0.0, 0.0, 0.0, 0.0, t23, t13, t13, 0.0, t13, t13 },
                { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, t13, t23, t13, t23, t13, t23, 0.0, t13, t13, t13 }
            };
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_N( const Matrix< DDRMat > & aXi,
                                                                                                              Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 4, "TET20 - eval_N: aXi not allocated or hat wrong size." );

            // unpack the tet coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );
            real zeta4 = aXi( 3 );

            // populate matrix with values
            aNXi.set_size( 1, 20 );
            aNXi( 0 ) = 0.5 * zeta1 * ( 3.0 * zeta1 - 1.0 ) * ( 3.0 * zeta1 - 2.0 );
            aNXi( 1 ) = 0.5 * zeta2 * ( 3.0 * zeta2 - 1.0 ) * ( 3.0 * zeta2 - 2.0 );
            aNXi( 2 ) = 0.5 * zeta3 * ( 3.0 * zeta3 - 1.0 ) * ( 3.0 * zeta3 - 2.0 );
            aNXi( 3 ) = 0.5 * zeta4 * ( 3.0 * zeta4 - 1.0 ) * ( 3.0 * zeta4 - 2.0 );

            aNXi( 4 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta1 - 1.0 );
            aNXi( 5 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta2 - 1.0 );

            aNXi( 6 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta2 - 1.0 );
            aNXi( 7 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta3 - 1.0 );

            aNXi( 8 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta1 - 1.0 );
            aNXi( 9 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta3 - 1.0 );

            aNXi( 10 ) = 4.5 * zeta1 * zeta4 * ( 3.0 * zeta1 - 1.0 );
            aNXi( 11 ) = 4.5 * zeta1 * zeta4 * ( 3.0 * zeta4 - 1.0 );

            aNXi( 12 ) = 4.5 * zeta2 * zeta4 * ( 3.0 * zeta2 - 1.0 );
            aNXi( 13 ) = 4.5 * zeta2 * zeta4 * ( 3.0 * zeta4 - 1.0 );

            aNXi( 14 ) = 4.5 * zeta3 * zeta4 * ( 3.0 * zeta3 - 1.0 );
            aNXi( 15 ) = 4.5 * zeta3 * zeta4 * ( 3.0 * zeta4 - 1.0 );

            aNXi( 16 ) = 27.0 * zeta1 * zeta2 * zeta3;
            aNXi( 17 ) = 27.0 * zeta1 * zeta2 * zeta4;
            aNXi( 18 ) = 27.0 * zeta2 * zeta3 * zeta4;
            aNXi( 19 ) = 27.0 * zeta1 * zeta3 * zeta4;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_dNdXi( const Matrix< DDRMat > & aXi,
                                                                                                                  Matrix< DDRMat > & adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 4, "TET20 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack the tet coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );
            real zeta4 = aXi( 3 );

            real zeta12 = std::pow( zeta1, 2 );
            real zeta22 = std::pow( zeta2, 2 );
            real zeta32 = std::pow( zeta3, 2 );
            real zeta42 = std::pow( zeta4, 2 );

            // populate output matrix
            adNdXi.set_size( 4, 20 );
            adNdXi( 0, 0 )  = 0.5 * ( 27.0 * zeta12 - 18.0 * zeta1 + 2.0 );
            adNdXi( 0, 1 )  = 0.0;
            adNdXi( 0, 2 )  = 0.0;
            adNdXi( 0, 3 )  = 0.0;
            adNdXi( 0, 4 )  = 4.5 * zeta2 * ( 6.0 * zeta1 - 1.0 );
            adNdXi( 0, 5 )  = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            adNdXi( 0, 6 )  = 0.0;
            adNdXi( 0, 7 )  = 0.0;
            adNdXi( 0, 8 )  = 4.5 * zeta3 * ( 6.0 * zeta1 - 1.0 );
            adNdXi( 0, 9 )  = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            adNdXi( 0, 10 ) = 4.5 * zeta4 * ( 6.0 * zeta1 - 1.0 );
            adNdXi( 0, 11 ) = 4.5 * zeta4 * ( 3.0 * zeta4 - 1.0 );
            adNdXi( 0, 12 )  = 0.0;
            adNdXi( 0, 13 )  = 0.0;
            adNdXi( 0, 14 )  = 0.0;
            adNdXi( 0, 15 )  = 0.0;
            adNdXi( 0, 16 ) = 27.0 * zeta2 * zeta3;
            adNdXi( 0, 17 ) = 27.0 * zeta2 * zeta4;
            adNdXi( 0, 18 )  = 0.0;
            adNdXi( 0, 19 ) = 27.0 * zeta3 * zeta4;

            adNdXi( 1, 0 )  = 0.0;
            adNdXi( 1, 1 )  = 0.5 * ( 27.0 * zeta22 - 18.0 * zeta2 + 2.0 );
            adNdXi( 1, 2 )  = 0.0;
            adNdXi( 1, 3 )  = 0.0;
            adNdXi( 1, 4 )  = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );
            adNdXi( 1, 5 )  = 4.5 * zeta1 * ( 6.0 * zeta2 - 1.0 );
            adNdXi( 1, 6 )  = 4.5 * zeta3 * ( 6.0 * zeta2 - 1.0 );
            adNdXi( 1, 7 )  = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            adNdXi( 1, 8 )  = 0.0;
            adNdXi( 1, 9 )  = 0.0;
            adNdXi( 1, 10 ) = 0.0;
            adNdXi( 1, 11 ) = 0.0;
            adNdXi( 1, 12 ) = 4.5 * zeta4 * ( 6.0 * zeta2 - 1.0 );
            adNdXi( 1, 13 ) = 4.5 * zeta4 * ( 3.0 * zeta4 - 1.0 );
            adNdXi( 1, 14 ) = 0.0;
            adNdXi( 1, 15 )  = 0.0;
            adNdXi( 1, 16 ) = 27.0 * zeta1 * zeta3;
            adNdXi( 1, 17 ) = 27.0 * zeta1 * zeta4;
            adNdXi( 1, 18 ) = 27.0 * zeta3 * zeta4;
            adNdXi( 1, 19 ) = 0.0;

            adNdXi( 2, 0 )  = 0.0;
            adNdXi( 2, 1 )  = 0.0;
            adNdXi( 2, 2 )  = 0.5 * ( 27.0 * zeta32 - 18.0 * zeta3 + 2.0 );
            adNdXi( 2, 3 )  = 0.0;
            adNdXi( 2, 4 )  = 0.0;
            adNdXi( 2, 5 )  = 0.0;
            adNdXi( 2, 6 )  = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            adNdXi( 2, 7 )  = 4.5 * zeta2 * ( 6.0 * zeta3 - 1.0 );
            adNdXi( 2, 8 )  = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );
            adNdXi( 2, 9 )  = 4.5 * zeta1 * ( 6.0 * zeta3 - 1.0 );
            adNdXi( 2, 10 )  = 0.0;
            adNdXi( 2, 11 )  = 0.0;
            adNdXi( 2, 12 )  = 0.0;
            adNdXi( 2, 13 )  = 0.0;
            adNdXi( 2, 14 ) = 4.5 * zeta4 * ( 6.0 * zeta3 - 1.0 );
            adNdXi( 2, 15 ) = 4.5 * zeta4 * ( 3.0 * zeta4 - 1.0 );
            adNdXi( 2, 16 ) = 27.0 * zeta1 * zeta2;
            adNdXi( 2, 17 )  = 0.0;
            adNdXi( 2, 18 ) = 27.0 * zeta2 * zeta4;
            adNdXi( 2, 19 ) = 27.0 * zeta1 * zeta4;

            adNdXi( 3, 0 )  = 0.0;
            adNdXi( 3, 1 )  = 0.0;
            adNdXi( 3, 2 )  = 0.0;
            adNdXi( 3, 3 )  = 0.5 * ( 27.0 * zeta42 - 18.0 * zeta4 + 2.0 );
            adNdXi( 3, 4 )  = 0.0;
            adNdXi( 3, 5 )  = 0.0;
            adNdXi( 3, 6 )  = 0.0;
            adNdXi( 3, 7 )  = 0.0;
            adNdXi( 3, 8 )  = 0.0;
            adNdXi( 3, 9 )  = 0.0;
            adNdXi( 3, 10 ) = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );
            adNdXi( 3, 11 ) = 4.5 * zeta1 * ( 6.0 * zeta4 - 1.0 );
            adNdXi( 3, 12 ) = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            adNdXi( 3, 13 ) = 4.5 * zeta2 * ( 6.0 * zeta4 - 1.0 );
            adNdXi( 3, 14 ) = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            adNdXi( 3, 15 ) = 4.5 * zeta3 * ( 6.0 * zeta4 - 1.0 );
            adNdXi( 3, 16 ) = 0.0;
            adNdXi( 3, 17 ) = 27.0 * zeta1 * zeta2;
            adNdXi( 3, 18 ) = 27.0 * zeta2 * zeta3;
            adNdXi( 3, 19 ) = 27.0 * zeta1 * zeta3;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi,
                                                                                                                    Matrix< DDRMat > & ad2NdXi2 ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 4, "TET20 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack the tet coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );
            real zeta4 = aXi( 3 );

            // populate output matrix
            ad2NdXi2.set_size( 10, 20, 0.0 );
            ad2NdXi2( 0, 0 )  = 9.0 * ( 3.0 * zeta1 - 1.0 );
            ad2NdXi2( 0, 4 )  = 27.0 * zeta2 ;
            ad2NdXi2( 0, 8 )  = 27.0 * zeta3;
            ad2NdXi2( 0, 10 ) = 27.0 * zeta4;

            ad2NdXi2( 1, 1 )  = 9.0 * ( 3.0 * zeta2 - 1.0 );
            ad2NdXi2( 1, 5 )  = 27.0 * zeta1;
            ad2NdXi2( 1, 6 )  = 27.0 * zeta3;
            ad2NdXi2( 1, 12 ) = 27.0 * zeta4 ;

            ad2NdXi2( 2, 2 )  = 9.0 * ( 3.0 * zeta3 - 1.0 );
            ad2NdXi2( 2, 7 )  = 27.0 * zeta2;
            ad2NdXi2( 2, 9 )  = 27.0 * zeta1;
            ad2NdXi2( 2, 14 ) = 27.0 * zeta4;

            ad2NdXi2( 3, 3 )  = 9.0 * ( 3.0 * zeta4 - 1.0 );
            ad2NdXi2( 3, 11 ) = 27.0 * zeta1;
            ad2NdXi2( 3, 13 ) = 27.0 * zeta2;
            ad2NdXi2( 3, 15 ) = 27.0 * zeta3;

            ad2NdXi2( 4, 14 ) = 4.5 * ( 6.0 * zeta3 - 1.0 );
            ad2NdXi2( 4, 15 ) = 4.5 * ( 6.0 * zeta4 - 1.0 );
            ad2NdXi2( 4, 18 ) = 27.0 * zeta2;
            ad2NdXi2( 4, 19 ) = 27.0 * zeta1;

            ad2NdXi2( 5, 12 ) = 4.5 * ( 6.0 * zeta2 - 1.0 );
            ad2NdXi2( 5, 13 ) = 4.5 * ( 6.0 * zeta4 - 1.0 );
            ad2NdXi2( 5, 17 ) = 27.0 * zeta1;
            ad2NdXi2( 5, 18 ) = 27.0 * zeta3;

            ad2NdXi2( 6, 10 ) = 4.5 * ( 6.0 * zeta1 - 1.0 );
            ad2NdXi2( 6, 11 ) = 4.5 * ( 6.0 * zeta4 - 1.0 );
            ad2NdXi2( 6, 17 ) = 27.0 * zeta2;
            ad2NdXi2( 6, 19 ) = 27.0 * zeta3;

            ad2NdXi2( 7, 6 )  = 4.5 * ( 6.0 * zeta2 - 1.0 );
            ad2NdXi2( 7, 7 )  = 4.5 * ( 6.0 * zeta3 - 1.0 );
            ad2NdXi2( 7, 16 ) = 27.0 * zeta1;
            ad2NdXi2( 7, 18 ) = 27.0 * zeta4;

            ad2NdXi2( 8, 8 )  = 4.5 * ( 6.0 * zeta1 - 1.0 );
            ad2NdXi2( 8, 9 )  = 4.5 * ( 6.0 * zeta3 - 1.0 );
            ad2NdXi2( 8, 16 ) = 27.0 * zeta2;
            ad2NdXi2( 8, 19 ) = 27.0 * zeta4;

            ad2NdXi2( 9, 4 )  = 4.5 * ( 6.0 * zeta1 - 1.0 );
            ad2NdXi2( 9, 5 )  = 4.5 * ( 6.0 * zeta2 - 1.0 );
            ad2NdXi2( 9, 16 ) = 27.0 * zeta3;
            ad2NdXi2( 9, 17 ) = 27.0 * zeta4;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_d3NdXi3( const Matrix< DDRMat > & aXi,
                                                                                                                    Matrix< DDRMat > & ad3NdXi3 ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( false, "TET20 - eval_d3NdXi3: 3rd order derivatives not implemented for this element." );

            ad3NdXi3.set_size( 1, 20, 0.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET20_HPP_ */