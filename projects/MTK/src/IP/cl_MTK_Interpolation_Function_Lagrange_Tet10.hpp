/*
 * cl_MTK_Interpolation_Function_Tet10.hpp
 *
 *  Created on: Apr 04, 2019
 *      Author: noel
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET10_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET10_HPP_

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
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 10 >::get_number_of_param_dimensions() const
        {
            return 4;
        }

//------------------------------------------------------------------------------

        template<>
        Interpolation_Order
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 10 >::get_interpolation_order() const
        {
            return Interpolation_Order::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 10 >::get_param_coords( Matrix< DDRMat > & aXiHat ) const
        {
            aXiHat =
            {
                { 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000,
                  0.000000000000000, 0.500000000000000, 0.500000000000000, 0.000000000000000, 0.000000000000000 },
                { 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000,
                  0.500000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000, 0.000000000000000 },
                { 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000,
                  0.500000000000000, 0.500000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000 },
                { 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000,
                  0.000000000000000, 0.000000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000 }
            };
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 10 >::eval_N( const Matrix< DDRMat > & aXi,
                                                                                                              Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TET10 - eval_N: aXi not allocated or hat wrong size." );

            // unpack zeta1, zeta2, zeta3, zeta4 from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );
            real zeta4 = aXi( 3 );

            // populate matrix with values
            aNXi.set_size( 1, 10 );
            aNXi( 0 ) = zeta1*(2.0*zeta1 - 1.0);
            aNXi( 1 ) = zeta2*(2.0*zeta2 - 1.0);
            aNXi( 2 ) = zeta3*(2.0*zeta3 - 1.0);
            aNXi( 3 ) = zeta4*(2.0*zeta4 - 1.0);
            aNXi( 4 ) = 4.0*zeta1*zeta2;
            aNXi( 5 ) = 4.0*zeta2*zeta3;
            aNXi( 6 ) = 4.0*zeta1*zeta3;
            aNXi( 7 ) = 4.0*zeta1*zeta4;
            aNXi( 8 ) = 4.0*zeta2*zeta4;
            aNXi( 9 ) = 4.0*zeta3*zeta4;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 10 >::eval_dNdXi( const Matrix< DDRMat > & aXi,
                                                                                                                  Matrix< DDRMat > & adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 4, "TET10 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack zeta1, zeta2, zeta3, zeta4 from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );
            real zeta4 = aXi( 3 );

            // populate output matrix
            adNdXi.set_size( 4, 10, 0.0 );
            adNdXi( 0, 0 ) = 4.0*zeta1 - 1.0;
            adNdXi( 0, 4 ) = 4.0*zeta2;
            adNdXi( 0, 6 ) = 4.0*zeta3;
            adNdXi( 0, 7 ) = 4.0*zeta4;

            adNdXi( 1, 1 ) = 4.0*zeta2 - 1.0;
            adNdXi( 1, 4 ) = 4.0*zeta1;
            adNdXi( 1, 5 ) = 4.0*zeta3;
            adNdXi( 1, 8 ) = 4.0*zeta4;

            adNdXi( 2, 2 ) = 4.0*zeta3 - 1.0;
            adNdXi( 2, 5 ) = 4.0*zeta2;
            adNdXi( 2, 6 ) = 4.0*zeta1;
            adNdXi( 2, 9 ) = 4.0*zeta4;

            adNdXi( 3, 3 ) = 4.0*zeta4 - 1.0;
            adNdXi( 3, 7 ) = 4.0*zeta1;
            adNdXi( 3, 8 ) = 4.0*zeta2;
            adNdXi( 3, 9 ) = 4.0*zeta3;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 10 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi,
                                                                                                                    Matrix< DDRMat > & ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 4, "TET10 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // populate output matrix
            ad2NdXi2.set_size( 10, 10, 0.0 );
            ad2NdXi2( 0, 0 ) = 4.0;
            ad2NdXi2( 1, 1 ) = 4.0;
            ad2NdXi2( 2, 2 ) = 4.0;
            ad2NdXi2( 3, 3 ) = 4.0;
            ad2NdXi2( 4, 9 ) = 4.0;
            ad2NdXi2( 5, 8 ) = 4.0;
            ad2NdXi2( 6, 7 ) = 4.0;
            ad2NdXi2( 7, 5 ) = 4.0;
            ad2NdXi2( 8, 6 ) = 4.0;
            ad2NdXi2( 9, 4 ) = 4.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 10 >::eval_d3NdXi3( const Matrix< DDRMat > & aXi,
                                                                                                                    Matrix< DDRMat > & ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( false, "TET10 - eval_d3NdXi3: 3rd order derivatives not implemented for this element." );

            ad3NdXi3.set_size( 1, 10, 0.0 );
        }
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_LAGRANGE_TET10_HPP_ */