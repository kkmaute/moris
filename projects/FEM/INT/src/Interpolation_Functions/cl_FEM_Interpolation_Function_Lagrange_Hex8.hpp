/*
 * cl_FEM_Interpolation_Function_Lagrange_Hex8.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_

#include "assert.h"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src

// checked against femdoc:
// - N
// - dNdxi
// - d2Ndxi2
// - d3Ndxi3 by FD

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::get_number_of_param_dimensions() const
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::LINEAR;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::get_param_coords( Matrix< DDRMat > & aXiHat ) const
        {
            aXiHat.set_size( 3, 8, 0.0 );

            aXiHat( 0, 0 ) = -1.000000;
            aXiHat( 1, 0 ) = -1.000000;
            aXiHat( 2, 0 ) = -1.000000;
            aXiHat( 0, 1 ) =  1.000000;
            aXiHat( 1, 1 ) = -1.000000;
            aXiHat( 2, 1 ) = -1.000000;
            aXiHat( 0, 2 ) =  1.000000;
            aXiHat( 1, 2 ) =  1.000000;
            aXiHat( 2, 2 ) = -1.000000;
            aXiHat( 0, 3 ) = -1.000000;
            aXiHat( 1, 3 ) =  1.000000;
            aXiHat( 2, 3 ) = -1.000000;
            aXiHat( 0, 4 ) = -1.000000;
            aXiHat( 1, 4 ) = -1.000000;
            aXiHat( 2, 4 ) =  1.000000;
            aXiHat( 0, 5 ) =  1.000000;
            aXiHat( 1, 5 ) = -1.000000;
            aXiHat( 2, 5 ) =  1.000000;
            aXiHat( 0, 6 ) =  1.000000;
            aXiHat( 1, 6 ) =  1.000000;
            aXiHat( 2, 6 ) =  1.000000;
            aXiHat( 0, 7 ) = -1.000000;
            aXiHat( 1, 7 ) =  1.000000;
            aXiHat( 2, 7 ) =  1.000000;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_N(
                const Matrix< DDRMat > & aXi,
                Matrix< DDRMat >       & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real    xi = aXi( 0 );
            real   eta = aXi( 1 );
            real  zeta = aXi( 2 );

            // populate output matrix
            aNXi.set_size(1,8);
            aNXi( 0 ) =  - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aNXi( 1 ) =    ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aNXi( 2 ) =  - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aNXi( 3 ) =    ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
            aNXi( 4 ) =    ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aNXi( 5 ) =  - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aNXi( 6 ) =    ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
            aNXi( 7 ) =  - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_dNdXi(
                const Matrix< DDRMat > & aXi,
                Matrix< DDRMat >       & adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real   xi = aXi( 0 );
            real  eta = aXi( 1 );
            real zeta = aXi( 2 );

            // populate output matrix
            adNdXi.set_size( 3, 8 );
            adNdXi( 0, 0 ) = -(  eta - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 0 ) = -(   xi - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 0 ) = -(  eta - 1 ) * (   xi - 1 ) * 0.125;

            adNdXi( 0, 1 ) =  (  eta - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 1 ) =  (   xi + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 1 ) =  (  eta - 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 2 ) = -(  eta + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 2 ) = -(   xi + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 2 ) = -(  eta + 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 3 ) =  (  eta + 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 1, 3 ) =  (   xi - 1 ) * ( zeta - 1 ) * 0.125;
            adNdXi( 2, 3 ) =  (  eta + 1 ) * (   xi - 1 ) * 0.125;

            adNdXi( 0, 4 ) =  (  eta - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 4 ) =  (   xi - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 4 ) =  (  eta - 1 ) * (   xi - 1 ) * 0.125;

            adNdXi( 0, 5 ) = -(  eta - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 5 ) = -(   xi + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 5 ) = -(  eta - 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 6 ) =  (  eta + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 6 ) =  (   xi + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 6 ) =  (  eta + 1 ) * (   xi + 1 ) * 0.125;

            adNdXi( 0, 7 ) = -(  eta + 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 1, 7 ) = -(   xi - 1 ) * ( zeta + 1 ) * 0.125;
            adNdXi( 2, 7 ) = -(  eta + 1 ) * (   xi - 1 ) * 0.125;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_d2NdXi2(
                const Matrix< DDRMat > & aXi,
                Matrix< DDRMat >       & ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real   xi = aXi( 0 );
            real  eta = aXi( 1 );
            real zeta = aXi( 2 );

            // populate output matrix
            ad2NdXi2.set_size( 6, 8 );
            ad2NdXi2( 0, 0 ) = 0.0;
            ad2NdXi2( 1, 0 ) = 0.0;
            ad2NdXi2( 2, 0 ) = 0.0;
            ad2NdXi2( 3, 0 ) = 0.125 * ( - xi + 1.0 );
            ad2NdXi2( 4, 0 ) = 0.125 * ( - eta + 1.0 );
            ad2NdXi2( 5, 0 ) = 0.125 * ( - zeta + 1.0 );

            ad2NdXi2( 0, 1 ) = 0.0;
            ad2NdXi2( 1, 1 ) = 0.0;
            ad2NdXi2( 2, 1 ) = 0.0;
            ad2NdXi2( 3, 1 ) = 0.125 * ( xi + 1.0 );
            ad2NdXi2( 4, 1 ) = 0.125 * ( eta - 1.0 );
            ad2NdXi2( 5, 1 ) = 0.125 * ( zeta - 1.0 );

            ad2NdXi2( 0, 2 ) = 0.0;
            ad2NdXi2( 1, 2 ) = 0.0;
            ad2NdXi2( 2, 2 ) = 0.0;
            ad2NdXi2( 3, 2 ) = 0.125 * ( - xi - 1.0 );
            ad2NdXi2( 4, 2 ) = 0.125 * ( - eta - 1.0 );
            ad2NdXi2( 5, 2 ) = 0.125 * ( - zeta + 1.0 );

            ad2NdXi2( 0, 3 ) = 0.0;
            ad2NdXi2( 1, 3 ) = 0.0;
            ad2NdXi2( 2, 3 ) = 0.0;
            ad2NdXi2( 3, 3 ) = 0.125 * ( xi - 1.0 );
            ad2NdXi2( 4, 3 ) = 0.125 * ( eta + 1.0 );
            ad2NdXi2( 5, 3 ) = 0.125 * ( zeta - 1.0 );

            ad2NdXi2( 0, 4 ) = 0.0;
            ad2NdXi2( 1, 4 ) = 0.0;
            ad2NdXi2( 2, 4 ) = 0.0;
            ad2NdXi2( 3, 4 ) = 0.125 * ( xi - 1.0 );
            ad2NdXi2( 4, 4 ) = 0.125 * ( eta - 1.0 );
            ad2NdXi2( 5, 4 ) = 0.125 * ( zeta + 1.0 );

            ad2NdXi2( 0, 5 ) = 0.0;
            ad2NdXi2( 1, 5 ) = 0.0;
            ad2NdXi2( 2, 5 ) = 0.0;
            ad2NdXi2( 3, 5 ) = 0.125 * ( - xi - 1.0 );
            ad2NdXi2( 4, 5 ) = 0.125 * ( - eta + 1.0 );
            ad2NdXi2( 5, 5 ) = 0.125 * ( - zeta - 1.0 );

            ad2NdXi2( 0, 6 ) = 0.0;
            ad2NdXi2( 1, 6 ) = 0.0;
            ad2NdXi2( 2, 6 ) = 0.0;
            ad2NdXi2( 3, 6 ) = 0.125 * ( xi + 1.0 );
            ad2NdXi2( 4, 6 ) = 0.125 * ( eta + 1.0 );
            ad2NdXi2( 5, 6 ) = 0.125 * ( zeta + 1.0 );

            ad2NdXi2( 0, 7 ) = 0.0;
            ad2NdXi2( 1, 7 ) = 0.0;
            ad2NdXi2( 2, 7 ) = 0.0;
            ad2NdXi2( 3, 7 ) = 0.125 * ( - xi + 1.0 );
            ad2NdXi2( 4, 7 ) = 0.125 * ( - eta - 1.0 );
            ad2NdXi2( 5, 7 ) = 0.125 * ( - zeta - 1.0 );
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_d3NdXi3(
                const Matrix< DDRMat > & aXi,
                Matrix< DDRMat >       & ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            // populate output matrix
            ad3NdXi3.set_size( 10, 8, 0.0 );
            ad3NdXi3( 9, 0 ) = - 0.125;
            ad3NdXi3( 9, 1 ) =   0.125;
            ad3NdXi3( 9, 2 ) = - 0.125;
            ad3NdXi3( 9, 3 ) =   0.125;
            ad3NdXi3( 9, 4 ) =   0.125;
            ad3NdXi3( 9, 5 ) = - 0.125;
            ad3NdXi3( 9, 6 ) =   0.125;
            ad3NdXi3( 9, 7 ) = - 0.125;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_ */
