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
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat( 3, 8 );

            tXiHat( 0, 0 ) = -1.000000;
            tXiHat( 1, 0 ) = -1.000000;
            tXiHat( 2, 0 ) = -1.000000;
            tXiHat( 0, 1 ) =  1.000000;
            tXiHat( 1, 1 ) = -1.000000;
            tXiHat( 2, 1 ) = -1.000000;
            tXiHat( 0, 2 ) =  1.000000;
            tXiHat( 1, 2 ) =  1.000000;
            tXiHat( 2, 2 ) = -1.000000;
            tXiHat( 0, 3 ) = -1.000000;
            tXiHat( 1, 3 ) =  1.000000;
            tXiHat( 2, 3 ) = -1.000000;
            tXiHat( 0, 4 ) = -1.000000;
            tXiHat( 1, 4 ) = -1.000000;
            tXiHat( 2, 4 ) =  1.000000;
            tXiHat( 0, 5 ) =  1.000000;
            tXiHat( 1, 5 ) = -1.000000;
            tXiHat( 2, 5 ) =  1.000000;
            tXiHat( 0, 6 ) =  1.000000;
            tXiHat( 1, 6 ) =  1.000000;
            tXiHat( 2, 6 ) =  1.000000;
            tXiHat( 0, 7 ) = -1.000000;
            tXiHat( 1, 7 ) =  1.000000;
            tXiHat( 2, 7 ) =  1.000000;
            return tXiHat;
       }

//------------------------------------------------------------------------------

         template<>
         Matrix< DDRMat >
         Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_N( const Matrix< DDRMat > & aXi ) const
         {
             // make sure that input is correct
             MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_N: aXi not allocated or hat wrong size." );

             // unpack xi and eta from input vector
             real    xi = aXi( 0 );
             real   eta = aXi( 1 );
             real  zeta = aXi( 2 );

             // populate output matrix
             Matrix< DDRMat > tN(1,8);
             tN( 0 ) =  - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
             tN( 1 ) =    ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
             tN( 2 ) =  - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
             tN( 3 ) =    ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
             tN( 4 ) =    ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
             tN( 5 ) =  - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
             tN( 6 ) =    ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
             tN( 7 ) =  - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
             return tN;

         }

//------------------------------------------------------------------------------

         template<>
         Matrix< DDRMat >
         Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
         {
             // make sure that input is correct
             MORIS_ASSERT( aXi.length() >= 3,
                           "HEX8 - eval_dNdXi: aXi not allocated or hat wrong size." );

             // unpack xi and eta from input vector
             auto   xi = aXi( 0 );
             auto  eta = aXi( 1 );
             auto zeta = aXi( 2 );

             // populate output matrix
             Matrix< DDRMat > tdNdXi(3,8);
             tdNdXi( 0, 0 ) = -(  eta - 1 ) * ( zeta - 1 ) * 0.125;
             tdNdXi( 1, 0 ) = -(   xi - 1 ) * ( zeta - 1 ) * 0.125;
             tdNdXi( 2, 0 ) = -(  eta - 1 ) * (   xi - 1 ) * 0.125;

             tdNdXi( 0, 1 ) =  (  eta - 1 ) * ( zeta - 1 ) * 0.125;
             tdNdXi( 1, 1 ) =  (   xi + 1 ) * ( zeta - 1 ) * 0.125;
             tdNdXi( 2, 1 ) =  (  eta - 1 ) * (   xi + 1 ) * 0.125;

             tdNdXi( 0, 2 ) = -(  eta + 1 ) * ( zeta - 1 ) * 0.125;
             tdNdXi( 1, 2 ) = -(   xi + 1 ) * ( zeta - 1 ) * 0.125;
             tdNdXi( 2, 2 ) = -(  eta + 1 ) * (   xi + 1 ) * 0.125;

             tdNdXi( 0, 3 ) =  (  eta + 1 ) * ( zeta - 1 ) * 0.125;
             tdNdXi( 1, 3 ) =  (   xi - 1 ) * ( zeta - 1 ) * 0.125;
             tdNdXi( 2, 3 ) =  (  eta + 1 ) * (   xi - 1 ) * 0.125;

             tdNdXi( 0, 4 ) =  (  eta - 1 ) * ( zeta + 1 ) * 0.125;
             tdNdXi( 1, 4 ) =  (   xi - 1 ) * ( zeta + 1 ) * 0.125;
             tdNdXi( 2, 4 ) =  (  eta - 1 ) * (   xi - 1 ) * 0.125;

             tdNdXi( 0, 5 ) = -(  eta - 1 ) * ( zeta + 1 ) * 0.125;
             tdNdXi( 1, 5 ) = -(   xi + 1 ) * ( zeta + 1 ) * 0.125;
             tdNdXi( 2, 5 ) = -(  eta - 1 ) * (   xi + 1 ) * 0.125;

             tdNdXi( 0, 6 ) =  (  eta + 1 ) * ( zeta + 1 ) * 0.125;
             tdNdXi( 1, 6 ) =  (   xi + 1 ) * ( zeta + 1 ) * 0.125;
             tdNdXi( 2, 6 ) =  (  eta + 1 ) * (   xi + 1 ) * 0.125;

             tdNdXi( 0, 7 ) = -(  eta + 1 ) * ( zeta + 1 ) * 0.125;
             tdNdXi( 1, 7 ) = -(   xi - 1 ) * ( zeta + 1 ) * 0.125;
             tdNdXi( 2, 7 ) = -(  eta + 1 ) * (   xi - 1 ) * 0.125;
             return tdNdXi;
         }

//------------------------------------------------------------------------------

         template<>
         Matrix< DDRMat >
         Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
         {
             // make sure that input is correct
             MORIS_ASSERT( aXi.length() >= 3,
                     "HEX8 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

             // unpack xi and eta from input vector
             auto   xi = aXi( 0 );
             auto  eta = aXi( 1 );
             auto zeta = aXi( 2 );

             // populate output matrix
             Matrix< DDRMat > td2NdXi2(6,8);
             td2NdXi2( 0, 0 ) = 0.0;
             td2NdXi2( 1, 0 ) = 0.0;
             td2NdXi2( 2, 0 ) = 0.0;
             td2NdXi2( 3, 0 ) = 0.125 * ( - xi + 1.0 );
             td2NdXi2( 4, 0 ) = 0.125 * ( - eta + 1.0 );
             td2NdXi2( 5, 0 ) = 0.125 * ( - zeta + 1.0 );

             td2NdXi2( 0, 1 ) = 0.0;
             td2NdXi2( 1, 1 ) = 0.0;
             td2NdXi2( 2, 1 ) = 0.0;
             td2NdXi2( 3, 1 ) = 0.125 * ( xi + 1.0 );
             td2NdXi2( 4, 1 ) = 0.125 * ( eta - 1.0 );
             td2NdXi2( 5, 1 ) = 0.125 * ( zeta - 1.0 );

             td2NdXi2( 0, 2 ) = 0.0;
             td2NdXi2( 1, 2 ) = 0.0;
             td2NdXi2( 2, 2 ) = 0.0;
             td2NdXi2( 3, 2 ) = 0.125 * ( - xi - 1.0 );
             td2NdXi2( 4, 2 ) = 0.125 * ( - eta - 1.0 );
             td2NdXi2( 5, 2 ) = 0.125 * ( - zeta + 1.0 );

             td2NdXi2( 0, 3 ) = 0.0;
             td2NdXi2( 1, 3 ) = 0.0;
             td2NdXi2( 2, 3 ) = 0.0;
             td2NdXi2( 3, 3 ) = 0.125 * ( xi - 1.0 );
             td2NdXi2( 4, 3 ) = 0.125 * ( eta + 1.0 );
             td2NdXi2( 5, 3 ) = 0.125 * ( zeta - 1.0 );

             td2NdXi2( 0, 4 ) = 0.0;
             td2NdXi2( 1, 4 ) = 0.0;
             td2NdXi2( 2, 4 ) = 0.0;
             td2NdXi2( 3, 4 ) = 0.125 * ( xi - 1.0 );
             td2NdXi2( 4, 4 ) = 0.125 * ( eta - 1.0 );
             td2NdXi2( 5, 4 ) = 0.125 * ( zeta + 1.0 );

             td2NdXi2( 0, 5 ) = 0.0;
             td2NdXi2( 1, 5 ) = 0.0;
             td2NdXi2( 2, 5 ) = 0.0;
             td2NdXi2( 3, 5 ) = 0.125 * ( - xi - 1.0 );
             td2NdXi2( 4, 5 ) = 0.125 * ( - eta + 1.0 );
             td2NdXi2( 5, 5 ) = 0.125 * ( - zeta - 1.0 );

             td2NdXi2( 0, 6 ) = 0.0;
             td2NdXi2( 1, 6 ) = 0.0;
             td2NdXi2( 2, 6 ) = 0.0;
             td2NdXi2( 3, 6 ) = 0.125 * ( xi + 1.0 );
             td2NdXi2( 4, 6 ) = 0.125 * ( eta + 1.0 );
             td2NdXi2( 5, 6 ) = 0.125 * ( zeta + 1.0 );

             td2NdXi2( 0, 7 ) = 0.0;
             td2NdXi2( 1, 7 ) = 0.0;
             td2NdXi2( 2, 7 ) = 0.0;
             td2NdXi2( 3, 7 ) = 0.125 * ( - xi + 1.0 );
             td2NdXi2( 4, 7 ) = 0.125 * ( - eta - 1.0 );
             td2NdXi2( 5, 7 ) = 0.125 * ( - zeta - 1.0 );
             return td2NdXi2;
         }

//------------------------------------------------------------------------------

         template<>
         Matrix< DDRMat >
         Interpolation_Function< mtk::Geometry_Type::HEX, Interpolation_Type::LAGRANGE, 3, 8 >::eval_d3NdXi3( const Matrix< DDRMat > & aXi ) const
         {
             // make sure that input is correct
             MORIS_ASSERT( aXi.length() >= 3,
                     "HEX8 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

             // populate output matrix
             Matrix< DDRMat > td3NdXi3(10,8,0.0);
             return td3NdXi3;
         }
//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX8_HPP_ */
