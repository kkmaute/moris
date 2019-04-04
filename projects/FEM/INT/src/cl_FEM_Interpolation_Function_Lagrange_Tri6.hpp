/*
 * cl_FEM_Interpolation_Function_Tri6.hpp
 *
 *  Created on: Apr 03, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TRI6_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TRI6_HPP_

#include "assert.h"

#include "typedefs.hpp" //MRS/COR/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src

#include "op_times.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 6 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 6 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat( 2, 6 );
            tXiHat( 0, 0 ) = 0.000000; tXiHat( 1, 0 ) = 0.000000;
            tXiHat( 0, 1 ) = 1.000000; tXiHat( 1, 1 ) = 0.000000;
            tXiHat( 0, 2 ) = 0.000000; tXiHat( 1, 2 ) = 1.000000;
            tXiHat( 0, 3 ) = 0.500000; tXiHat( 1, 3 ) = 0.000000;
            tXiHat( 0, 4 ) = 0.500000; tXiHat( 1, 4 ) = 0.500000;
            tXiHat( 0, 5 ) = 0.000000; tXiHat( 1, 5 ) = 0.500000;
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 6 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "eval_shape: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            // get the triangular coordinates
            real zeta1 = 1 - xi - eta;
            real zeta2 = xi;
            real zeta3 = eta;

            // populate matrix with values
            Matrix< DDRMat > tN( 1, 6 );
            tN( 0 ) = zeta1 * ( 2.0 * zeta1 - 1.0 );
            tN( 1 ) = zeta2 * ( 2.0 * zeta2 - 1.0 );
            tN( 2 ) = zeta3 * ( 2.0 * zeta3 - 1.0 );
            tN( 3 ) = 4.0 * zeta1 * zeta2;
            tN( 4 ) = 4.0 * zeta2 * zeta3;
            tN( 5 ) = 4.0 * zeta3 * zeta1;
            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 6 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "eval_shape: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            // get the triangular coordinates
            real zeta1 = 1 - xi - eta;
            real zeta2 = xi;
            real zeta3 = eta;

            // populate output matrix
//            Matrix< DDRMat > tdNdXi( 2, 6 );
//            tdNdXi( 0, 0 ) = 1.0 - 4.0 * zeta1;
//            tdNdXi( 1, 0 ) = 1.0 - 4.0 * zeta1;
//
//            tdNdXi( 0, 1 ) = 4.0 * zeta2 - 1.0;
//            tdNdXi( 1, 1 ) = 0.0;
//
//            tdNdXi( 0, 2 ) = 0.0;
//            tdNdXi( 1, 2 ) = 4.0 * zeta3 - 1.0;
//
//            tdNdXi( 0, 3 ) = 4.0 * ( zeta1 - zeta2 );
//            tdNdXi( 1, 3 ) = -4.0 * zeta2;
//
//            tdNdXi( 0, 4 ) = 4.0 * zeta3;
//            tdNdXi( 1, 4 ) = 4.0 * zeta2;
//
//            tdNdXi( 0, 5 ) = -4.0 * zeta3;
//            tdNdXi( 1, 5 ) = 4.0 * ( zeta1 - zeta3 );
//            print( tdNdXi, "tdNdXi" );

            Matrix< DDRMat > tdZetadXi = {{ -1.0, 1.0, 0.0 },
                                          { -1.0, 0.0, 1.0 }};

            Matrix< DDRMat > tdNdZeta( 3, 6, 0.0 );
            tdNdZeta( 0, 0 ) = 4 * zeta1 - 1.0;
            tdNdZeta( 1, 1 ) = 4 * zeta2 - 1.0;
            tdNdZeta( 2, 2 ) = 4 * zeta3 - 1.0;

            tdNdZeta( 0, 3 ) = 4 * zeta2;
            tdNdZeta( 1, 3 ) = 4 * zeta1;

            tdNdZeta( 1, 4 ) = 4 * zeta3;
            tdNdZeta( 2, 4 ) = 4 * zeta2;

            tdNdZeta( 0, 5 ) = 4 * zeta3;
            tdNdZeta( 2, 5 ) = 4 * zeta1;

            return tdZetadXi * tdNdZeta;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 6 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "eval_shape: aXi not allocated or hat wrong size." );

            // populate output matrix
//            Matrix< DDRMat > td2NdXi2( 3, 6, 0.0);
//            td2NdXi2( 0, 0 ) = 4.0;
//            td2NdXi2( 1, 0 ) = 4.0;
//            td2NdXi2( 2, 0 ) = 4.0;
//
//            td2NdXi2( 0, 1 ) = 4.0;
//            td2NdXi2( 1, 1 ) = 0.0;
//            td2NdXi2( 2, 1 ) = 0.0;
//
//            td2NdXi2( 0, 2 ) = 0.0;
//            td2NdXi2( 1, 2 ) = 4.0;
//            td2NdXi2( 2, 2 ) = 0.0;
//
//            td2NdXi2( 0, 3 ) = -8.0;
//            td2NdXi2( 1, 3 ) =  0.0;
//            td2NdXi2( 2, 3 ) = -4.0;
//
//            td2NdXi2( 0, 4 ) = 0.0;
//            td2NdXi2( 1, 4 ) = 0.0;
//            td2NdXi2( 2, 4 ) = 4.0;
//
//            td2NdXi2( 0, 5 ) =  0.0;
//            td2NdXi2( 1, 5 ) = -8.0;
//            td2NdXi2( 2, 5 ) = -4.0;
//            print( td2NdXi2, "td2NdXi2" );

            //test
            Matrix< DDRMat > tdNdZeta2( 6, 6, 0.0 );
            tdNdZeta2( 0, 0 ) = 4.0;
            tdNdZeta2( 1, 1 ) = 4.0;
            tdNdZeta2( 2, 2 ) = 4.0;
            tdNdZeta2( 5, 3 ) = 4.0;
            tdNdZeta2( 3, 4 ) = 4.0;
            tdNdZeta2( 4, 5 ) = 4.0;

            Matrix< DDRMat > tdZetadXi2 = {{ 1.0, 1.0, 0.0, 0.0,  0.0, -2.0 },
                                           { 1.0, 0.0, 1.0, 0.0, -2.0,  0.0 },
                                           { 1.0, 0.0, 0.0, 1.0, -1.0, -1.0 }};

            return tdZetadXi2 * tdNdZeta2;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TRI6_HPP_ */
