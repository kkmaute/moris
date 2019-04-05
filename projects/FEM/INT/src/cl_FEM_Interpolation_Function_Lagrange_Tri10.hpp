/*
 * cl_FEM_Interpolation_Function_Tri10.hpp
 *
 *  Created on: Apr 04, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TRI10_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TRI10_HPP_

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
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 10 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 10 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat( 2, 10 );
            tXiHat( 0, 0 ) = 0.000000000000000; tXiHat( 1, 0 ) = 0.000000000000000;
            tXiHat( 0, 1 ) = 1.000000000000000; tXiHat( 1, 1 ) = 0.000000000000000;
            tXiHat( 0, 2 ) = 0.000000000000000; tXiHat( 1, 2 ) = 1.000000000000000;
            tXiHat( 0, 3 ) = 0.333333333333334; tXiHat( 1, 3 ) = 0.000000000000000;
            tXiHat( 0, 4 ) = 0.666666666666667; tXiHat( 1, 4 ) = 0.000000000000000;
            tXiHat( 0, 5 ) = 0.666666666666667; tXiHat( 1, 5 ) = 0.333333333333334;
            tXiHat( 0, 6 ) = 0.333333333333334; tXiHat( 1, 6 ) = 0.666666666666667;
            tXiHat( 0, 7 ) = 0.000000000000000; tXiHat( 1, 7 ) = 0.333333333333334;
            tXiHat( 0, 8 ) = 0.000000000000000; tXiHat( 1, 8 ) = 0.666666666666667;
            tXiHat( 0, 9 ) = 0.333333333333334; tXiHat( 1, 9 ) = 0.333333333333334;

            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 10 >::eval_N( const Matrix< DDRMat > & aXi ) const
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
            Matrix< DDRMat > tN( 1, 10 );
            tN( 0 ) = 0.5 * ( 3.0 * zeta1 - 1.0 ) * ( 3.0 * zeta1 - 2.0 ) * zeta1;
            tN( 1 ) = 0.5 * ( 3.0 * zeta2 - 1.0 ) * ( 3.0 * zeta2 - 2.0 ) * zeta2;
            tN( 2 ) = 0.5 * ( 3.0 * zeta3 - 1.0 ) * ( 3.0 * zeta3 - 2.0 ) * zeta3;
            tN( 3 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta1 - 1.0 );
            tN( 4 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            tN( 5 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta2 - 1.0 );
            tN( 6 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            tN( 7 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta1 - 1.0 );
            tN( 8 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            tN( 9 ) = 27.0 * zeta1 * zeta2 * zeta3;
            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 10 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
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

            real zeta12 = std::pow( zeta1, 2 );
            real zeta22 = std::pow( zeta2, 2 );
            real zeta32 = std::pow( zeta3, 2 );

            // populate output matrix
            Matrix< DDRMat > tdZetadXi = {{ -1.0, 1.0, 0.0 },
                                          { -1.0, 0.0, 1.0 }};

            Matrix< DDRMat > tdNdZeta( 3, 10, 0.0 );
            tdNdZeta( 0, 0 ) = 0.5 * ( 27.0 * zeta12 - 18.0 * zeta1 + 2.0 );
            tdNdZeta( 1, 1 ) = 0.5 * ( 27.0 * zeta22 - 18.0 * zeta2 + 2.0 );
            tdNdZeta( 2, 2 ) = 0.5 * ( 27.0 * zeta32 - 18.0 * zeta3 + 2.0 );

            tdNdZeta( 0, 3 ) = 4.5 * zeta2 * ( 6.0 * zeta1 - 1.0 );
            tdNdZeta( 1, 3 ) = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );

            tdNdZeta( 0, 4 ) = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            tdNdZeta( 1, 4 ) = 4.5 * zeta1 * ( 6.0 * zeta2 - 1.0 );

            tdNdZeta( 1, 5 ) = 4.5 * zeta3 * ( 6.0 * zeta2 - 1.0 );
            tdNdZeta( 2, 5 ) = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );

            tdNdZeta( 1, 6 ) = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            tdNdZeta( 2, 6 ) = 4.5 * zeta2 * ( 6.0 * zeta3 - 1.0 );

            tdNdZeta( 0, 7 ) = 4.5 * zeta3 * ( 6.0 * zeta1 - 1.0 );
            tdNdZeta( 2, 7 ) = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );

            tdNdZeta( 0, 8 ) = 4.5 * zeta3 * ( 3.0 * zeta2 - 1.0 );
            tdNdZeta( 2, 8 ) = 4.5 * zeta1 * ( 6.0 * zeta2 - 1.0 );

            tdNdZeta( 0, 9 ) = 27.0 * zeta2 * zeta3;
            tdNdZeta( 1, 9 ) = 27.0 * zeta1 * zeta3;
            tdNdZeta( 2, 9 ) = 27.0 * zeta1 * zeta2;
            return tdZetadXi * tdNdZeta;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 10 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
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
            Matrix< DDRMat > tdNdZeta2( 6, 10, 0.0 );
            tdNdZeta2( 0, 0 ) = 9.0 * ( 3.0 * zeta1 - 1.0 );
            tdNdZeta2( 1, 1 ) = 9.0 * ( 3.0 * zeta2 - 1.0 );
            tdNdZeta2( 2, 2 ) = 9.0 * ( 3.0 * zeta3 - 1.0 );

            tdNdZeta2( 0, 3 ) = 27.0 * zeta2;
            tdNdZeta2( 5, 3 ) = 4.5 * ( 6 * zeta1 - 1 );

            tdNdZeta2( 1, 4 ) = 27.0 * zeta1;
            tdNdZeta2( 5, 4 ) = 4.5 * ( 6 * zeta2 - 1 );

            tdNdZeta2( 1, 5 ) = 27.0 * zeta3;
            tdNdZeta2( 3, 5 ) = 4.5 * ( 6 * zeta2 - 1 );

            tdNdZeta2( 2, 6 ) = 27.0 * zeta2;
            tdNdZeta2( 3, 6 ) = 4.5 * ( 6 * zeta3 - 1 );

            tdNdZeta2( 0, 7 ) = 27.0 * zeta3;
            tdNdZeta2( 4, 7 ) = 4.5 * ( 6 * zeta1 - 1 );

            tdNdZeta2( 2, 8 ) = 27.0 * zeta1;
            tdNdZeta2( 4, 8 ) = 4.5 * ( 6 * zeta3 - 1 );

            tdNdZeta2( 3, 9 ) = 27.0 * zeta1;
            tdNdZeta2( 4, 9 ) = 27.0 * zeta2;
            tdNdZeta2( 4, 9 ) = 27.0 * zeta3;

            Matrix< DDRMat > tdZetadXi2 = {{ 1.0, 1.0, 0.0, 0.0,  0.0, -2.0 },
                                           { 1.0, 0.0, 1.0, 0.0, -2.0,  0.0 },
                                           { 1.0, 0.0, 0.0, 1.0, -1.0, -1.0 }};

            return tdZetadXi2 * tdNdZeta2;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TRI10_HPP_ */
