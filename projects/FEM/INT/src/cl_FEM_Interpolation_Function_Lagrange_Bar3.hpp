/*
 * cl_FEM_Interpolation_Function_Lagrange_Bar3.hpp
 *
 *  Created on: Jul 13, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR3_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR3_HPP_

#include "assert.h"

#include "cl_FEM_Interpolation_Matrix.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src
namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3  >::get_interpolation_order()  const
        {
            return mtk::Interpolation_Order::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3  >::get_param_coords(
                Matrix< DDRMat > & aXihat )  const
        {
            aXihat.set_size( 1, 3 );

            aXihat( 0 ) = -1.000000;
            aXihat( 1 ) =  1.000000;
            aXihat( 2 ) =  0.000000;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3  >::eval_N(
                      Interpolation_Matrix  & aN,
                const Matrix< DDRMat > & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                    "eval_N: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( aN.n_cols() == 2,
                    "eval_N: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( aN.n_rows() == 1,
                    "eval_N: aN not allocated or hat wrong size." );

            auto xi = aXi( 0 );
            auto xi2 = std::pow( xi , 2 );

            aN( 0 ) = 0.5 * ( xi2 - xi );
            aN( 1 ) = 0.5 * ( xi2 + xi );
            aN( 2 ) = 1.0 - xi2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3  >::eval_dNdXi(
                      Interpolation_Matrix  & adNdXi,
                const Matrix< DDRMat > & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                    "eval_dNdXi: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( adNdXi.n_cols() == 3,
                    "eval_dNdXi: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( adNdXi.n_rows() == 1,
                    "eval_dNdXi: aN not allocated or hat wrong size." );

            auto xi = aXi( 0 );

            adNdXi( 0 ) =   xi - 0.5;
            adNdXi( 1 ) =   xi + 0.5;
            adNdXi( 2 ) = - 2.0 * xi;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3  >::eval_d2NdXi2(
                Interpolation_Matrix        & ad2NdXi2,
                const Matrix< DDRMat > & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                    "ad2NdXi2: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( ad2NdXi2.n_cols() == 3,
                    "ad2NdXi2: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( ad2NdXi2.n_rows() == 1,
                    "ad2NdXi2: aN not allocated or hat wrong size." );

            ad2NdXi2( 0 ) =   1.0;
            ad2NdXi2( 1 ) =   1.0;
            ad2NdXi2( 2 ) =  -2.0;
        }
//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR3_HPP_ */
