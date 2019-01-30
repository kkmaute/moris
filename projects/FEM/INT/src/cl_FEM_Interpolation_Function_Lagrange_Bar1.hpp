/*
 * cl_FEM_Interpolation_Function_Lagrange_Bar1.hpp
 *
 *  Created on: Jul 13, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_

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
    Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 1 >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::CONSTANT;
    }

//------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 1  >::get_param_coords(
            Matrix< DDRMat > & aXihat ) const
    {
        aXihat.set_size( 1, 1 );

        aXihat( 0 ) =  0.000000;
    }


//------------------------------------------------------------------------------
        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 1  >::eval_N(
                      Interpolation_Matrix  & aN,
                const Matrix< DDRMat > 		& aXi
        )  const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                    "eval_N: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( aN.n_cols() == 1,
                    "eval_N: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( aN.n_rows() == 1,
                    "eval_N: aN not allocated or hat wrong size." );

            aN( 0 ) = 1.0;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 1  >::eval_dNdXi(
                      Interpolation_Matrix  & adNdXi,
                const Matrix< DDRMat > & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                    "eval_dNdXi: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( adNdXi.n_cols() == 1,
                    "eval_dNdXi: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( adNdXi.n_rows() == 1,
                    "eval_dNdXi: aN not allocated or hat wrong size." );


            adNdXi( 0 ) = 0.0;

        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 1  >::eval_d2NdXi2(
                      Interpolation_Matrix  & ad2NdXi2,
                const Matrix< DDRMat > & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                    "ad2NdXi2: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( ad2NdXi2.n_cols() == 1,
                    "ad2NdXi2: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( ad2NdXi2.n_rows() == 1,
                    "ad2NdXi2: aN not allocated or hat wrong size." );

            ad2NdXi2( 0 ) =  0.0;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

//------------------------------------------------------------------------------
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_ */
