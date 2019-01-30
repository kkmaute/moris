/*
 * cl_FEM_Interpolation_Function_Quad4.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD4_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD4_HPP_

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
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 4 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 4  >::get_param_coords(
                Matrix< DDRMat > & aXihat ) const
        {
            aXihat.set_size( 2, 4 );

            aXihat( 0, 0 ) = -1.000000;
            aXihat( 1, 0 ) = -1.000000;
            aXihat( 0, 1 ) =  1.000000;
            aXihat( 1, 1 ) = -1.000000;
            aXihat( 0, 2 ) =  1.000000;
            aXihat( 1, 2 ) =  1.000000;
            aXihat( 0, 3 ) = -1.000000;
            aXihat( 1, 3 ) =  1.000000;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 4  >::eval_N(
                  Interpolation_Matrix  & aN,
            const Matrix< DDRMat > 		& aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( aN.n_cols() == 4,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( aN.n_rows() == 1,
                    "eval_shape: aN not allocated or hat wrong size." );


            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // populate matrix with values
            aN( 0 ) = ( ( 1.0 - xi ) * ( 1.0 - eta ) ) * 0.25;
            aN( 1 ) = ( ( 1.0 + xi ) * ( 1.0 - eta ) ) * 0.25;
            aN( 2 ) = ( ( 1.0 + xi ) * ( 1.0 + eta ) ) * 0.25;
            aN( 3 ) = ( ( 1.0 - xi ) * ( 1.0 + eta ) ) * 0.25;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 4  >::eval_dNdXi(
                       Interpolation_Matrix & adNdXi,
                const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( adNdXi.n_cols() == 4,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( adNdXi.n_rows() == 2,
                    "eval_shape: aN not allocated or hat wrong size." );


            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // populate output matrix

            adNdXi( 0, 0 ) =  0.25 * ( eta - 1.0 );
            adNdXi( 1, 0 ) =  0.25 * ( xi - 1.0 );

            adNdXi( 0, 1 ) = -0.25 * ( eta - 1.0 );
            adNdXi( 1, 1 ) = -0.25 * ( xi + 1.0 );

            adNdXi( 0, 2 ) =  0.25 * ( eta + 1.0 );
            adNdXi( 1, 2 ) =  0.25 * ( xi + 1.0 );

            adNdXi( 0, 3 ) = -0.25 * ( eta + 1.0 );
            adNdXi( 1, 3 ) = -0.25 * ( xi - 1.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 4  >::eval_d2NdXi2(
                      Interpolation_Matrix  & ad2NdXi2,
                const Matrix< DDRMat > & aXi ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( ad2NdXi2.n_cols() == 4,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( ad2NdXi2.n_rows() == 3,
                    "eval_shape: aN not allocated or hat wrong size." );


            // populate output matrix

            ad2NdXi2( 0, 0 ) =  0.0;
            ad2NdXi2( 1, 0 ) =  0.0;
            ad2NdXi2( 2, 0 ) =  0.25;

            ad2NdXi2( 0, 1 ) =  0.0;
            ad2NdXi2( 1, 1 ) =  0.0;
            ad2NdXi2( 2, 1 ) = -0.25;

            ad2NdXi2( 0, 2 ) =  0.0;
            ad2NdXi2( 1, 2 ) =  0.0;
            ad2NdXi2( 2, 2 ) =  0.25;

            ad2NdXi2( 0, 3 ) =  0.0;
            ad2NdXi2( 1, 3 ) =  0.0;
            ad2NdXi2( 2, 3 ) = -0.25;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD4_HPP_ */
