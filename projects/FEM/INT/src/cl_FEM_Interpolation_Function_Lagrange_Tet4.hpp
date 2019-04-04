/*
 * cl_FEM_Interpolation_Function_Tet4.hpp
 *
 *  Created on: Apr 04, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TET4_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TET4_HPP_

#include "assert.h"

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
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 4 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 4 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat( 3, 4 );
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 4 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "eval_shape: aXi not allocated or hat wrong size." );

//            // unpack xi and eta from input vector
//            real  xi  = aXi( 0 );
//            real eta  = aXi( 1 );
//            real zeta = aXi( 2 );

            // populate matrix with values
            Matrix< DDRMat > tN( 1, 4 );
            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 4 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "eval_shape: aXi not allocated or hat wrong size." );

//            // unpack xi and eta from input vector
//            real  xi  = aXi( 0 );
//            real eta  = aXi( 1 );
//            real zeta = aXi( 2 );


            // populate output matrix
            Matrix< DDRMat > tdNdXi( 3, 4 );
            return tdNdXi;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 4 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "eval_shape: aXi not allocated or hat wrong size." );

//            // unpack xi and eta from input vector
//            real  xi  = aXi( 0 );
//            real eta  = aXi( 1 );
//            real zeta = aXi( 2 );

            // populate output matrix
            Matrix< DDRMat > td2NdXi2( 6, 4 );
            return td2NdXi2;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TET4_HPP_ */
