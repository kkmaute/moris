/*
 * cl_FEM_Interpolation_Function_Lagrange_Bar3.hpp
 *
 *  Created on: Jul 13, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR3_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR3_HPP_

#include "assert.h"

//#include "cl_FEM_Interpolation_Matrix.hpp"
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
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3 >::get_interpolation_order()  const
        {
            return mtk::Interpolation_Order::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat(1,3);
            tXiHat( 0 ) = -1.000000;
            tXiHat( 1 ) =  1.000000;
            tXiHat( 2 ) =  0.000000;
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                          "eval_N: aXi not allocated or hat wrong size." );

            auto xi = aXi( 0 );
            auto xi2 = std::pow( xi , 2 );

            Matrix< DDRMat > tN(1,3);
            tN( 0 ) = 0.5 * ( xi2 - xi );
            tN( 1 ) = 0.5 * ( xi2 + xi );
            tN( 2 ) = 1.0 - xi2;
            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                          "eval_dNdXi: aXi not allocated or hat wrong size." );

            auto xi = aXi( 0 );
            Matrix< DDRMat > tdNdXi(1,3);
            tdNdXi( 0 ) =   xi - 0.5;
            tdNdXi( 1 ) =   xi + 0.5;
            tdNdXi( 2 ) = - 2.0 * xi;
            return tdNdXi;

        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 3 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                          "ad2NdXi2: aXi not allocated or hat wrong size." );

            Matrix< DDRMat > td2NdXi2(1,3);
            td2NdXi2( 0 ) =   1.0;
            td2NdXi2( 1 ) =   1.0;
            td2NdXi2( 2 ) =  -2.0;
            return td2NdXi2;
        }
//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR3_HPP_ */
