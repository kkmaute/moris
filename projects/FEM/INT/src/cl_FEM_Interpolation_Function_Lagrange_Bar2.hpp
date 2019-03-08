/*
 * cl_FEM_Interpolation_Function_Lagrange_Bar2.hpp
 *
 *  Created on: Jul 13, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR2_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR2_HPP_

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
    Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 2 >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::LINEAR;
    }

//------------------------------------------------------------------------------

    template<>
    Matrix< DDRMat >
    Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 2 >::get_param_coords() const
    {
        Matrix< DDRMat > tXiHat(1,2);
        tXiHat( 0 ) = -1.000000;
        tXiHat( 1 ) =  1.000000;
        return tXiHat;
    }

//------------------------------------------------------------------------------

    template<>
    Matrix< DDRMat >
    Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 2 >::eval_N( const Matrix< DDRMat > & aXi) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 1,
                      "eval_N: aXi not allocated or hat wrong size." );

        auto xi = aXi( 0 );

        Matrix< DDRMat > tN(1,2);
        tN( 0 ) = 0.5*(1.0-xi);
        tN( 1 ) = 0.5*(1.0+xi);
        return tN;
    }

//------------------------------------------------------------------------------

    template<>
    Matrix< DDRMat >
    Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 2 >::eval_dNdXi( const Matrix< DDRMat > & aXi )  const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 1,
                      "eval_dNdXi: aXi not allocated or hat wrong size." );

        Matrix< DDRMat > tdNdXi(1,2);
        tdNdXi( 0 ) = -0.5;
        tdNdXi( 1 ) =  0.5;
        return tdNdXi;

    }

//------------------------------------------------------------------------------

    template<>
    Matrix< DDRMat >
    Interpolation_Function< Interpolation_Type::LAGRANGE, 1, 2 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 1,
                      "ad2NdXi2: aXi not allocated or hat wrong size." );

        Matrix< DDRMat > td2dXi2(1,2);
        td2dXi2( 0 ) =  0.0;
        td2dXi2( 1 ) =  0.0;
        return td2dXi2;
    }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

//------------------------------------------------------------------------------
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR2_HPP_ */
