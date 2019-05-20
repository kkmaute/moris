/*
 * cl_FEM_Interpolation_Function_Lagrange_Bar1.hpp
 *
 *  Created on: Jul 13, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_

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
        uint
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1 >::get_number_of_param_dimensions() const
        {
            return 1;
        }

//------------------------------------------------------------------------------

    template<>
    mtk::Interpolation_Order
    Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1 >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::CONSTANT;
    }

//------------------------------------------------------------------------------

    template<>
    Matrix< DDRMat >
    Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1 >::get_param_coords() const
    {
        Matrix< DDRMat > tXiHat(1,1);
        tXiHat( 0 ) =  0.000000;
        return tXiHat;
    }

//------------------------------------------------------------------------------

        template<>
        Matrix < DDRMat >
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1  >::eval_N(
                const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                    "LINE1 - eval_N: aXi not allocated or hat wrong size." );

            Matrix< DDRMat > tN(1,1);
            tN( 0 ) = 1.0;
            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix < DDRMat >
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                          "LINE1 - eval_dNdXi: aXi not allocated or hat wrong size." );

            Matrix< DDRMat > tdNdXi(1,1);
            tdNdXi( 0 ) = 0.0;
            return tdNdXi;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1  >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                          "LINE1 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            Matrix< DDRMat > td2NdXi2(1,1);
            td2NdXi2( 0 ) =  0.0;
            return td2NdXi2;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 1  >::eval_d3NdXi3( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                          "LINE1 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            Matrix< DDRMat > td3NdXi3(1,1);
            td3NdXi3( 0 ) =  0.0;
            return td3NdXi3;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

//------------------------------------------------------------------------------
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_ */
