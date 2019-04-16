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
        uint
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::get_number_of_param_dimensions() const
        {
            return 3;
        }

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat =
            {
                { 1.000000000000000, 0.000000000000000, 0.000000000000000,
                  0.500000000000000, 0.000000000000000, 0.500000000000000 },
                { 0.000000000000000, 1.000000000000000, 0.000000000000000,
                  0.500000000000000, 0.500000000000000, 0.000000000000000 },
                { 0.000000000000000, 0.000000000000000, 1.000000000000000,
                  0.000000000000000, 0.500000000000000, 0.500000000000000 }
            };
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TRI6 - eval_N: aXi not allocated or hat wrong size." );

            // unpack  the triangular coordinates input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

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
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "TRI6 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack  the triangular coordinates input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

            // populate output matrix
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

            return tdNdZeta;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 6 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TRI6 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // populate output matrix
            Matrix< DDRMat > td2NdZeta2( 6, 6, 0.0 );
            td2NdZeta2( 0, 0 ) = 4.0;
            td2NdZeta2( 1, 1 ) = 4.0;
            td2NdZeta2( 2, 2 ) = 4.0;
            td2NdZeta2( 5, 3 ) = 4.0;
            td2NdZeta2( 3, 4 ) = 4.0;
            td2NdZeta2( 4, 5 ) = 4.0;

            return td2NdZeta2;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TRI6_HPP_ */
