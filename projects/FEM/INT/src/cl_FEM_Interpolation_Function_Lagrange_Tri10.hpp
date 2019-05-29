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
        uint
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::get_number_of_param_dimensions() const
        {
            return 3;
        }

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::get_param_coords() const
        {
            real t13 = 1.0/3.0;
            real t23 = 2.0/3.0;

            Matrix< DDRMat > tXiHat =
            {
                { 1.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t13, t23, t13 },
                { 0.0, 1.0, 0.0, t13, t23, t23, t13, 0.0, 0.0, t13 },
                { 0.0, 0.0, 1.0, 0.0, 0.0, t13, t23, t23, t13, t13 }
            };

            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TRI10 - eval_N: aXi not allocated or hat wrong size." );

            // unpack the triangular coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

            // populate matrix with values
            Matrix< DDRMat > tN( 1, 10 );
            tN( 0 ) = 0.5 * ( 3.0 * zeta1 - 1.0 ) * ( 3.0 * zeta1 - 2.0 ) * zeta1;
            tN( 1 ) = 0.5 * ( 3.0 * zeta2 - 1.0 ) * ( 3.0 * zeta2 - 2.0 ) * zeta2;
            tN( 2 ) = 0.5 * ( 3.0 * zeta3 - 1.0 ) * ( 3.0 * zeta3 - 2.0 ) * zeta3;
            tN( 3 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta1 - 1.0 );
            tN( 4 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            tN( 5 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta2 - 1.0 );
            tN( 6 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            tN( 7 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            tN( 8 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta1 - 1.0 );
            tN( 9 ) = 27.0 * zeta1 * zeta2 * zeta3;
            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "TRI10 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack the triangular coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

            real zeta12 = std::pow( zeta1, 2 );
            real zeta22 = std::pow( zeta2, 2 );
            real zeta32 = std::pow( zeta3, 2 );

            // populate output matrix
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

            tdNdZeta( 0, 7 ) = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            tdNdZeta( 2, 7 ) = 4.5 * zeta1 * ( 6.0 * zeta3 - 1.0 );

            tdNdZeta( 0, 8 ) = 4.5 * zeta3 * ( 6.0 * zeta1 - 1.0 );
            tdNdZeta( 2, 8 ) = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );

            tdNdZeta( 0, 9 ) = 27.0 * zeta2 * zeta3;
            tdNdZeta( 1, 9 ) = 27.0 * zeta1 * zeta3;
            tdNdZeta( 2, 9 ) = 27.0 * zeta1 * zeta2;
            return tdNdZeta;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "TRI10 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack the triangular coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );

            // populate output matrix
            Matrix< DDRMat > td2NdZeta2( 6, 10, 0.0 );
            td2NdZeta2( 0, 0 ) = 9.0 * ( 3.0 * zeta1 - 1.0 );
            td2NdZeta2( 1, 1 ) = 9.0 * ( 3.0 * zeta2 - 1.0 );
            td2NdZeta2( 2, 2 ) = 9.0 * ( 3.0 * zeta3 - 1.0 );

            td2NdZeta2( 0, 3 ) = 27.0 * zeta2;
            td2NdZeta2( 5, 3 ) = 4.5 * ( 6.0 * zeta1 - 1 );

            td2NdZeta2( 1, 4 ) = 27.0 * zeta1;
            td2NdZeta2( 5, 4 ) = 4.5 * ( 6.0 * zeta2 - 1 );

            td2NdZeta2( 1, 5 ) = 27.0 * zeta3;
            td2NdZeta2( 3, 5 ) = 4.5 * ( 6.0 * zeta2 - 1 );

            td2NdZeta2( 2, 6 ) = 27.0 * zeta2;
            td2NdZeta2( 3, 6 ) = 4.5 * ( 6.0 * zeta3 - 1 );

            td2NdZeta2( 2, 7 ) = 27.0 * zeta1;
            td2NdZeta2( 4, 7 ) = 4.5 * ( 6.0 * zeta3 - 1 );

            td2NdZeta2( 0, 8 ) = 27.0 * zeta3;
            td2NdZeta2( 4, 8 ) = 4.5 * ( 6.0 * zeta1 - 1 );

            td2NdZeta2( 3, 9 ) = 27.0 * zeta1;
            td2NdZeta2( 4, 9 ) = 27.0 * zeta2;
            td2NdZeta2( 5, 9 ) = 27.0 * zeta3;

            return td2NdZeta2;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TRI, Interpolation_Type::LAGRANGE, 2, 10 >::eval_d3NdXi3( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( false, "TRI10 - eval_d3NdXi3: 3rd order derivatives not implemented for this element." );

            Matrix< DDRMat > td3NdXi3(1,10,0.0);
            return td3NdXi3;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TRI10_HPP_ */
