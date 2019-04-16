/*
 * cl_FEM_Interpolation_Function_Tet20.hpp
 *
 *  Created on: Apr 04, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TET20_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TET20_HPP_

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
        Interpolation_Function< mtk::Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::get_number_of_param_dimensions() const
        {
            return 4;
        }

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< mtk::Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::get_param_coords() const
        {
            real t13 = 1.0/3.0;
            real t23 = 2.0/3.0;

            Matrix< DDRMat > tXiHat =
            {
                { 1.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t23, t13, t23, t13, 0.0, 0.0, 0.0, 0.0, t13, t13, 0.0, t13 },
                { 0.0, 1.0, 0.0, 0.0, t13, t23, t23, t13, 0.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t13, t13, t13, 0.0 },
                { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, t13, t23, t13, t23, 0.0, 0.0, 0.0, 0.0, t23, t13, t13, 0.0, t13, t13 },
                { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, t13, t23, t13, t23, t13, t23, 0.0, t13, t13, t13 }
            };
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 4, "TET20 - eval_N: aXi not allocated or hat wrong size." );

            // unpack the tet coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );
            real zeta4 = aXi( 3 );

            // populate matrix with values
            Matrix< DDRMat > tN( 1, 20 );
            tN( 0 ) = 0.5 * zeta1 * ( 3.0 * zeta1 - 1.0 ) * ( 3.0 * zeta1 - 2.0 );
            tN( 1 ) = 0.5 * zeta2 * ( 3.0 * zeta2 - 1.0 ) * ( 3.0 * zeta2 - 2.0 );
            tN( 2 ) = 0.5 * zeta3 * ( 3.0 * zeta3 - 1.0 ) * ( 3.0 * zeta3 - 2.0 );
            tN( 3 ) = 0.5 * zeta4 * ( 3.0 * zeta4 - 1.0 ) * ( 3.0 * zeta4 - 2.0 );

            tN( 4 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta1 - 1.0 );
            tN( 5 ) = 4.5 * zeta1 * zeta2 * ( 3.0 * zeta2 - 1.0 );

            tN( 6 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta2 - 1.0 );
            tN( 7 ) = 4.5 * zeta2 * zeta3 * ( 3.0 * zeta3 - 1.0 );

            tN( 8 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta1 - 1.0 );
            tN( 9 ) = 4.5 * zeta1 * zeta3 * ( 3.0 * zeta3 - 1.0 );

            tN( 10 ) = 4.5 * zeta1 * zeta4 * ( 3.0 * zeta1 - 1.0 );
            tN( 11 ) = 4.5 * zeta1 * zeta4 * ( 3.0 * zeta4 - 1.0 );

            tN( 12 ) = 4.5 * zeta2 * zeta4 * ( 3.0 * zeta2 - 1.0 );
            tN( 13 ) = 4.5 * zeta2 * zeta4 * ( 3.0 * zeta4 - 1.0 );

            tN( 14 ) = 4.5 * zeta3 * zeta4 * ( 3.0 * zeta3 - 1.0 );
            tN( 15 ) = 4.5 * zeta3 * zeta4 * ( 3.0 * zeta4 - 1.0 );

            tN( 16 ) = 27.0 * zeta1 * zeta2 * zeta3;
            tN( 17 ) = 27.0 * zeta1 * zeta2 * zeta4;
            tN( 18 ) = 27.0 * zeta2 * zeta3 * zeta4;
            tN( 19 ) = 27.0 * zeta1 * zeta3 * zeta4;

            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 4, "TET20 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack the tet coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );
            real zeta4 = aXi( 3 );

            real zeta12 = std::pow( zeta1, 2 );
            real zeta22 = std::pow( zeta2, 2 );
            real zeta32 = std::pow( zeta3, 2 );
            real zeta42 = std::pow( zeta4, 2 );

            // populate output matrix
            Matrix< DDRMat > tdNdXi( 4, 20, 0.0 );
            tdNdXi( 0, 0 )  = 0.5 * ( 27.0 * zeta12 - 18.0 * zeta1 + 2.0 );
            tdNdXi( 0, 4 )  = 4.5 * zeta2 * ( 6.0 * zeta1 - 1.0 );
            tdNdXi( 0, 5 )  = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            tdNdXi( 0, 8 )  = 4.5 * zeta3 * ( 6.0 * zeta1 - 1.0 );
            tdNdXi( 0, 9 )  = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            tdNdXi( 0, 10 ) = 4.5 * zeta4 * ( 6.0 * zeta1 - 1.0 );
            tdNdXi( 0, 11 ) = 4.5 * zeta4 * ( 3.0 * zeta4 - 1.0 );
            tdNdXi( 0, 16 ) = 27.0 * zeta2 * zeta3;
            tdNdXi( 0, 17 ) = 27.0 * zeta2 * zeta4;
            tdNdXi( 0, 19 ) = 27.0 * zeta3 * zeta4;

            tdNdXi( 1, 1 )  = 0.5 * ( 27.0 * zeta22 - 18.0 * zeta2 + 2.0 );
            tdNdXi( 1, 4 )  = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );
            tdNdXi( 1, 5 )  = 4.5 * zeta1 * ( 6.0 * zeta2 - 1.0 );
            tdNdXi( 1, 6 )  = 4.5 * zeta3 * ( 6.0 * zeta2 - 1.0 );
            tdNdXi( 1, 7 )  = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            tdNdXi( 1, 12 ) = 4.5 * zeta4 * ( 6.0 * zeta2 - 1.0 );
            tdNdXi( 1, 13 ) = 4.5 * zeta4 * ( 3.0 * zeta4 - 1.0 );
            tdNdXi( 1, 16 ) = 27.0 * zeta1 * zeta3;
            tdNdXi( 1, 17 ) = 27.0 * zeta1 * zeta4;
            tdNdXi( 1, 18 ) = 27.0 * zeta3 * zeta4;

            tdNdXi( 2, 2 )  = 0.5 * ( 27.0 * zeta32 - 18.0 * zeta3 + 2.0 );
            tdNdXi( 2, 6 )  = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            tdNdXi( 2, 7 )  = 4.5 * zeta2 * ( 6.0 * zeta3 - 1.0 );
            tdNdXi( 2, 8 )  = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );
            tdNdXi( 2, 9 )  = 4.5 * zeta1 * ( 6.0 * zeta3 - 1.0 );
            tdNdXi( 2, 14 ) = 4.5 * zeta4 * ( 6.0 * zeta3 - 1.0 );
            tdNdXi( 2, 15 ) = 4.5 * zeta4 * ( 3.0 * zeta4 - 1.0 );
            tdNdXi( 2, 16 ) = 27.0 * zeta1 * zeta2;
            tdNdXi( 2, 18 ) = 27.0 * zeta2 * zeta4;
            tdNdXi( 2, 19 ) = 27.0 * zeta1 * zeta4;

            tdNdXi( 3, 3 )  = 0.5 * ( 27.0 * zeta42 - 18.0 * zeta4 + 2.0 );
            tdNdXi( 3, 10 ) = 4.5 * zeta1 * ( 3.0 * zeta1 - 1.0 );
            tdNdXi( 3, 11 ) = 4.5 * zeta1 * ( 6.0 * zeta4 - 1.0 );
            tdNdXi( 3, 12 ) = 4.5 * zeta2 * ( 3.0 * zeta2 - 1.0 );
            tdNdXi( 3, 13 ) = 4.5 * zeta2 * ( 6.0 * zeta4 - 1.0 );
            tdNdXi( 3, 14 ) = 4.5 * zeta3 * ( 3.0 * zeta3 - 1.0 );
            tdNdXi( 3, 15 ) = 4.5 * zeta3 * ( 6.0 * zeta4 - 1.0 );
            tdNdXi( 3, 17 ) = 27.0 * zeta1 * zeta2;
            tdNdXi( 3, 18 ) = 27.0 * zeta2 * zeta3;
            tdNdXi( 3, 19 ) = 27.0 * zeta1 * zeta3;

            return tdNdXi;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::TET, Interpolation_Type::LAGRANGE, 3, 20 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 4, "TET20 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack the tet coordinates from input vector
            real zeta1 = aXi( 0 );
            real zeta2 = aXi( 1 );
            real zeta3 = aXi( 2 );
            real zeta4 = aXi( 3 );

            // populate output matrix
            Matrix< DDRMat > td2NdXi2( 10, 20, 0.0 );
            td2NdXi2( 0, 0 )  = 9.0 * ( 3.0 * zeta1 - 1.0 );
            td2NdXi2( 0, 4 )  = 27.0 * zeta2 ;
            td2NdXi2( 0, 8 )  = 27.0 * zeta3;
            td2NdXi2( 0, 10 ) = 27.0 * zeta4;

            td2NdXi2( 1, 1 )  = 9.0 * ( 3.0 * zeta2 - 1.0 );
            td2NdXi2( 1, 5 )  = 27.0 * zeta1;
            td2NdXi2( 1, 6 )  = 27.0 * zeta3;
            td2NdXi2( 1, 12 ) = 27.0 * zeta4 ;

            td2NdXi2( 2, 2 )  = 9.0 * ( 3.0 * zeta3 - 1.0 );
            td2NdXi2( 2, 7 )  = 27.0 * zeta2;
            td2NdXi2( 2, 9 )  = 27.0 * zeta1;
            td2NdXi2( 2, 14 ) = 27.0 * zeta4;

            td2NdXi2( 3, 3 )  = 9.0 * ( 3.0 * zeta4 - 1.0 );
            td2NdXi2( 3, 11 ) = 27.0 * zeta1;
            td2NdXi2( 3, 13 ) = 27.0 * zeta2;
            td2NdXi2( 3, 15 ) = 27.0 * zeta3;

            td2NdXi2( 4, 14 ) = 4.5 * ( 6.0 * zeta3 - 1.0 );
            td2NdXi2( 4, 15 ) = 4.5 * ( 6.0 * zeta4 - 1.0 );
            td2NdXi2( 4, 18 ) = 27.0 * zeta2;
            td2NdXi2( 4, 19 ) = 27.0 * zeta1;

            td2NdXi2( 5, 12 ) = 4.5 * ( 6.0 * zeta2 - 1.0 );
            td2NdXi2( 5, 13 ) = 4.5 * ( 6.0 * zeta4 - 1.0 );
            td2NdXi2( 5, 17 ) = 27.0 * zeta1;
            td2NdXi2( 5, 18 ) = 27.0 * zeta3;

            td2NdXi2( 6, 10 ) = 4.5 * ( 6.0 * zeta1 - 1.0 );
            td2NdXi2( 6, 11 ) = 4.5 * ( 6.0 * zeta4 - 1.0 );
            td2NdXi2( 6, 17 ) = 27.0 * zeta2;
            td2NdXi2( 6, 19 ) = 27.0 * zeta3;

            td2NdXi2( 7, 6 )  = 4.5 * ( 6.0 * zeta2 - 1.0 );
            td2NdXi2( 7, 7 )  = 4.5 * ( 6.0 * zeta3 - 1.0 );
            td2NdXi2( 7, 16 ) = 27.0 * zeta1;
            td2NdXi2( 7, 18 ) = 27.0 * zeta4;

            td2NdXi2( 8, 8 )  = 4.5 * ( 6.0 * zeta1 - 1.0 );
            td2NdXi2( 8, 9 )  = 4.5 * ( 6.0 * zeta3 - 1.0 );
            td2NdXi2( 8, 16 ) = 27.0 * zeta2;
            td2NdXi2( 8, 19 ) = 27.0 * zeta4;

            td2NdXi2( 9, 4 )  = 4.5 * ( 6.0 * zeta1 - 1.0 );
            td2NdXi2( 9, 5 )  = 4.5 * ( 6.0 * zeta2 - 1.0 );
            td2NdXi2( 9, 16 ) = 27.0 * zeta3;
            td2NdXi2( 9, 17 ) = 27.0 * zeta4;

            return td2NdXi2;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_TET20_HPP_ */
