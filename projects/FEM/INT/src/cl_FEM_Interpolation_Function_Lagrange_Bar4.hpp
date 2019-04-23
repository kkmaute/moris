/*
 * cl_FEM_Interpolation_Function_Lagrange_Bar4.hpp
 *
 *  Created on: Apr 11, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR4_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR4_HPP_

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
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::get_number_of_param_dimensions() const
        {
            return 1;
        }

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::get_interpolation_order()  const
        {
            return mtk::Interpolation_Order::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat(1,4);
            tXiHat( 0 ) = -1.000000;
            tXiHat( 1 ) =  1.000000;
            tXiHat( 2 ) = -1.0/3.0;
            tXiHat( 3 ) =  1.0/3.0;
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                          "LINE4 - eval_N: aXi not allocated or hat wrong size." );

            real xi = aXi( 0 );
            real t116 = 1.0 / 16.0;

            Matrix< DDRMat > tN(1,4);
            tN( 0 ) = t116 * ( 3.0 * xi + 1.0 ) * ( 3.0 * xi - 1.0 ) * ( 1.0 - xi );
            tN( 1 ) = t116 * ( 3.0 * xi + 1.0 ) * ( 3.0 * xi - 1.0 ) * ( 1.0 + xi );
            tN( 2 ) = 9.0 * t116 * ( xi + 1.0 ) * ( 3.0 * xi - 1.0 ) * ( xi - 1.0 );
            tN( 3 ) = 9.0 * t116 * ( xi + 1.0 ) * ( 3.0 * xi + 1.0 ) * ( 1.0 - xi );
            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                          "LINE4 - eval_dNdXi: aXi not allocated or hat wrong size." );

            real xi  = aXi( 0 );
            real xi2 = std::pow( xi, 2 );
            real t116 = 1.0 / 16.0;

            Matrix< DDRMat > tdNdXi(1,4);
            tdNdXi( 0 ) = t116 * ( -27.0 * xi2 + 18.0 * xi + 1.0 );
            tdNdXi( 1 ) = t116 * (  27.0 * xi2 + 18.0 * xi - 1.0 );
            tdNdXi( 2 ) = 9.0 * t116 * (  9.0 * xi2 - 2.0 * xi - 3.0 );
            tdNdXi( 3 ) = 9.0 * t116 * ( -9.0 * xi2 - 2.0 * xi + 3.0 );
            return tdNdXi;

        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::LINE, Interpolation_Type::LAGRANGE, 1, 4 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1,
                          "LINE4 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            real xi  = aXi( 0 );
            real t98 = 9.0 / 8.0;

            Matrix< DDRMat > td2NdXi2(1,4);
            td2NdXi2( 0 ) =  t98 * ( 1.0 - 3.0 * xi );
            td2NdXi2( 1 ) =  t98 * ( 1.0 + 3.0 * xi );
            td2NdXi2( 2 ) =  t98 * ( 9.0 * xi - 1.0 );
            td2NdXi2( 2 ) = -t98 * ( 9.0 * xi + 1.0 );
            return td2NdXi2;
        }
//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR4_HPP_ */
