/*
 * cl_FEM_Interpolation_Function_Quad4.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD4_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD4_HPP_

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
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::get_number_of_param_dimensions() const
        {
            return 2;
        }

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat( 2, 4 );
            tXiHat( 0, 0 ) = -1.000000;
            tXiHat( 1, 0 ) = -1.000000;
            tXiHat( 0, 1 ) =  1.000000;
            tXiHat( 1, 1 ) = -1.000000;
            tXiHat( 0, 2 ) =  1.000000;
            tXiHat( 1, 2 ) =  1.000000;
            tXiHat( 0, 3 ) = -1.000000;
            tXiHat( 1, 3 ) =  1.000000;
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD4 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // populate matrix with values
            Matrix< DDRMat> tN(1,4);
            tN( 0 ) = ( ( 1.0 - xi ) * ( 1.0 - eta ) ) * 0.25;
            tN( 1 ) = ( ( 1.0 + xi ) * ( 1.0 - eta ) ) * 0.25;
            tN( 2 ) = ( ( 1.0 + xi ) * ( 1.0 + eta ) ) * 0.25;
            tN( 3 ) = ( ( 1.0 - xi ) * ( 1.0 + eta ) ) * 0.25;
            return tN;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "QUAD4 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // populate output matrix
            Matrix< DDRMat > tdNdXi(2,4);
            tdNdXi( 0, 0 ) =  0.25 * ( eta - 1.0 );
            tdNdXi( 1, 0 ) =  0.25 * ( xi - 1.0 );

            tdNdXi( 0, 1 ) = -0.25 * ( eta - 1.0 );
            tdNdXi( 1, 1 ) = -0.25 * ( xi + 1.0 );

            tdNdXi( 0, 2 ) =  0.25 * ( eta + 1.0 );
            tdNdXi( 1, 2 ) =  0.25 * ( xi + 1.0 );

            tdNdXi( 0, 3 ) = -0.25 * ( eta + 1.0 );
            tdNdXi( 1, 3 ) = -0.25 * ( xi - 1.0 );
            return tdNdXi;
        }


//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                          "QUAD4 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // populate output matrix
            Matrix< DDRMat > td2NdXi2(3,4);
            td2NdXi2( 0, 0 ) =  0.0;
            td2NdXi2( 1, 0 ) =  0.0;
            td2NdXi2( 2, 0 ) =  0.25;

            td2NdXi2( 0, 1 ) =  0.0;
            td2NdXi2( 1, 1 ) =  0.0;
            td2NdXi2( 2, 1 ) = -0.25;

            td2NdXi2( 0, 2 ) =  0.0;
            td2NdXi2( 1, 2 ) =  0.0;
            td2NdXi2( 2, 2 ) =  0.25;

            td2NdXi2( 0, 3 ) =  0.0;
            td2NdXi2( 1, 3 ) =  0.0;
            td2NdXi2( 2, 3 ) = -0.25;
            return td2NdXi2;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD4_HPP_ */
