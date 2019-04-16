/*
 * cl_FEM_Interpolation_Function_Lagrange_Quad9.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD9_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD9_HPP_

#include "assert.h"

//#include "cl_FEM_Interpolation_Matrix.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src

// does not exist in femdoc

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::get_number_of_param_dimensions() const
        {
            return 2;
        }

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::QUADRATIC;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::get_param_coords() const
        {
            Matrix< DDRMat > tXiHat( 2, 9 );
            tXiHat( 0, 0 ) = -1.000000;
            tXiHat( 1, 0 ) = -1.000000;
            tXiHat( 0, 1 ) =  1.000000;
            tXiHat( 1, 1 ) = -1.000000;
            tXiHat( 0, 2 ) =  1.000000;
            tXiHat( 1, 2 ) =  1.000000;
            tXiHat( 0, 3 ) = -1.000000;
            tXiHat( 1, 3 ) =  1.000000;
            tXiHat( 0, 4 ) =  0.000000;
            tXiHat( 1, 4 ) = -1.000000;
            tXiHat( 0, 5 ) =  1.000000;
            tXiHat( 1, 5 ) =  0.000000;
            tXiHat( 0, 6 ) =  0.000000;
            tXiHat( 1, 6 ) =  1.000000;
            tXiHat( 0, 7 ) = -1.000000;
            tXiHat( 1, 7 ) =  0.000000;
            tXiHat( 0, 8 ) =  0.000000;
            tXiHat( 1, 8 ) =  0.000000;
            return tXiHat;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::eval_N(const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                          "QUAD9 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // often used constants
            auto    c = xi * eta * 0.25;
            auto  xi2 = std::pow(  xi, 2 );
            auto eta2 = std::pow( eta, 2 );

            // populate output matrix
            Matrix < DDRMat > tN(1,9);
            tN( 0 ) = ( c * ( eta - 1.0 ) * (xi - 1.0) );
            tN( 1 ) = ( c * ( eta - 1.0 ) * (xi + 1.0) );
            tN( 2 ) = ( c * ( eta + 1.0 ) * (xi + 1.0) );
            tN( 3 ) = ( c * ( eta + 1.0 ) * (xi - 1.0) );
            tN( 4 ) = ( eta * ( 1.0 - xi2 ) * ( eta - 1.0 ) ) * 0.5;
            tN( 5 ) = ( xi * ( 1.0 - eta2)*( xi + 1.0 ) )*0.5;
            tN( 6 ) = ( eta * (1.0 - xi2)*( eta + 1.0 ) )*0.5;
            tN( 7 ) = ( xi*( 1.0 - eta2 )*( xi - 1.0 ) )*0.5;
            tN( 8 ) = ( eta2 - 1.0 )*( xi2 - 1.0 );
            return tN;

        }
//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                          "QUAD9 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // often used constants
            auto    c = xi*eta;
            auto  xi2 = std::pow(  xi, 2 );
            auto eta2 = std::pow( eta, 2 );

            // populate output matrix
            Matrix< DDRMat > tdNdXi(2,9);
            tdNdXi( 0, 0 ) = ( eta * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            tdNdXi( 1, 0 ) =  ( xi * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) ) * 0.25;

            tdNdXi( 0, 1 ) = ( eta * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            tdNdXi( 1, 1 ) =  ( xi * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) ) * 0.25;

            tdNdXi( 0, 2 ) = ( eta * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            tdNdXi( 1, 2 ) =  ( xi * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) ) * 0.25;

            tdNdXi( 0, 3 ) = ( eta * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            tdNdXi( 1, 3 ) =  ( xi * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) ) * 0.25;

            tdNdXi( 0, 4 ) =  - c * ( eta - 1.0 );
            tdNdXi( 1, 4 ) =  -( ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) ) * 0.5;

            tdNdXi( 0, 5 ) = -( ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.5;
            tdNdXi( 1, 5 ) = - c * ( xi + 1.0 );

            tdNdXi( 0, 6 ) = - c * ( eta + 1.0 );
            tdNdXi( 1, 6 ) =  -( ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) ) * 0.5;

            tdNdXi( 0, 7 ) = -( ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.5;
            tdNdXi( 1, 7 ) = - c * ( xi - 1.0 );

            tdNdXi( 0, 8 ) = 2.0 * xi * ( eta2 - 1.0 );
            tdNdXi( 1, 8 ) = 2.0 * eta * ( xi2 - 1.0 );
            return tdNdXi;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 9 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                          "QUAD9 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto  xi = aXi( 0 );
            auto eta = aXi( 1 );

            // often used constants
            auto  xi2 = std::pow(  xi, 2 );
            auto eta2 = std::pow( eta, 2 );

            // populate output matrix
            Matrix< DDRMat > td2NdXi2(3,9);
            td2NdXi2( 0, 0 ) = ( eta * ( eta - 1.0 ) ) * 0.5;
            td2NdXi2( 1, 0 ) = ( xi * ( xi - 1.0 ) ) * 0.5;
            td2NdXi2( 2, 0 ) = ( ( 2.0 * eta - 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.25;

            td2NdXi2( 0, 1 ) = ( eta * ( eta - 1.0 ) ) * 0.5;
            td2NdXi2( 1, 1 ) = ( xi * ( xi + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 1 ) = ( ( 2.0 * eta - 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.25;

            td2NdXi2( 0, 2 ) = ( eta * ( eta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 2 ) = ( xi * ( xi + 1.0 ) ) * 0.5;
            td2NdXi2( 2, 2 ) = ( ( 2.0 * eta + 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.25;

            td2NdXi2( 0, 3 ) = ( eta * ( eta + 1.0 ) ) * 0.5;
            td2NdXi2( 1, 3 ) = ( xi * ( xi - 1.0 ) ) * 0.5;
            td2NdXi2( 2, 3 ) = ( ( 2.0 * eta + 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.25;

            td2NdXi2( 0, 4 ) = -eta * ( eta - 1.0 );
            td2NdXi2( 1, 4 ) = 1.0 - xi2;
            td2NdXi2( 2, 4 ) = -xi * ( 2.0 * eta - 1.0 );

            td2NdXi2( 0, 5 ) = 1.0 - eta2;
            td2NdXi2( 1, 5 ) = -xi * ( xi + 1.0 );
            td2NdXi2( 2, 5 ) = -eta * ( 2.0 * xi + 1.0 );

            td2NdXi2( 0, 6 ) = -eta * ( eta + 1.0 );
            td2NdXi2( 1, 6 ) = 1.0 - xi2;
            td2NdXi2( 2, 6 ) = -xi * ( 2.0 * eta + 1.0 );

            td2NdXi2( 0, 7 ) = 1.0 - eta2;
            td2NdXi2( 1, 7 ) = -xi * ( xi - 1.0 );
            td2NdXi2( 2, 7 ) = -eta * ( 2.0 * xi - 1.0 );

            td2NdXi2( 0, 8 ) = 2.0 * eta2 - 2.0;
            td2NdXi2( 1, 8 ) = 2.0 * xi2 - 2.0;
            td2NdXi2( 2, 8 ) = 4.0 * eta * xi;
            return td2NdXi2;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD9_HPP_ */
