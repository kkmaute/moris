/*
 * cl_FEM_Interpolation_Function_Lagrange_Hex27.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX27_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX27_HPP_

#include "assert.h"

#include "cl_FEM_Interpolation_Matrix.hpp"
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
    Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 27  >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::QUADRATIC;
    }

//------------------------------------------------------------------------------

    template<>
    void
    Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 27  >::get_param_coords(
            Matrix< DDRMat > & aXihat ) const
    {
        aXihat.set_size( 3, 27 );

        aXihat( 0,  0 ) = -1.000000;
        aXihat( 1,  0 ) = -1.000000;
        aXihat( 2,  0 ) = -1.000000;
        aXihat( 0,  1 ) =  1.000000;
        aXihat( 1,  1 ) = -1.000000;
        aXihat( 2,  1 ) = -1.000000;
        aXihat( 0,  2 ) =  1.000000;
        aXihat( 1,  2 ) =  1.000000;
        aXihat( 2,  2 ) = -1.000000;
        aXihat( 0,  3 ) = -1.000000;
        aXihat( 1,  3 ) =  1.000000;
        aXihat( 2,  3 ) = -1.000000;
        aXihat( 0,  4 ) = -1.000000;
        aXihat( 1,  4 ) = -1.000000;
        aXihat( 2,  4 ) =  1.000000;
        aXihat( 0,  5 ) =  1.000000;
        aXihat( 1,  5 ) = -1.000000;
        aXihat( 2,  5 ) =  1.000000;
        aXihat( 0,  6 ) =  1.000000;
        aXihat( 1,  6 ) =  1.000000;
        aXihat( 2,  6 ) =  1.000000;
        aXihat( 0,  7 ) = -1.000000;
        aXihat( 1,  7 ) =  1.000000;
        aXihat( 2,  7 ) =  1.000000;
        aXihat( 0,  8 ) =  0.000000;
        aXihat( 1,  8 ) = -1.000000;
        aXihat( 2,  8 ) = -1.000000;
        aXihat( 0,  9 ) =  1.000000;
        aXihat( 1,  9 ) =  0.000000;
        aXihat( 2,  9 ) = -1.000000;
        aXihat( 0, 10 ) =  0.000000;
        aXihat( 1, 10 ) =  1.000000;
        aXihat( 2, 10 ) = -1.000000;
        aXihat( 0, 11 ) = -1.000000;
        aXihat( 1, 11 ) =  0.000000;
        aXihat( 2, 11 ) = -1.000000;
        aXihat( 0, 12 ) = -1.000000;
        aXihat( 1, 12 ) = -1.000000;
        aXihat( 2, 12 ) =  0.000000;
        aXihat( 0, 13 ) =  1.000000;
        aXihat( 1, 13 ) = -1.000000;
        aXihat( 2, 13 ) =  0.000000;
        aXihat( 0, 14 ) =  1.000000;
        aXihat( 1, 14 ) =  1.000000;
        aXihat( 2, 14 ) =  0.000000;
        aXihat( 0, 15 ) = -1.000000;
        aXihat( 1, 15 ) =  1.000000;
        aXihat( 2, 15 ) =  0.000000;
        aXihat( 0, 16 ) =  0.000000;
        aXihat( 1, 16 ) = -1.000000;
        aXihat( 2, 16 ) =  1.000000;
        aXihat( 0, 17 ) =  1.000000;
        aXihat( 1, 17 ) =  0.000000;
        aXihat( 2, 17 ) =  1.000000;
        aXihat( 0, 18 ) =  0.000000;
        aXihat( 1, 18 ) =  1.000000;
        aXihat( 2, 18 ) =  1.000000;
        aXihat( 0, 19 ) = -1.000000;
        aXihat( 1, 19 ) =  0.000000;
        aXihat( 2, 19 ) =  1.000000;
        aXihat( 0, 20 ) =  0.000000;
        aXihat( 1, 20 ) =  0.000000;
        aXihat( 2, 20 ) =  0.000000;
        aXihat( 0, 21 ) =  0.000000;
        aXihat( 1, 21 ) =  0.000000;
        aXihat( 2, 21 ) = -1.000000;
        aXihat( 0, 22 ) =  0.000000;
        aXihat( 1, 22 ) =  0.000000;
        aXihat( 2, 22 ) =  1.000000;
        aXihat( 0, 23 ) = -1.000000;
        aXihat( 1, 23 ) =  0.000000;
        aXihat( 2, 23 ) =  0.000000;
        aXihat( 0, 24 ) =  1.000000;
        aXihat( 1, 24 ) =  0.000000;
        aXihat( 2, 24 ) =  0.000000;
        aXihat( 0, 25 ) =  0.000000;
        aXihat( 1, 25 ) = -1.000000;
        aXihat( 2, 25 ) =  0.000000;
        aXihat( 0, 26 ) =  0.000000;
        aXihat( 1, 26 ) =  1.000000;
        aXihat( 2, 26 ) =  0.000000;
    }

//------------------------------------------------------------------------------
        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 27  >::eval_N(
                      Interpolation_Matrix  & aN,
                const Matrix< DDRMat > & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( aN.n_cols() == 27,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( aN.n_rows() == 1,
                    "eval_shape: aN not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto    xi = aXi( 0 );
            auto   eta = aXi( 1 );
            auto  zeta = aXi( 2 );

            // often used constants
            auto    xi2 = std::pow(   xi, 2 );
            auto   eta2 = std::pow(  eta, 2 );
            auto  zeta2 = std::pow( zeta, 2 );

            auto a = -0.25 * eta * zeta;
            auto b = -0.25 * xi * zeta;
            auto c = -0.25 * xi * eta;
            auto d = 0.125 * xi * eta * zeta;

            // populate output matrix
            aN(  0 ) = d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
            aN(  1 ) = d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
            aN(  2 ) = d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
            aN(  3 ) = d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
            aN(  4 ) = d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
            aN(  5 ) = d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
            aN(  6 ) = d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
            aN(  7 ) = d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
            aN(  8 ) = a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
            aN(  9 ) = b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
            aN( 10 ) = a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
            aN( 11 ) = b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
            aN( 12 ) = c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );
            aN( 13 ) = c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );
            aN( 14 ) = c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );
            aN( 15 ) = c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );
            aN( 16 ) = a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
            aN( 17 ) = b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
            aN( 18 ) = a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
            aN( 19 ) = b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
            aN( 20 ) = -( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 );
            aN( 21 ) = ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            aN( 22 ) = ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            aN( 23 ) = ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.5;
            aN( 24 ) = ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.5;
            aN( 25 ) = ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.5;
            aN( 26 ) = ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.5;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 27  >::eval_dNdXi(
                      Interpolation_Matrix  & adNdXi,
                const Matrix< DDRMat > & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( adNdXi.n_cols() == 27,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( adNdXi.n_rows() == 3,
                    "eval_shape: aN not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // often used constants
            auto   xi2 = std::pow(   xi, 2 );
            auto  eta2 = std::pow(  eta, 2 );
            auto zeta2 = std::pow( zeta, 2 );

            auto a = 0.125*eta*zeta;
            auto b = 0.125*xi*zeta;
            auto c = 0.125*xi*eta;
            auto d = -0.5*xi*eta*zeta;

            adNdXi( 0,  0 ) = a * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
            adNdXi( 1,  0 ) = b * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
            adNdXi( 2,  0 ) = c * ( 2.0 * zeta - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );


            adNdXi( 0,  1 ) = a * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
            adNdXi( 1,  1 ) = b * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
            adNdXi( 2,  1 ) = c * ( 2.0 * zeta - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );

            adNdXi( 0,  2 ) = a * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
            adNdXi( 1,  2 ) = b * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
            adNdXi( 2,  2 ) = c * ( 2.0 * zeta - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );

            adNdXi( 0,  3 ) = a * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
            adNdXi( 1,  3 ) = b * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
            adNdXi( 2,  3 ) = c * ( 2.0 * zeta - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );

            adNdXi( 0,  4 ) = a * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
            adNdXi( 1,  4 ) = b * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
            adNdXi( 2,  4 ) = c * ( 2.0 * zeta + 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );

            adNdXi( 0,  5 ) = a * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
            adNdXi( 1,  5 ) = b * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
            adNdXi( 2,  5 ) = c * ( 2.0 * zeta + 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );

            adNdXi( 0,  6 ) = a * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
            adNdXi( 1,  6 ) = b * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
            adNdXi( 2,  6 ) = c * ( 2.0 * zeta + 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );

            adNdXi( 0,  7 ) = a * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
            adNdXi( 1,  7 ) = b * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
            adNdXi( 2,  7 ) = c * ( 2.0 * zeta + 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );

            adNdXi( 0,  8 ) = d * ( eta - 1.0 ) * ( zeta - 1.0 );
            adNdXi( 1,  8 ) = - ( zeta * ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            adNdXi( 2,  8 ) = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( eta - 1.0 ) ) * 0.25;

            adNdXi( 0,  9 ) = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            adNdXi( 1,  9 ) = d * ( xi + 1.0 ) * ( zeta - 1.0 );
            adNdXi( 2,  9 ) = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( xi + 1.0 ) ) * 0.25;


            adNdXi( 0, 10 ) = d * ( eta + 1.0 ) * ( zeta - 1.0 );
            adNdXi( 1, 10 ) =  - ( zeta * ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            adNdXi( 2, 10 ) = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( eta + 1.0 ) ) * 0.25;

            adNdXi( 0, 11 ) = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            adNdXi( 1, 11 ) = d * ( xi - 1.0 ) * ( zeta - 1.0 );
            adNdXi( 2, 11 ) = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( xi - 1.0 ) ) * 0.25;

            adNdXi( 0, 12 ) = - ( eta * ( 2.0 * xi - 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            adNdXi( 1, 12 ) = - ( xi * ( 2.0 * eta - 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
            adNdXi( 2, 12 ) = d * ( eta - 1.0 ) * ( xi - 1.0 );

            adNdXi( 0, 13 ) = - ( eta * ( 2.0 * xi + 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
            adNdXi( 1, 13 ) = - ( xi * ( 2.0 * eta - 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
            adNdXi( 2, 13 ) = d * ( eta - 1.0 ) * ( xi + 1.0 );

            adNdXi( 0, 14 ) = - ( eta * ( 2.0 * xi + 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            adNdXi( 1, 14 ) = - ( xi * ( 2.0 * eta + 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
            adNdXi( 2, 14 ) = d * ( eta + 1.0 ) * ( xi + 1.0 );

            adNdXi( 0, 15 ) = - ( eta * ( 2.0 * xi - 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
            adNdXi( 1, 15 ) = - ( xi * ( 2.0 * eta + 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
            adNdXi( 2, 15 ) = d * ( eta + 1.0 ) * ( xi - 1.0 );

            adNdXi( 0, 16 ) = d * ( eta - 1.0 ) * ( zeta + 1.0 );
            adNdXi( 1, 16 ) = - ( zeta * ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            adNdXi( 2, 16 ) = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( eta - 1.0 ) ) * 0.25;


            adNdXi( 0, 17 ) = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            adNdXi( 1, 17 ) = d * ( xi + 1.0 ) * ( zeta + 1.0 );
            adNdXi( 2, 17 ) = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( xi + 1.0 ) ) * 0.25;

            adNdXi( 0, 18 ) = d * ( eta + 1.0 ) * ( zeta + 1.0 );
            adNdXi( 1, 18 ) = - ( zeta * ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            adNdXi( 2, 18 ) = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( eta + 1.0 ) ) * 0.25;

            adNdXi( 0, 19 ) = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            adNdXi( 1, 19 ) = d * ( xi - 1.0 ) * ( zeta + 1.0 );
            adNdXi( 2, 19 ) = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( xi - 1.0 ) ) * 0.25;

            adNdXi( 0, 20 ) = - 2.0 * xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 );
            adNdXi( 1, 20 ) = - 2.0 * eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 );
            adNdXi( 2, 20 ) = - 2.0 * zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 );

            adNdXi( 0, 21 ) = 8.0 * b * ( eta2 - 1.0 ) * ( zeta - 1.0 );
            adNdXi( 1, 21 ) = 8.0 * a * ( xi2 - 1.0 ) * ( zeta - 1.0 );
            adNdXi( 2, 21 ) = ( ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.5;

            adNdXi( 0, 22 ) = 8.0 * b * ( eta2 - 1.0 ) * ( zeta + 1.0 );
            adNdXi( 1, 22 ) = 8.0 * a * ( xi2 - 1.0 ) * ( zeta + 1.0 );
            adNdXi( 2, 22 ) = ( ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.5;

            adNdXi( 0, 23 ) =( ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
            adNdXi( 1, 23 ) = 8.0 * c * ( zeta2 - 1.0 ) * ( xi - 1.0 );
            adNdXi( 2, 23 ) = 8.0 * b * ( eta2 - 1.0 ) * ( xi - 1.0 );

            adNdXi( 0, 24 ) =( ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
            adNdXi( 1, 24 ) = 8.0 * c * ( zeta2 - 1.0 ) * ( xi + 1.0 );
            adNdXi( 2, 24 ) = 8.0 * b * ( eta2 - 1.0 ) * ( xi + 1.0 );

            adNdXi( 0, 25 ) = 8.0 * c * ( zeta2 - 1.0 ) * ( eta - 1.0 );
            adNdXi( 1, 25 ) = ( ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
            adNdXi( 2, 25 ) = 8.0 * a * ( xi2 - 1.0 ) * ( eta - 1.0 );

            adNdXi( 0, 26 ) = 8.0 * c * ( zeta2 - 1.0 ) * ( eta + 1.0 );
            adNdXi( 1, 26 ) = ( ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
            adNdXi( 2, 26 ) = 8.0 * a * ( xi2 - 1.0 ) * ( eta + 1.0 );

        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< Interpolation_Type::LAGRANGE, 3, 27  >::eval_d2NdXi2(
                      Interpolation_Matrix  & ad2NdXi2,
                const Matrix< DDRMat > & aXi
        ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3,
                    "eval_shape: aXi not allocated or hat wrong size." );

            // make sure that output array has correct number of columns
            MORIS_ASSERT( ad2NdXi2.n_cols() == 27,
                    "eval_shape: aN not allocated or hat wrong size." );

            // make sure that output array has correct number of rows
            MORIS_ASSERT( ad2NdXi2.n_rows() == 6,
                    "eval_shape: aN not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto   xi = aXi( 0 );
            auto  eta = aXi( 1 );
            auto zeta = aXi( 2 );

            // often used constants
            auto   xi2 = std::pow(   xi, 2 );
            auto  eta2 = std::pow(  eta, 2 );
            auto zeta2 = std::pow( zeta, 2 );

            auto a = eta * zeta;
            auto b = xi * zeta;
            auto c = xi * eta;
            auto d = 2 * xi* eta * zeta;

            ad2NdXi2( 0,  0 ) = ( a * (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 1,  0 ) = ( b * (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 2,  0 ) = ( c * (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            ad2NdXi2( 3,  0 ) = (   xi * ( 2.0 *  eta - 1.0 ) * ( 2.0 * zeta - 1.0 ) * (   xi - 1.0 ) ) * 0.125;
            ad2NdXi2( 4,  0 ) = (  eta * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta - 1.0 ) * (  eta - 1.0 ) ) * 0.125;
            ad2NdXi2( 5,  0 ) = ( zeta * ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.125;

            ad2NdXi2( 0,  1 ) = ( a * (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 1,  1 ) = ( b * (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 2,  1 ) = ( c * (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            ad2NdXi2( 3,  1 ) = (   xi * ( 2.0 *  eta - 1.0 ) * ( 2.0 * zeta - 1.0 ) * (   xi + 1.0 ) ) * 0.125;
            ad2NdXi2( 4,  1 ) = (  eta * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta - 1.0 ) * (  eta - 1.0 ) ) * 0.125;
            ad2NdXi2( 5,  1 ) = ( zeta * ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.125;

            ad2NdXi2( 0,  2 ) = ( a * (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 1,  2 ) = ( b * (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 2,  2 ) = ( c * (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            ad2NdXi2( 3,  2 ) = (   xi * ( 2.0 *  eta + 1.0 ) * ( 2.0 * zeta - 1.0 ) * (   xi + 1.0 ) ) * 0.125;
            ad2NdXi2( 4,  2 ) = (  eta * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta - 1.0 ) * (  eta + 1.0 ) ) * 0.125;
            ad2NdXi2( 5,  2 ) = ( zeta * ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.125;

            ad2NdXi2( 0,  3 ) = ( a * (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 1,  3 ) = ( b * (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 2,  3 ) = ( c * (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            ad2NdXi2( 3,  3 ) = (   xi * ( 2.0 *  eta + 1.0 ) * ( 2.0 * zeta - 1.0 ) * (   xi - 1.0 ) ) * 0.125;
            ad2NdXi2( 4,  3 ) = (  eta * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta - 1.0 ) * (  eta + 1.0 ) ) * 0.125;
            ad2NdXi2( 5,  3 ) = ( zeta * ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.125;

            ad2NdXi2( 0,  4 ) = ( a * (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 1,  4 ) = ( b * (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 2,  4 ) = ( c * (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            ad2NdXi2( 3,  4 ) = (   xi * ( 2.0 *  eta - 1.0 ) * ( 2.0 * zeta + 1.0 ) * (   xi - 1.0 ) ) * 0.125;
            ad2NdXi2( 4,  4 ) = (  eta * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta + 1.0 ) * (  eta - 1.0 ) ) * 0.125;
            ad2NdXi2( 5,  4 ) = ( zeta * ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.125;

            ad2NdXi2( 0,  5 ) = ( a * (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 1,  5 ) = ( b * (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 2,  5 ) = ( c * (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            ad2NdXi2( 3,  5 ) = (   xi * ( 2.0 *  eta - 1.0 ) * ( 2.0 * zeta + 1.0 ) * (   xi + 1.0 ) ) * 0.125;
            ad2NdXi2( 4,  5 ) = (  eta * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta + 1.0 ) * (  eta - 1.0 ) ) * 0.125;
            ad2NdXi2( 5,  5 ) = ( zeta * ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.125;

            ad2NdXi2( 0,  6 ) = ( a * (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 1,  6 ) = ( b * (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 2,  6 ) = ( c * (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.25;
            ad2NdXi2( 3,  6 ) = (   xi * ( 2.0 *  eta + 1.0 ) * ( 2.0 * zeta + 1.0 ) * (   xi + 1.0 ) ) * 0.125;
            ad2NdXi2( 4,  6 ) = (  eta * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta + 1.0 ) * (  eta + 1.0 ) ) * 0.125;
            ad2NdXi2( 5,  6 ) = ( zeta * ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.125;

            ad2NdXi2( 0,  7 ) = ( a * (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 1,  7 ) = ( b * (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 2,  7 ) = ( c * (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.25;
            ad2NdXi2( 3,  7 ) = (   xi * ( 2.0 *  eta + 1.0 ) * ( 2.0 * zeta + 1.0 ) * (   xi - 1.0 ) ) * 0.125;
            ad2NdXi2( 4,  7 ) = (  eta * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta + 1.0 ) * (  eta + 1.0 ) ) * 0.125;
            ad2NdXi2( 5,  7 ) = ( zeta * ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.125;

            ad2NdXi2( 0,  8 ) = - ( a * (  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1,  8 ) = - ( zeta * (   xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 2,  8 ) = - (  eta * (   xi2 - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 3,  8 ) = - ( ( 2.0 *  eta - 1.0 ) * (   xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 4,  8 ) = - ( c * ( 2.0 * zeta - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 5,  8 ) = - ( b * ( 2.0 *  eta - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;

            ad2NdXi2( 0,  9 ) = - ( zeta * (  eta2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1,  9 ) = - ( b * (   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 2,  9 ) = - (   xi * (  eta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 3,  9 ) = - ( c * ( 2.0 * zeta - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 4,  9 ) = - ( (  eta2 - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 5,  9 ) = - ( a * ( 2.0 *   xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 10 ) = - ( a * (  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 10 ) = - ( zeta * (   xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 10 ) = - (  eta * (   xi2 - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 10 ) = - ( ( 2.0 *  eta + 1.0 ) * (   xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 4, 10 ) = - ( c * ( 2.0 * zeta - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 10 ) = - ( b * ( 2.0 *  eta + 1.0 ) * ( zeta - 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 11 ) = - ( zeta * (  eta2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 11 ) = - ( b * (   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 11 ) = - (   xi * (  eta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 11 ) = - ( c * ( 2.0 * zeta - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 11 ) = - ( (  eta2 - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.25;
            ad2NdXi2( 5, 11 ) = - ( a * ( 2.0 *   xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 12 ) = - (  eta * ( zeta2 - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 12 ) = - (   xi * ( zeta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 12 ) = - ( c * (  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 12 ) = - ( b * ( 2.0 *  eta - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 12 ) = - ( a * ( 2.0 *   xi - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 12 ) = - ( ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 13 ) = - (  eta * ( zeta2 - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 13 ) = - (   xi * ( zeta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 13 ) = - ( c * (  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 13 ) = - ( b * ( 2.0 *  eta - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 13 ) = - ( a * ( 2.0 *   xi + 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 13 ) = - ( ( 2.0 *  eta - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta2 - 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 14 ) = - (  eta * ( zeta2 - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 14 ) = - (   xi * ( zeta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 14 ) = - ( c * (  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 14 ) = - ( b * ( 2.0 *  eta + 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 14 ) = - ( a * ( 2.0 *   xi + 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 14 ) = - ( ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( zeta2 - 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 15 ) = - (  eta * ( zeta2 - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 15 ) = - (   xi * ( zeta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 15 ) = - ( c * (  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 15 ) = - ( b * ( 2.0 *  eta + 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 15 ) = - ( a * ( 2.0 *   xi - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 15 ) = - ( ( 2.0 *  eta + 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.25;

            ad2NdXi2( 0, 16 ) = - ( a * (  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 16 ) = - ( zeta * (   xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 16 ) = - (  eta * (   xi2 - 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 16 ) = - ( ( 2.0 *  eta - 1.0 ) * (   xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 4, 16 ) =  - ( c * ( 2.0 * zeta + 1.0 ) * (  eta - 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 16 ) =  - ( b * ( 2.0 *  eta - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 17 ) = - ( zeta * (  eta2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 17 ) = - ( b * (   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 17 ) = - (   xi * (  eta2 - 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 17 ) = - ( c * ( 2.0 * zeta + 1.0 ) * (   xi + 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 17 ) =  - ( (  eta2 - 1.0 ) * ( 2.0 *   xi + 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 5, 17 ) =  - ( a * ( 2.0 *   xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 18 ) = - ( a * (  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 18 ) = - ( zeta * (   xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 18 ) = - (  eta * (   xi2 - 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 18 ) = - ( ( 2.0 *  eta + 1.0 ) * (   xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 4, 18 ) = - ( c * ( 2.0 * zeta + 1.0 ) * (  eta + 1.0 ) ) * 0.5;
            ad2NdXi2( 5, 18 ) = - ( b * ( 2.0 *  eta + 1.0 ) * ( zeta + 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 19 ) = - ( zeta * (  eta2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 1, 19 ) = - ( b * (   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            ad2NdXi2( 2, 19 ) = - (   xi * (  eta2 - 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 3, 19 ) = - ( c * ( 2.0 * zeta + 1.0 ) * (   xi - 1.0 ) ) * 0.5;
            ad2NdXi2( 4, 19 ) =  - ( (  eta2 - 1.0 ) * ( 2.0 *   xi - 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.25;
            ad2NdXi2( 5, 19 ) =  - ( a * ( 2.0 *   xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;

            ad2NdXi2( 0, 20 ) = - 2.0 * (  eta2 - 1.0 ) * ( zeta2 - 1.0 );
            ad2NdXi2( 1, 20 ) = - 2.0 * (   xi2 - 1.0 ) * ( zeta2 - 1.0 );
            ad2NdXi2( 2, 20 ) = - 2.0 * (  eta2 - 1.0 ) * (   xi2 - 1.0 );
            ad2NdXi2( 3, 20 ) = - 4.0 * a * (   xi2 - 1.0 );
            ad2NdXi2( 4, 20 ) = - 4.0 * b * (  eta2 - 1.0 );
            ad2NdXi2( 5, 20 ) = - 4.0 * c * ( zeta2 - 1.0 );

            ad2NdXi2( 0, 21 ) = zeta * (  eta2 - 1.0 ) * ( zeta - 1.0 );
            ad2NdXi2( 1, 21 ) = zeta * (   xi2 - 1.0 ) * ( zeta - 1.0 );
            ad2NdXi2( 2, 21 ) = (  eta2 - 1.0 ) * (   xi2 - 1.0 );
            ad2NdXi2( 3, 21 ) =  eta * (   xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 );
            ad2NdXi2( 4, 21 ) =   xi * (  eta2 - 1.0 ) * ( 2.0 * zeta - 1.0 );
            ad2NdXi2( 5, 21 ) = d * ( zeta - 1.0 );

            ad2NdXi2( 0, 22 ) = zeta * (  eta2 - 1.0 ) * ( zeta + 1.0 );
            ad2NdXi2( 1, 22 ) = zeta * (   xi2 - 1.0 ) * ( zeta + 1.0 );
            ad2NdXi2( 2, 22 ) = (  eta2 - 1.0 ) * (   xi2 - 1.0 );
            ad2NdXi2( 3, 22 ) =  eta * (   xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 );
            ad2NdXi2( 4, 22 ) =   xi * (  eta2 - 1.0 ) * ( 2.0 * zeta + 1.0 );
            ad2NdXi2( 5, 22 ) = d * ( zeta + 1.0 );

            ad2NdXi2( 0, 23 ) = (  eta2 - 1.0 ) * ( zeta2 - 1.0 );
            ad2NdXi2( 1, 23 ) =   xi * ( zeta2 - 1.0 ) * (   xi - 1.0 );
            ad2NdXi2( 2, 23 ) =   xi * (  eta2 - 1.0 ) * (   xi - 1.0 );
            ad2NdXi2( 3, 23 ) = d * (   xi - 1.0 );
            ad2NdXi2( 4, 23 ) = zeta * (  eta2 - 1.0 ) * ( 2.0 *   xi - 1.0 );
            ad2NdXi2( 5, 23 ) =  eta * ( 2.0 *   xi - 1.0 ) * ( zeta2 - 1.0 );

            ad2NdXi2( 0, 24 ) = (  eta2 - 1.0 ) * ( zeta2 - 1.0 );
            ad2NdXi2( 1, 24 ) =   xi * ( zeta2 - 1.0 ) * (   xi + 1.0 );
            ad2NdXi2( 2, 24 ) =   xi * (  eta2 - 1.0 ) * (   xi + 1.0 );
            ad2NdXi2( 3, 24 ) = d * (   xi + 1.0 );
            ad2NdXi2( 4, 24 ) = zeta * (  eta2 - 1.0 ) * ( 2.0 *   xi + 1.0 );
            ad2NdXi2( 5, 24 ) =  eta * ( 2.0 *   xi + 1.0 ) * ( zeta2 - 1.0 );

            ad2NdXi2( 0, 25 ) =  eta * ( zeta2 - 1.0 ) * (  eta - 1.0 );
            ad2NdXi2( 1, 25 ) = (   xi2 - 1.0 ) * ( zeta2 - 1.0 );
            ad2NdXi2( 2, 25 ) =  eta * (   xi2 - 1.0 ) * (  eta - 1.0 );
            ad2NdXi2( 3, 25 ) = zeta * ( 2.0 *  eta - 1.0 ) * (   xi2 - 1.0 );
            ad2NdXi2( 4, 25 ) = d * (  eta - 1.0 );
            ad2NdXi2( 5, 25 ) =   xi * ( 2.0 *  eta - 1.0 ) * ( zeta2 - 1.0 );

            ad2NdXi2( 0, 26 ) =  eta * ( zeta2 - 1.0 ) * (  eta + 1.0 );
            ad2NdXi2( 1, 26 ) = (   xi2 - 1.0 ) * ( zeta2 - 1.0 );
            ad2NdXi2( 2, 26 ) =  eta * (   xi2 - 1.0 ) * (  eta + 1.0 );
            ad2NdXi2( 3, 26 ) = zeta * ( 2.0 *  eta + 1.0 ) * (   xi2 - 1.0 );
            ad2NdXi2( 4, 26 ) = d * (  eta + 1.0 );
            ad2NdXi2( 5, 26 ) =   xi * ( 2.0 *  eta + 1.0 ) * ( zeta2 - 1.0 );
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */



#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_HEX27_HPP_ */
