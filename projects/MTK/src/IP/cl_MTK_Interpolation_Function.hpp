/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_HPP_

#include "assert.hpp"

#include "cl_MTK_Enums.hpp"                          //MTK/src
#include "cl_MTK_Interpolation_Function_Base.hpp"    //MTK/src
#include "cl_Matrix.hpp"                             //LINALG/src

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        /**
         * shape function templated class
         * G : Geometry
         * T : Type
         * N : Dimension
         * B : Number of Basis
         */
        template< Geometry_Type G, Interpolation_Type T, uint N, uint B >
        class Interpolation_Function : public Interpolation_Function_Base
        {
            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            // default constructor
            Interpolation_Function(){};

            //------------------------------------------------------------------------------

            /**
             * evaluates the shape function at a given point
             *
             * @param[ out ] aN  shape function as
             *                   ( 1 x <number of nodes> )
             *
             * @param[ in ]  aXi parameter coordinates
             *                   ( <number of dimensions>  x 1 )
             */
            void eval_N(
                    const Matrix< DDRMat >& aXi,
                    Matrix< DDRMat >&       aNXi ) const;

            //------------------------------------------------------------------------------

            /**
             * calculates the first derivative of the shape function
             * in parameter space
             *
             * @param[ out ] adNdXi ( <number of dimensions> x <number of nodes> )
             *
             * @param[ in ] aXi    point where function is evaluated
             *                     ( <number of dimensions>  x 1 )
             *
             */

            void eval_dNdXi(
                    const Matrix< DDRMat >& aXi,
                    Matrix< DDRMat >&       adNdXi ) const;

            //------------------------------------------------------------------------------

            /**
             * calculates the second derivative of the shape function
             * in parameter space
             *
             * @param[ in ] ad2NdXi2 ( <number of 2nd order derivatives> x <number of nodes> )
             *
             * @param[ in ] aXi    point where function is evaluated
             *                     ( <number of dimensions>  x 1 )
             *
             * @param[ in ] ad2NdXi2 ( 3 or 6 x <number of nodes> )
             *                     matrix is ordered as follows:
             *             2D      1. row: d2N over dXi_1 dXi_1
             *                     2. row: d2N over dXi_2 dXi_2
             *                     3. row: d2N over dXi_1 dXi_2
             *
             *             3D      1. row: d2N over dXi_1 dXi_1
             *                     2. row: d2N over dXi_2 dXi_2
             *                     3. row: d2N over dXi_3 dXi_3
             *                     4. row: d2N over dXi_2 dXi_3
             *                     5. row: d2N over dXi_1 dXi_3
             *                     6. row: d2N over dXi_1 dXi_2
             */
            void eval_d2NdXi2(
                    const Matrix< DDRMat >& aXi,
                    Matrix< DDRMat >&       ad2NdXi2 ) const;

            //------------------------------------------------------------------------------

            /**
             * calculates the third derivatives of the shape function
             * in parameter space
             *
             * @param[ in ] ad3NdXi3 ( <number of 3rd order derivatives> x <number of nodes> )
             *
             * @param[ in ] aXi    point where function is evaluated
             *                     ( <number of dimensions>  x 1 )
             *
             */
            void eval_d3NdXi3(
                    const Matrix< DDRMat >& aXi,
                    Matrix< DDRMat >&       ad3NdXi3 ) const;

            //------------------------------------------------------------------------------

            /**
             * returns a matrix containing the parameter coordinates
             * < number of dimensions * number of basis >
             */
            void get_param_coords( Matrix< DDRMat >& aXiHat ) const;

            //------------------------------------------------------------------------------

            /**
             * returns the number of bases for this shape function
             */
            uint
            get_number_of_bases() const
            {
                return B;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the number of dimensions for this shape function
             */
            uint
            get_number_of_dimensions() const
            {
                return N;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the number of dimensions for this shape function
             */
            uint get_number_of_param_dimensions() const;

            //------------------------------------------------------------------------------

            /**
             * returns the interpolation order
             */
            Interpolation_Order get_interpolation_order() const;

            //------------------------------------------------------------------------------

            /**
             * returns the interpolation type
             */
            Interpolation_Type
            get_interpolation_type() const
            {
                return T;
            }

            //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------

        template< Geometry_Type G, Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< G, T, N, B >::eval_N(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       aNXi ) const
        {
            MORIS_ERROR( false, "eval_N not implemented for this interpolation function" );
        }

        //------------------------------------------------------------------------------

        template< Geometry_Type G, Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< G, T, N, B >::eval_dNdXi(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       adNdXi ) const
        {
            MORIS_ERROR( false, "eval_dNdXi not implemented for this interpolation function" );
        }

        //------------------------------------------------------------------------------

        template< Geometry_Type G, Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< G, T, N, B >::eval_d2NdXi2(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad2NdXi2 ) const
        {
            MORIS_ERROR( false, "eval_d2NdXi2 not implemented for this interpolation function" );
        }

        //------------------------------------------------------------------------------

        template< Geometry_Type G, Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< G, T, N, B >::eval_d3NdXi3(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       ad3NdXi3 ) const
        {
            MORIS_ERROR( false, "eval_d3NdXi3 not implemented for this interpolation function" );
        }

        //------------------------------------------------------------------------------

        template< Geometry_Type G, Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< G, T, N, B >::get_param_coords( Matrix< DDRMat >& aXiHat ) const
        {
            MORIS_ERROR( false, "get_param_coords not implemented for this interpolation function" );
        }

        //------------------------------------------------------------------------------

        template< Geometry_Type G, Interpolation_Type T, uint N, uint B >
        Interpolation_Order
        Interpolation_Function< G, T, N, B >::get_interpolation_order() const
        {
            MORIS_ERROR( false, "get_interpolation_order implemented for this interpolation function" );
            return Interpolation_Order::UNDEFINED;
        }

        //------------------------------------------------------------------------------

        template< Geometry_Type G, Interpolation_Type T, uint N, uint B >
        uint
        Interpolation_Function< G, T, N, B >::get_number_of_param_dimensions() const
        {
            MORIS_ERROR( false, "get_number_of_param_dimensions - not implemented for this interpolation function" );
            return 0;
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_HPP_ */
