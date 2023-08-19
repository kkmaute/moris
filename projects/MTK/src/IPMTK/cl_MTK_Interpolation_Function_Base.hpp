/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Base.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_BASE_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_BASE_HPP_

#include "cl_MTK_Enums.hpp"    //MTK/src
#include "cl_MTK_Enums.hpp"    //MTK/src
#include "cl_Matrix.hpp"       //LINALG/src

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------

        /**
         * shape function base class
         */
        class Interpolation_Function_Base
        {

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Interpolation_Function_Base(){};

            //------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            virtual ~Interpolation_Function_Base(){};

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
            virtual void eval_N(
                    const Matrix< DDRMat >& aXi,
                    Matrix< DDRMat >&       aNXi ) const = 0;

            //------------------------------------------------------------------------------

            /**
             * calculates the first derivative of the shape function
             * in parameter space
             *
             * @param[ in ] adNdXi ( <number of dimensions> x <number of nodes> )
             *
             * @param[ in ] aXi    point where function is evaluated
             *                     ( <number of dimensions>  x 1 )
             *
             */
            virtual void eval_dNdXi(
                    const Matrix< DDRMat >& aXi,
                    Matrix< DDRMat >&       adNdXi ) const = 0;

            //------------------------------------------------------------------------------

            /**
             * calculates the second derivative of the shape function
             * in parameter space
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
            virtual void eval_d2NdXi2(
                    const Matrix< DDRMat >& aXi,
                    Matrix< DDRMat >&       ad2NdXi2 ) const = 0;

            //------------------------------------------------------------------------------

            /**
             * calculates the third derivatives of the shape function
             * in parameter space
             *
             * @param[ in ] aXi    point where function is evaluated
             *                     ( <number of dimensions>  x 1 )
             *
             * @param[ out ] ad3NdXi3 ( <number of dimensions> x <number of nodes> )
             */
            virtual void eval_d3NdXi3(
                    const Matrix< DDRMat >& aXi,
                    Matrix< DDRMat >&       ad3NdXi3 ) const = 0;

            //------------------------------------------------------------------------------

            /**
             * returns a matrix containing the parameter coordinates
             * < number of dimensions * number of basis >
             */
            virtual void get_param_coords( Matrix< DDRMat >& aXiHat ) const = 0;

            //------------------------------------------------------------------------------

            /**
             * returns the number of bases for this shape function
             */
            virtual uint get_number_of_bases() const = 0;

            //------------------------------------------------------------------------------

            /**
             * returns the number of dimensions for this shape function
             */
            virtual uint get_number_of_dimensions() const = 0;

            //------------------------------------------------------------------------------

            /**
             * returns the number of parametric dimensions for this shape function
             */
            virtual uint get_number_of_param_dimensions() const = 0;

            //------------------------------------------------------------------------------

            /**
             * returns the interpolation order
             */
            virtual Interpolation_Order get_interpolation_order() const = 0;

            //------------------------------------------------------------------------------

            /**
             * returns the interpolation type
             */
            virtual Interpolation_Type get_interpolation_type() const = 0;

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------

    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_BASE_HPP_ */

