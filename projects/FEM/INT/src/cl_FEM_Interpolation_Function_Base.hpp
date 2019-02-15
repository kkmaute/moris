/*
 * cl_FEM_Interpolation_Function.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_BASE_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_BASE_HPP_

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_Matrix.hpp"   //LINALG/src

namespace moris
{
    namespace fem
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
            virtual Matrix< DDRMat > eval_N( const Matrix< DDRMat > & aXi ) const = 0;

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
            virtual Matrix< DDRMat > eval_dNdXi( const Matrix< DDRMat > & aXi ) const = 0;

//------------------------------------------------------------------------------

            /**
             * calculates the second derivative of the shape function
             * in parameter space
             *
             * @param[ out ] ad2NdXi2 ( <number of dimensions> x <number of nodes> )
             *
             * @param[ in ] aXi    point where function is evaluated
             *                     ( <number of dimensions>  x 1 )
             *
             */
            virtual Matrix< DDRMat > eval_d2NdXi2 ( const Matrix< DDRMat > & aXi ) const = 0;

//------------------------------------------------------------------------------

            /**
             * returns a matrix containing the parameter coordinates
             * < number of dimensions * number of basis >
             */
            virtual Matrix< DDRMat > get_param_coords() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the number of basis for this shape function
             */
            virtual uint get_number_of_bases() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the number of dimensions for this shape function
             */
            virtual uint get_number_of_dimensions() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the interpolation order
             */
            virtual mtk::Interpolation_Order get_interpolation_order() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the interpolation type
             */
            virtual Interpolation_Type get_interpolation_type() const = 0;

//------------------------------------------------------------------------------

        };
//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_BASE_HPP_ */
