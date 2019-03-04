/*
 * cl_FEM_Interpolation_Function.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_HPP_

#include "assert.hpp"

#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src
#include "cl_Matrix.hpp"   //LINALG/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        /**
         * shape function templated class
         * T : Type
         * N : Dimension
         * B : Number of Basis
         */
        template< Interpolation_Type T, uint N, uint B >
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
            Matrix< DDRMat > eval_N( const Matrix< DDRMat > & aXi ) const;

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

            Matrix< DDRMat > eval_dNdXi( const Matrix< DDRMat > & aXi ) const;

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
            Matrix< DDRMat > eval_d2NdXi2 ( const Matrix< DDRMat > & aXi ) const;

//------------------------------------------------------------------------------

            /**
             * returns a matrix containing the parameter coordinates
             * < number of dimensions * number of basis >
             */
            Matrix< DDRMat > get_param_coords() const;

//------------------------------------------------------------------------------

            /**
             * returns the number of bases for this shape function
             */
            uint get_number_of_bases() const
            {
                return B;
            }

//------------------------------------------------------------------------------

            /**
             * returns the number of dimensions for this shape function
             */
            uint get_number_of_dimensions() const
            {
                return N;
            }

//------------------------------------------------------------------------------

            /**
             * returns the interpolation order
             */
            mtk::Interpolation_Order get_interpolation_order() const;

//------------------------------------------------------------------------------

            /**
             * returns the interpolation type
             */
            Interpolation_Type get_interpolation_type() const
            {
                return T;
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        Matrix< DDRMat > Interpolation_Function< T, N, B>::eval_N( const Matrix< DDRMat > & aXi ) const
        {
            MORIS_ERROR( false, "eval_N not implemented for this interpolation function" );
            Matrix< DDRMat > aEmpty;
            return aEmpty;
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        Matrix< DDRMat > Interpolation_Function< T, N, B>::eval_dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            MORIS_ERROR( false, "eval_dNdXi not implemented for this interpolation function" );
            Matrix< DDRMat > aEmpty;
            return aEmpty;
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        Matrix< DDRMat > Interpolation_Function< T, N, B>::eval_d2NdXi2 ( const Matrix< DDRMat > & aXi ) const
        {
            MORIS_ERROR( false, "eval_d2NdXi2 not implemented for this interpolation function" );
            Matrix< DDRMat > aEmpty;
            return aEmpty;
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        Matrix< DDRMat > Interpolation_Function< T, N, B>::get_param_coords() const
        {
            MORIS_ERROR( false, "get_param_coords not implemented for this interpolation function" );
            Matrix< DDRMat > aEmpty;
            return aEmpty;
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        mtk::Interpolation_Order Interpolation_Function< T, N, B>::get_interpolation_order() const
        {
            MORIS_ERROR( false, "get_interpolation_order implemented for this interpolation function" );
            return mtk::Interpolation_Order::UNDEFINED;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_HPP_ */
