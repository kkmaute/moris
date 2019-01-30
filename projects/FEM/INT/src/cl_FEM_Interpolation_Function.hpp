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
        template< Interpolation_Type T, uint N, uint B  >
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
            void
            eval_N(       Interpolation_Matrix  & aN,
                    const Matrix< DDRMat > & aXi  ) const;

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
            void
            eval_dNdXi(        Interpolation_Matrix  & adNdXi,
                         const Matrix< DDRMat >     & aXi ) const;

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
            void
            eval_d2NdXi2 (
                           Interpolation_Matrix & ad2NdXi2,
                    const Matrix< DDRMat >     & aXi ) const;

//------------------------------------------------------------------------------

            /**
             * returns a matrix containing the parameter coordinates
             * < number of dimensions * number of basis >
             */
            void
            get_param_coords( Matrix< DDRMat > & aXihat ) const;

//------------------------------------------------------------------------------

            /**
             * returns the number of basis for this shape function
             */
            uint
            get_number_of_basis() const
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
             * returns the interpolation order
             */
            mtk::Interpolation_Order
            get_interpolation_order() const;

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

            /**
             * creates am interpolation matrix
             *
             * @param[ in ]  aDerivativeInSpace, 0, 1 or 2
             * @param[ in ]  aDerivativeInTime   0, 1 or 2
             */
            Interpolation_Matrix
            create_matrix(
                    const uint & aNumberOfFields,
                    const uint & aDerivativeInSpace,
                    const uint & aDerivativeInTime ) const;

//------------------------------------------------------------------------------

            /**
             * creates a pointer to a new interpolation matrix
             *
             * @param[ in ]  aDerivativeInSpace, 0, 1 or 2
             * @param[ in ]  aDerivativeInTime   0, 1 or 2
             */
            Interpolation_Matrix *
            create_matrix_pointer(
                    const uint & aNumberOfFields,
                    const uint & aDerivativeInSpace,
                    const uint & aDerivativeInTime ) const;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            get_matrix_size(
                          uint & aNumberOfRows,
                          uint & aNumberOfCols,
                    const uint & aDerivativeInSpace,
                    const uint & aDerivativeInTime ) const;
        };

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< T, N, B>::eval_N(       Interpolation_Matrix  & aN,
                                                  const Matrix< DDRMat > & aXi  ) const
        {
            MORIS_ERROR( false,
                "eval_N not implemented for this interpolation function" );
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< T, N, B>::eval_dNdXi(
                      Interpolation_Matrix  & adNdXi,
                const Matrix< DDRMat > & aXi  ) const
        {
            MORIS_ERROR( false,
                "eval_dNdXi not implemented for this interpolation function" );
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< T, N, B>::eval_d2NdXi2 (
                       Interpolation_Matrix & ad2NdXi2,
                const Matrix< DDRMat > & aXi ) const
        {
            MORIS_ERROR( false,
                "eval_d2NdXi2 not implemented for this interpolation function" );
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< T, N, B>::get_param_coords( Matrix< DDRMat > & aXihat ) const
        {
            MORIS_ERROR( false,
                "get_param_coords not implemented for this interpolation function" );
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        mtk::Interpolation_Order
        Interpolation_Function< T, N, B>::get_interpolation_order() const
        {
            MORIS_ERROR( false,
                    "get_interpolation_order implemented for this interpolation function" );
            return mtk::Interpolation_Order::UNDEFINED;
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        Interpolation_Matrix
        Interpolation_Function< T, N, B >::create_matrix(
                            const uint & aNumberOfFields,
                            const uint & aDerivativeInSpace,
                            const uint & aDerivativeInTime ) const
        {
            uint tNumberOfRows;
            uint tNumberOfCols;

            // determine number of rows and cols
            this->get_matrix_size(
                    tNumberOfRows,
                    tNumberOfCols,
                    aDerivativeInSpace,
                    aDerivativeInTime );

            // return new matrix
            return Interpolation_Matrix(
                    aDerivativeInSpace,
                    aDerivativeInTime,
                    tNumberOfRows,
                    tNumberOfCols );
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        Interpolation_Matrix *
        Interpolation_Function< T, N, B >::create_matrix_pointer(
                const uint & aNumberOfFields,
                const uint & aDerivativeInSpace,
                const uint & aDerivativeInTime ) const
        {
            uint tNumberOfRows;
            uint tNumberOfCols;

            // determine number of rows and cols
            this->get_matrix_size(
                    tNumberOfRows,
                    tNumberOfCols,
                    aDerivativeInSpace,
                    aDerivativeInTime );

            // return new matrix
            return new Interpolation_Matrix(
                    aDerivativeInSpace,
                    aDerivativeInTime,
                    tNumberOfRows,
                    tNumberOfCols );
        }

//------------------------------------------------------------------------------

        template< Interpolation_Type T, uint N, uint B >
        void
        Interpolation_Function< T, N, B >::get_matrix_size(
                uint & aNumberOfRows,
                uint & aNumberOfCols,
                const uint & aDerivativeInSpace,
                const uint & aDerivativeInTime ) const
       {
            // determine number of rows
            switch( aDerivativeInSpace )
            {
                case( 0 ) :
                {
                    aNumberOfRows = 1;
                    break;
                }
                case( 1 ) :
                {
                    aNumberOfRows = N;
                    break;
                }
                case( 2 ) :
                {
                    uint tSecondDeriv[ 3 ] = { 1, 3, 6 };
                    aNumberOfRows = tSecondDeriv[ N-1 ];
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown number of dimensions.");
                    aNumberOfRows = 0;
                    break;
                }
            }

            aNumberOfCols = B;
            // determine number of columns
            /*switch( aCoeffsSwitch )
            {
                case( 0 ) :    // evaluated property
                {
                    aNumberOfCols = 1;
                    break;
                }
                case( 1 ) :    // coefficients
                {
                    aNumberOfCols = B;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown aCoeffsSwitch");
                    aNumberOfCols = 0;
                    break;
                }
            }*/
       }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_HPP_ */
