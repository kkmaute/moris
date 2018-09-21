/*
 * cl_FEM_Matrix.hpp
 *
 *  Created on: Jul 13, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_MATRIX_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_MATRIX_HPP_

#include <string>
#include <utility>
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"   //LINALG/src
#include "linalg_typedefs.hpp"
#include "op_times.hpp" //LINALG/src
#include "op_plus.hpp"  //LINALG/src
#include "op_minus.hpp" //LINALG/src
#include "fn_trans.hpp" //LINALG/src
#include "fn_det.hpp"   //LINALG/src
#include "fn_inv.hpp"   //LINALG/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        // forward declaration needed for interpolation class
        class Interpolator;

//------------------------------------------------------------------------------

        class Interpolation_Matrix
        {
            //! pointer to interpolaition function object
            Interpolator * mInterpolator = nullptr;

            //! space derivative flag set during construction
            const uint mSpaceFlag;

            //! time derivative flag set during construction
            const uint mTimeFlag;

            //! coefficient flag set during construction
            //const uint mCoeffFlag;

            //! matrix that contains data
            Matrix< DDRMat > mData;

            //! pointer to function that evaluates Matrix
            void
            ( *mEvaluate )
                    ( Interpolator         * aInterpolator,
                      Interpolation_Matrix * aMatrix,
                      const Matrix< DDRMat >    & aPoint );


//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            /**
             *  constructor
             */
            
            Interpolation_Matrix(
                    const uint   & aSpaceFlag,
                    const uint   & aTimeFlag,
                    const uint   & aNumberOfRows,
                    const uint   & aNumberOfCols );
            
            Interpolation_Matrix(
                    Interpolator * aInterpolator,
                    const uint   & aSpaceFlag,
                    const uint   & aTimeFlag,
                    const uint   & aNumberOfRows,
                    const uint   & aNumberOfCols );

//             Interpolation_Matrix(
//                     const uint & aSpaceFlag,
//                     const uint & aTimeFlag,
//                     const uint & aCoeffFlag,
//                     const uint & aNumberOfRows,
//                     const uint & aNumberOfCols ) :
//                         mSpaceFlag( aSpaceFlag ),
//                         mTimeFlag( aTimeFlag ),
//                         mCoeffFlag( aCoeffFlag )
//             {
//             mData.set_size(  aNumberOfRows, aNumberOfCols );
// 
//             //this->assign_evaluation_function();
//             }


//------------------------------------------------------------------------------

            /**
             * alternative constructor using moris::mat
             */
            
                Interpolation_Matrix(
                    const uint        & aSpaceFlag,
                    const uint        & aTimeFlag,
                    const Matrix< DDRMat > & aData );

//------------------------------------------------------------------------------

            /**
             * default destructor
             */
            ~Interpolation_Matrix(){};

//------------------------------------------------------------------------------

            /**
             * sets size of member matrix
             */
            void
            set_size(
                    const uint & aNumberOfRows,
                    const uint & aNumberOfCols )
            {
                mData.set_size( aNumberOfRows, aNumberOfCols );
            }

//------------------------------------------------------------------------------

            /**
             * sets size of member matrix
             */
            void
            set_size(
                    const uint & aNumberOfRows,
                    const uint & aNumberOfCols,
                    const real & aValue )
            {
                mData.set_size( aNumberOfRows, aNumberOfCols, aValue );
            }
//------------------------------------------------------------------------------

            /**
             *  expose data object
             */
            DDRMat &
            matrix_data() //-> decltype( mData.matrix_data() )
            {
                return mData.matrix_data();
            }

//------------------------------------------------------------------------------

            /**
             *  expose data object ( const version )
             */
            const DDRMat &
            matrix_data() const // -> decltype( mData.matrix_data() )
            {
                return mData.matrix_data();
            }
//------------------------------------------------------------------------------

            /**
             *  expose matrix object
             */
            Matrix< DDRMat > &
            matrix()
            {
                return mData;
            }

//------------------------------------------------------------------------------


            /**
             *  expose matrix object ( const version )
             */
            const Matrix< DDRMat > &
            matrix() const
            {
                return mData;
            }

//------------------------------------------------------------------------------

            /**
             * allows access to an entry in a vector
             *
             * @param[ in ] aI  index in vector
             */
            auto
            operator()( const uint & aI )
                -> decltype( mData( aI ) )
            {
                return mData( aI );
            }

//------------------------------------------------------------------------------

            /**
             * allows access to an entry in a vector ( const version )
             *
             * @param[ in ] aI  index in vector
             */
            auto
            operator()( const uint & aI ) const
                -> decltype( mData( aI ) )
            {
                return mData( aI );
            }

//------------------------------------------------------------------------------

            /**
             * allows access to an entry in a matrix
             *
             * @param[ in ] aI  row index
             *
             * @param[ in ] aJ  column index
             *
             */
            auto
            operator()( const uint & aI, const uint & aJ )
                -> decltype( mData( aI, aJ ) )
            {
                return mData( aI, aJ );
            }

//------------------------------------------------------------------------------

            /**
             * allows access to an entry in a matrix ( const version )
             *
             * @param[ in ] aI  row index
             *
             * @param[ in ] aJ  column index
             *
             */
            auto
            operator()( const uint & aI, const uint & aJ ) const
                -> decltype( mData( aI, aJ ) )
            {
                return mData( aI, aJ );
            }

//------------------------------------------------------------------------------

            /**
             * returns the length of a vector
             */
            auto
            length() const -> decltype( mData.length() )
            {
                return mData.length();
            }

//------------------------------------------------------------------------------

            /**
             * returns the number of rows of a matrix
             */
            auto
            n_rows() const -> decltype ( mData.n_rows() )
            {
                return mData.n_rows();
            }

//------------------------------------------------------------------------------
            /**
             * returns the number of cols of a matrix
             */
            auto
            n_cols() const -> decltype ( mData.n_cols() )
            {
                return mData.n_cols();
            }

//------------------------------------------------------------------------------

            /**
             * calls the print routine of the moris::Mat
             */
//            void
//            print( const std::string & aVarName = std::string() )
//            {
//                mData.print( aVarName );
//            }

//------------------------------------------------------------------------------

            /**
             * returns the sum of the data ( needed for testing of N unity )
             */
            real
            sum() const
            {

                real aSum = 0.0;
                uint tLength = mData.length();

                for( uint k=0; k<tLength; ++k )
                {
                    aSum += mData( k );
                }

                return aSum;
            }

//------------------------------------------------------------------------------

            /**
             * returns the space derivative flag
             */
            auto get_space_flag() const
                -> decltype ( mSpaceFlag )
            {
                return mSpaceFlag;
            }

//------------------------------------------------------------------------------

            /**
             * returns the time derivative flag
             */
            auto get_time_flag() const
                -> decltype ( mTimeFlag )
            {
                return mTimeFlag;
            }

//------------------------------------------------------------------------------

            /**
             * returns the flag for the coefficients
             */
           /* auto get_coeff_flag() const
                -> decltype ( mCoeffFlag )
            {
                return mCoeffFlag;
            } */

//------------------------------------------------------------------------------

            /**
             * evaluates matrix with respect to linked function
             */
           /* void
            evaluate(
                    const Matrix< DDRMat >           & aPoint )
            {
                // call linked function
                ( *mEvaluate )( aFunction, this, aPoint );
            } */

//------------------------------------------------------------------------------

            /**
             * returns a pointer to the linked interpolation function
             */
            auto
            get_interpolator() -> decltype( mInterpolator )
            {
                return mInterpolator;
            }

//------------------------------------------------------------------------------

            /**
             * called by the creator from the interpolator
             */
            void
            assign_interpolator_and_function( Interpolator * aInterpolator );

//------------------------------------------------------------------------------
            /**
             * evaluates the matrix at given point
             */
            void
            compute( const Matrix< DDRMat > & aPoint );

//------------------------------------------------------------------------------

            /**
             * evaluates the matrix at given integration point
             */
            void
            compute( const uint & aPoint );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
//  free operators
//------------------------------------------------------------------------------


        auto
        operator*(  Interpolation_Matrix & aA,
                    Matrix< DDRMat >     & aB )
            ->  decltype( aA.matrix_data() * aB.matrix_data() )

        {
            return aA.matrix_data() * aB.matrix_data();
        }

//------------------------------------------------------------------------------


        auto
        operator*( Matrix< DDRMat >      & aA,
                   Interpolation_Matrix  & aB )
            ->  decltype( aA.matrix_data() * aB.matrix_data() )

        {
            return aA.matrix_data() * aB.matrix_data();
        }

//------------------------------------------------------------------------------

        auto
        operator*(  const Interpolation_Matrix  & aA,
                    const Interpolation_Matrix  & aB )
        ->  decltype( aA.matrix_data() * aB.matrix_data() )

        {
            return aA.matrix_data() * aB.matrix_data();
        }

//------------------------------------------------------------------------------


        auto
        operator*(  const Interpolation_Matrix  & aA,
                    const Interpolation_Matrix  * aB )
        ->  decltype( aA.matrix_data() * aB->matrix_data() )

        {
            return aA.matrix_data() * aB->matrix_data();
        }

//------------------------------------------------------------------------------

        auto
        operator*(  const Interpolation_Matrix  * aA,
                    const Matrix< DDRMat >      & aB )
        ->  decltype( aA->matrix_data() * aB.matrix_data() )

        {
            return aA->matrix_data() * aB.matrix_data();
        }

//------------------------------------------------------------------------------

        auto
        operator*(  const Matrix< DDRMat >  & aA,
                    const Interpolation_Matrix  * aB )
        ->  decltype( aA.matrix_data() * aB->matrix_data() )

        {
            return aA.matrix_data() * aB->matrix_data() ;
        }

//------------------------------------------------------------------------------

        /**
         * calculates the determinant of a matrix
         * @param[ in ] aA   matrix to process
         */
        auto
        det( Interpolation_Matrix & aA ) -> decltype( det( aA.matrix() ) )
        {
            return det( aA.matrix() );
        }


//------------------------------------------------------------------------------

        /**
         * inverts a matrix
         * @param[ in ] aA   matrix to process
         */
        auto
        inv( Interpolation_Matrix & aA ) -> decltype( inv( aA.matrix() ) )
        {
            return inv( aA.matrix() ) ;
        }

//------------------------------------------------------------------------------

        /**
         * transposes a matrix
         *
         * warning: pointers to functions are not copied
         * @param[ in ] aA   matrix to process
         */
        Interpolation_Matrix
        trans( Interpolation_Matrix & aA )
        {
            return Interpolation_Matrix(
                    aA.get_space_flag(),
                    aA.get_time_flag(),
                    trans( aA.matrix() ) );
        }


//------------------------------------------------------------------------------

        /**
         * transposes a matrix
         * @param[ in ] aA   matrix to process
         */
        Interpolation_Matrix
        trans( Interpolation_Matrix * aA )
        {
            // dereference pointer
            return trans( * aA );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_INTERPOLATION_MATRIX_HPP_ */
