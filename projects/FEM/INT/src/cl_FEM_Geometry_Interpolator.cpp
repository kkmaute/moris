
#include "cl_FEM_Geometry_Interpolator.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Geometry_Interpolator::Geometry_Interpolator( const Interpolation_Rule & aInterpolationRule )
        {
            // create member pointer to space interpolation function
            mSpaceInterpolation = aInterpolationRule.create_space_interpolation_function();

            // create member pointer to time interpolation function
            mTimeInterpolation  = aInterpolationRule.create_time_interpolation_function();

            // set default xHat
            mXHat.set_size( mSpaceInterpolation->get_number_of_bases(),
                            mSpaceInterpolation->get_number_of_dimensions(), 0.0);

            // set default tHat
            mTHat.set_size( mTimeInterpolation->get_number_of_bases(),
                            mTimeInterpolation->get_number_of_dimensions(), 0.0);

            // set pointers for second derivative depending on space and time dimensions
            this->set_function_pointers();
        }

//------------------------------------------------------------------------------

        Geometry_Interpolator::~Geometry_Interpolator()
        {
            // delete interpolation functions
            if( mSpaceInterpolation != NULL )
             {
                 delete mSpaceInterpolation;
             }
             if( mTimeInterpolation != NULL )
             {
                 delete mTimeInterpolation;
             }
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::set_coeff( const Matrix< DDRMat > & aXHat,
                                               const Matrix< DDRMat > & aTHat )
        {
            //check the space coefficients input size
            MORIS_ASSERT( ( ( aXHat.n_cols() == mXHat.n_cols() ) && ( aXHat.n_rows() == mXHat.n_rows() )),
                          " Geometry_Interpolator::set_coeff - Wrong input size (aXHat). ");

            // set the space coefficients
            mXHat = aXHat;

            //check the time coefficients input size
            MORIS_ASSERT( ( ( aTHat.n_cols() == mTHat.n_cols() ) && ( aTHat.n_rows() == mTHat.n_rows() )),
                           " Geometry_Interpolator::set_coeff - Wrong input size (aTHat). ");

            // set the time coefficients
            mTHat = aTHat;
        }

//------------------------------------------------------------------------------

        Matrix < DDRMat > Geometry_Interpolator::NXi( const Matrix< DDRMat > & aXi ) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> tN = mSpaceInterpolation->eval_N( aXi );
            return tN;
         }

//------------------------------------------------------------------------------

         Matrix < DDRMat > Geometry_Interpolator::NTau( const Matrix< DDRMat > & aTau ) const
         {
             // pass data through interpolation function
             Matrix <DDRMat> tN = mTimeInterpolation->eval_N( aTau );
             return tN;
         }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> tdNdXi = mSpaceInterpolation->eval_dNdXi( aXi );
            return tdNdXi;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::dNdTau( const Matrix< DDRMat > & aTau) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> tdNdTau = mTimeInterpolation->eval_dNdXi( aTau );
            return tdNdTau;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> td2NdXi2 = mSpaceInterpolation->eval_d2NdXi2( aXi );
            return td2NdXi2;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::d2NdTau2( const Matrix< DDRMat > & aTau ) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> td2NdTau2 = mTimeInterpolation->eval_d2NdXi2( aTau );
            return td2NdTau2;
        }
//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::space_jacobian( const Matrix< DDRMat > & adNdXi ) const
        {
            Matrix< DDRMat > tJt = adNdXi * mXHat ;
            return tJt;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::time_jacobian( const Matrix< DDRMat > & adNdTau ) const
        {
            Matrix< DDRMat > tJt = adNdTau * mTHat ;
            return tJt;
        }

//------------------------------------------------------------------------------

        real Geometry_Interpolator::det_J( const Matrix< DDRMat > & aParamPoint )
        {
            // get space dimension
            uint tNSpaceDim = mSpaceInterpolation->get_number_of_dimensions();

            // get tXi and tTau
            Matrix< DDRMat > tXi( tNSpaceDim, 1, 0.0 );
            for ( moris::uint Ik = 0; Ik < tNSpaceDim; Ik++ )
            {
                // set input values
                tXi( Ik ) = aParamPoint( Ik );
            }
            //tXi( { 0, tNSpaceDim-1 }, { 0, 0 } ) = aParamPoint( { 0, tNSpaceDim-1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, aParamPoint( tNSpaceDim ) );

            // get the space jacobian
            Matrix< DDRMat > tdNSpacedXi = this->dNdXi( tXi );
            Matrix< DDRMat > tSpaceJt    = this->space_jacobian( tdNSpacedXi );

            // get the time Jacobian
            Matrix< DDRMat > tdNTimedTau = this->dNdTau( tTau );
            Matrix< DDRMat > tTimeJt     = this->time_jacobian( tdNTimedTau );

            // compute the determinant of the space time Jacobian
            return det( tSpaceJt ) * det( tTimeJt );

        }
//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::valx( const Matrix< DDRMat > & aXi )
        {
//            // evaluate the space time shape functions at Xi, Tau
//            Matrix< DDRMat > tN = this->NXi( aXi );

            //evaluate the field
            return this->NXi( aXi ) * mXHat ;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::valt( const Matrix< DDRMat > & aTau )
        {
//            // evaluate the space time shape functions at Xi, Tau
//            Matrix< DDRMat > tN = this->NTau( aTau );

            //evaluate the field
            return this->NTau( aTau ) * mTHat ;
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::space_jacobian_and_matrices_for_second_derivatives(
                      Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & adNdXi,
                const Matrix< DDRMat > & ad2NdXi2) const
        {
            // evaluate transposed of geometry Jacobian
            aJt = this->space_jacobian( adNdXi );

            // call calculator for second derivatives
            this->mSecondDerivativeMatricesSpace( aJt,
                                                  aKt,
                                                  aLt,
                                                  ad2NdXi2,
                                                  mXHat);
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::time_jacobian_and_matrices_for_second_derivatives(
                      Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & adNdTau,
                const Matrix< DDRMat > & ad2NdTau2) const
        {
            // evaluate transposed of geometry Jacobian
            aJt = this->time_jacobian( adNdTau );

            // call calculator for second derivatives
            this->mSecondDerivativeMatricesTime( aJt,
                                                 aKt,
                                                 aLt,
                                                 ad2NdTau2,
                                                 mTHat);
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_1d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat)
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 1, 1 );
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_2d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 3, 3 );
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
            aLt( 1, 0 ) = std::pow( aJt( 1, 0 ), 2 );
            aLt( 2, 0 ) = aJt( 0 , 0 ) * aJt( 1 , 0 );

            aLt( 0, 1 ) = std::pow( aJt( 0, 1 ), 2 );
            aLt( 1, 1 ) = std::pow( aJt( 1, 1 ), 2 );
            aLt( 2, 1 ) = aJt( 0 , 1 ) * aJt( 1, 1 );

            aLt( 0, 2 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 1 );
            aLt( 1, 2 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 1 );
            aLt( 2, 2 ) = aJt( 0, 0 )* aJt( 1, 1 ) +  aJt( 0, 1 ) * aJt( 1, 0 );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_3d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 6, 6 );
            for( uint j=0; j<3; ++j )
            {
                aLt( 0, j ) = std::pow( aJt( 0, j ), 2 );
                aLt( 1, j ) = std::pow( aJt( 1, j ), 2 );
                aLt( 2, j ) = std::pow( aJt( 2, j ), 2 );
                aLt( 3, j ) = aJt( 1 , j ) * aJt( 2 , j );
                aLt( 4, j ) = aJt( 0 , j ) * aJt( 2 , j );
                aLt( 5, j ) = aJt( 0 , j ) * aJt( 1 , j );
            }

            aLt( 0, 3 ) = 2.0 * aJt( 0, 1 ) * aJt( 0, 2 );
            aLt( 1, 3 ) = 2.0 * aJt( 1, 1 ) * aJt( 1, 2 );
            aLt( 2, 3 ) = 2.0 * aJt( 2, 1 ) * aJt( 2, 2 );
            aLt( 3, 3 ) = aJt( 1, 1 ) * aJt( 2, 2 ) + aJt( 2, 1 ) * aJt( 1, 2 );
            aLt( 4, 3 ) = aJt( 0, 1 ) * aJt( 2, 2 ) + aJt( 2, 1 ) * aJt( 0, 2 );
            aLt( 5, 3 ) = aJt( 0, 1 ) * aJt( 1, 2 ) + aJt( 1, 1 ) * aJt( 0, 2 );

            aLt( 0, 4 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 2 );
            aLt( 1, 4 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 2 );
            aLt( 2, 4 ) = 2.0 * aJt( 2, 0 ) * aJt( 2, 2 );
            aLt( 3, 4 ) = aJt( 1, 0 ) * aJt( 2, 2 ) + aJt( 2, 0 ) * aJt( 1, 2 );
            aLt( 4, 4 ) = aJt( 0, 0 ) * aJt( 2, 2 ) + aJt( 2, 0 ) * aJt( 0, 2 );
            aLt( 5, 4 ) = aJt( 0, 0 ) * aJt( 1, 2 ) + aJt( 1, 0 ) * aJt( 0, 2 );

            aLt( 0, 5 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 1 );
            aLt( 1, 5 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 1 );
            aLt( 2, 5 ) = 2.0 * aJt( 2, 0 ) * aJt( 2, 1 );
            aLt( 3, 5 ) = aJt( 1, 0 ) * aJt( 2, 1 ) + aJt( 2, 0 ) * aJt( 1, 1 );
            aLt( 4, 5 ) = aJt( 0, 0 ) * aJt( 2, 1 ) + aJt( 2, 0 ) * aJt( 0, 1 );
            aLt( 5, 5 ) = aJt( 0, 0 ) * aJt( 1, 1 ) + aJt( 1, 0 ) * aJt( 0, 1 );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::set_function_pointers()
        {
            // get number of dimensions and set pointer to function
            // for second space derivative
            switch ( mSpaceInterpolation->get_number_of_dimensions() )
            {
                case( 1 ) :
                {
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_1d;
                    break;
                }
                case( 2 ) :
                {
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_2d;
                    break;
                }
                case( 3 ) :
                {
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_3d;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, " Geometry_Interpolator::set_function_pointers - unknown number of dimensions. " );
                    break;
                }
            }

            // get number of dimensions and set pointer to function
            // for second time derivative
            switch ( mTimeInterpolation->get_number_of_dimensions() )
            {
                case( 1 ) :
                {
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_1d;
                    break;
                }
                case( 2 ) :
                {
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_2d;
                    break;
                }
                case( 3 ) :
                {
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_3d;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, " Geometry_Interpolator::set_function_pointers - unknown number of dimensions. " );
                    break;
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
