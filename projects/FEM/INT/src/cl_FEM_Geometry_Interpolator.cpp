
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Element.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

//        Geometry_Interpolator::Geometry_Interpolator(
//                Element                   * aElement,
//                const Interpolation_Rule  & aInterpolationRule  )
    Geometry_Interpolator::Geometry_Interpolator( const Interpolation_Rule  & aInterpolationRule  )
        {
            // create a factory
            Interpolation_Function_Factory tFactory;

            // make sure that rule is legal for this element

            //MORIS_ERROR( aElement->get_geometry_type() == aInterpolationRule.get_geometry_type(),
            //         "chosen interpolation rule not allowed for this element" );

            // create member pointer to interpolation function ( actually only space)
            //mInterpolation = aInterpolationRule.create_space_time_interpolation_function();
            mInterpolation = aInterpolationRule.create_space_interpolation_function();   //FIXME

            // set pointers for second derivative depending on space dimensions
            this->set_function_pointers();
        }

//------------------------------------------------------------------------------

        Geometry_Interpolator::~Geometry_Interpolator()
        {
            // delete interpolation function
            delete mInterpolation;
        }

//------------------------------------------------------------------------------

        uint Geometry_Interpolator::get_number_of_basis() const
        {
            return mInterpolation->get_number_of_basis();
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_N(       Interpolation_Matrix  & aN,
                                            const Matrix< DDRMat >      & aXi  ) const
        {
            // pass data through interpolation function
            this->mInterpolation->eval_N( aN, aXi );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_dNdXi(      Interpolation_Matrix  & adNdXi,
                                                const Matrix< DDRMat >     & aXi  ) const
        {
            // pass data through interpolation function
            this->mInterpolation->eval_dNdXi( adNdXi, aXi );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_d2NdXi2(       Interpolation_Matrix  & ad2NdXi2,
                                                  const Matrix< DDRMat >      & aXi  ) const
        {
            // pass data through interpolation function
            this->mInterpolation->eval_d2NdXi2( ad2NdXi2, aXi );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_jacobian(
                  Matrix< DDRMat > 		& aJt,
            const Interpolation_Matrix  & adNdXi,
            const Matrix< DDRMat > 		& aXhat ) const
        {
            aJt = adNdXi.matrix_data() * aXhat.matrix_data() ;
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_jacobian_and_matrices_for_second_derivatives(
                      Matrix< DDRMat > 		& aJt,
                      Matrix< DDRMat > 		& aKt,
                      Matrix< DDRMat > 		& aLt,
                const Interpolation_Matrix  & adNdXi,
                const Interpolation_Matrix  & ad2NdXi2,
                const Matrix< DDRMat >   	& aXhat ) const
        {
            // evaluate transposed of geometry Jacobian
            this->eval_jacobian( aJt, adNdXi, aXhat );

            // call calculator for second derivatives
            this->mSecondDerivativeMatrices(
                    aJt,
                    aKt,
                    aLt,
                    adNdXi,
                    ad2NdXi2,
                    aXhat );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_1d(
                const Matrix< DDRMat > 		& aJt,
                      Matrix< DDRMat > 		& aKt,
                      Matrix< DDRMat > 		& aLt,
                const Interpolation_Matrix 	& adNdXi,
                const Interpolation_Matrix 	& ad2NdXi2,
                const Matrix< DDRMat > 		& aXhat )
        {
            // help matrix K
            aKt = ad2NdXi2.matrix_data() * aXhat.matrix_data();

            // help matrix L
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_2d(
                const Matrix< DDRMat > 		& aJt,
                      Matrix< DDRMat > 		& aKt,
                      Matrix< DDRMat > 		& aLt,
                const Interpolation_Matrix  & adNdXi,
                const Interpolation_Matrix  & ad2NdXi2,
                const Matrix< DDRMat > 		& aXhat )
        {
            // help matrix K
            aKt = ad2NdXi2.matrix_data() * aXhat.matrix_data();

            // help matrix L
            //aLt.resize(3,3);
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
            aLt( 1, 0 ) = std::pow( aJt( 1, 0 ), 2 );
            aLt( 2, 0 ) = aJt( 0 , 0 ) * aJt( 1 , 0 );

            aLt( 0, 1 ) = std::pow( aJt( 0, 1 ), 2 );
            aLt( 1, 1 ) = std::pow( aJt( 1, 1 ), 2 );
            aLt( 2, 0 ) = aJt( 0 , 1 ) * aJt( 1, 1 );

            aLt( 0, 2 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 1 );
            aLt( 1, 2 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 1 );
            aLt( 2, 2 ) = aJt( 0, 0 )* aJt( 1, 1 ) +  aJt( 0, 1 ) * aJt( 1, 0 );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_3d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Interpolation_Matrix  & adNdXi,
                const Interpolation_Matrix  & ad2NdXi2,
                const Matrix< DDRMat > & aXhat )
        {
            // help matrix K
            aKt = ad2NdXi2.matrix_data() * aXhat.matrix_data();

            // help matrix L
            //aLt.resize(5,5);
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

        Interpolation_Matrix * Geometry_Interpolator::create_matrix_pointer(
                const uint & aDerivativeInSpace,
                const uint & aDerivativeInTime ) const
        {
            // pass through to member function
           return this->mInterpolation->create_matrix_pointer(
                    1,
                    aDerivativeInSpace,
                    aDerivativeInTime );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::set_function_pointers()
        {
            // get number of dimensions and set pointer to function
            // for second derivative
            switch ( mInterpolation->get_number_of_dimensions() )
            {
                case( 1 ) :
                {
                    mSecondDerivativeMatrices = this->eval_matrices_for_second_derivative_1d;
                    break;
                }
                case( 2 ) :
                {
                    mSecondDerivativeMatrices = this->eval_matrices_for_second_derivative_2d;
                    break;
                }
                case( 3 ) :
                {
                    mSecondDerivativeMatrices = this->eval_matrices_for_second_derivative_3d;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown number of dimensions" );
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * returns the order of the interpolation
         */
        mtk::Interpolation_Order Geometry_Interpolator::get_interpolation_order() const
        {
            return mInterpolation->get_interpolation_order();
        }

//------------------------------------------------------------------------------

        /**
         * returns the order of the interpolation
         */
        Interpolation_Type Geometry_Interpolator::get_interpolation_type() const
        {
            return mInterpolation->get_interpolation_type();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
