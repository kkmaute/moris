
#include "fn_inv.hpp" //LNA/src
#include "fn_det.hpp" //LNA/src
#include "op_times.hpp" //LNA/src
#include "op_equal_equal.hpp" //LNA/src
#include "cl_FEM_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Element.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Interpolator::Interpolator(
                Element                     * aElement,
                const uint                  & aNumberOfFields,
                const Interpolation_Rule    & aFieldInterpolationRule,
                const Interpolation_Rule    & aGeometryInterpolationRule,
                const Integration_Rule      & aIntegrationRule )
                    : mNumberOfFields( aNumberOfFields )

        {
            // create interpolation functions
            if( aFieldInterpolationRule.has_two_rules() )
            {
                // create space function
                mSpaceInterpolation
                    = aFieldInterpolationRule.create_space_interpolation_function();

                // create time function
                mTimeInterpolation
                    = aFieldInterpolationRule.create_time_interpolation_function();

                // set matrix creator
                mMatrixCreator = mSpaceInterpolation;
            }
            else
            {
                mSpaceTimeInterpolation
                    = aFieldInterpolationRule.create_space_time_interpolation_function();

                // set matrix creator
                mMatrixCreator = mSpaceTimeInterpolation;
            }

            // create interpolator
            mGeometryInterpolator = new Geometry_Interpolator(
                    aElement,
                    aGeometryInterpolationRule );

            // get node coordinates
            mNodeCoords = aElement->get_node_coords();

            // number of nodes expected by geometry interpolator
            auto tNumberOfGeometryBasis
                = mGeometryInterpolator->get_number_of_basis();

            uint tNumberOfNodes = mNodeCoords.n_rows();
            uint tNumberOfDimensions = mNodeCoords.n_cols();

            // test if node coords vector needs to be chopped to linear
            if( tNumberOfNodes != tNumberOfGeometryBasis )
            {
                if ( mGeometryInterpolator->get_interpolation_order() == mtk::Interpolation_Order::LINEAR )
                {
                    mNodeCoords.resize( tNumberOfGeometryBasis, tNumberOfDimensions );
                    mIsoparametricFlag = false;
                }
                else
                {
                    MORIS_ERROR( false,
                            "Geometry interpolation must be either linear or of identical order as Field interpolation.");
                }
            }
            else if ( mGeometryInterpolator->get_interpolation_type()
                    == aFieldInterpolationRule.get_type_in_space() )
            {
                mIsoparametricFlag = true;
            }
            else
            {
                mIsoparametricFlag = false;
            }

            if( ! mIsoparametricFlag )
            {
                mGdNdXi = mGeometryInterpolator->create_matrix_pointer( 1, 0, 1 );
            }
            else
            {
                mGdNdXi = mMatrixCreator->create_matrix_pointer( 1, 1, 0, 1 );
            }

            // create integrator
            mIntegrator = new Integrator( aIntegrationRule );

            // store integration points in memory
            mIntegrationPoints = mIntegrator->get_points();

            // store integration weights in memory
            mIntegrationWeights = mIntegrator->get_weights();


            // set dimension of last point
            mLastPointJt.set_size(
                    tNumberOfDimensions,
                    1,
                    1e12); // <- put an value in there which is never used
                           //    by an integration point

        }
//------------------------------------------------------------------------------

        Interpolator::~Interpolator()
        {
            if( mSpaceInterpolation != NULL )
            {
                delete mSpaceInterpolation;
            }
            if( mTimeInterpolation != NULL )
            {
                delete mTimeInterpolation;
            }
            if( mSpaceTimeInterpolation != NULL )
            {
                delete mSpaceTimeInterpolation;
            }

            delete mGeometryInterpolator;

            if ( mIntegrator != NULL )
            {
                delete mIntegrator;
            }

            if( mGdNdXi != NULL )
            {
                delete mGdNdXi;
            }

        }

//------------------------------------------------------------------------------

        uint
        Interpolator::get_number_of_dofs()
        {
            // this needs to be changed if interpolation in space and time
            return mMatrixCreator->get_number_of_basis();
        }

//------------------------------------------------------------------------------

        uint
        Interpolator::get_number_of_integration_points()
        {
            // test if integrator was initialized
            if ( mIntegrator != NULL )
            {
                return mIntegrator->get_number_of_points();
            }
            else
            {
                return 0;
            }
        }

//------------------------------------------------------------------------------

        Interpolation_Matrix
        Interpolator::create_matrix( const uint & aDerivativeInSpace,
                                     const uint & aDerivativeInTime,
                                     const uint & aCoeffsType ) const
        {
            // pass through to member function
            return this->mMatrixCreator->create_matrix(
                    mNumberOfFields,
                    aDerivativeInSpace,
                    aDerivativeInTime,
                    aCoeffsType );
        }

//------------------------------------------------------------------------------

        // This needs to be split into several subroutines.
        // Not sure about the many switches though.
        void
        Interpolator::evaluate_matrix(
                Interpolation_Matrix & aMatrix,
                const uint           & aPoint )
        {
            // pass integration coordinate to other function
            this->evaluate_matrix(
                    aMatrix,
                    mIntegrationPoints.cols( aPoint, aPoint ) );
        }

//------------------------------------------------------------------------------

        real
        Interpolator::get_integration_weight( const uint & aPoint )
        {
            return mIntegrationWeights( aPoint );
        }

//------------------------------------------------------------------------------

        // This needs to be split into several subroutines.
        // Not sure about the many switches though.
        void
        Interpolator::evaluate_matrix(
                Interpolation_Matrix & aMatrix,
                const Mat< real >    & aPoint )
        {
            switch ( aMatrix.get_space_flag() )
            {
                case( 0 ) :
                {
                    switch ( aMatrix.get_time_flag() )
                    {
                        case( 0 ) :
                        {
                            switch ( aMatrix.get_coeff_flag() )
                            {
                                case( 1 ) :
                                {
                                    mMatrixCreator->eval_N( aMatrix, aPoint );
                                    break;
                                }
                                default :
                                {
                                    MORIS_ERROR( false, "evaluate_matrix: coeffs not available");
                                    break;
                                }
                            }
                            break;
                        }
                        default :
                        {
                            MORIS_ERROR( false, "evaluate_matrix: time derivative not available");
                            break;
                        }
                    }
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "evaluate_matrix: space derivative not available");
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Interpolator::eval_dNdx(
                Interpolation_Matrix & aMatrix,
                const Mat< real >    & aPoint )
        {
            mMatrixCreator->eval_dNdXi( aMatrix, aPoint );

            // test if element is isoparametric
            if ( mIsoparametricFlag )
            {
                // calculate geometry jacobian
                mGeometryInterpolator->eval_jacobian( mJt, aMatrix, mNodeCoords );
            }
            else
            {
                // calculate geometry derivative
                mGeometryInterpolator->eval_dNdXi( *mGdNdXi, aPoint );

                // calculate geometry Jacobian
                mGeometryInterpolator->eval_jacobian( mJt, *mGdNdXi, mNodeCoords );
            }

            // transform output matrix
            aMatrix.data() = inv( mJt ) * aMatrix.data();

            // remember point
            mLastPointJt = aPoint;
        }

//------------------------------------------------------------------------------

        real
        Interpolator::get_det_J( const uint & aPoint )
        {
            // pass to other function
            return this->get_det_J( mIntegrationPoints.cols( aPoint, aPoint ) );
        }

//------------------------------------------------------------------------------

        real
        Interpolator::get_det_J( const Mat< real > & aPoint )
        {
            // test if Jacobi matrix is up to date
            if ( ( aPoint == mLastPointJt ).min() == 0 )
            {
                // calculate derivative
                if ( mIsoparametricFlag )
                {
                    mMatrixCreator->eval_dNdXi( *mGdNdXi, aPoint );
                }
                else
                {
                    mGeometryInterpolator->eval_dNdXi( *mGdNdXi, aPoint );
                }

                // calculate geometry Jacobian
                mGeometryInterpolator->eval_jacobian( mJt, *mGdNdXi, mNodeCoords );

                // remember point
                mLastPointJt = aPoint;
            }
            return det( mJt );
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

