//
//
//#include "cl_Matrix.hpp" //LNA/src
//#include "linalg_typedefs.hpp" //LNA/src
//#include "fn_inv.hpp" //LNA/src
//#include "fn_det.hpp" //LNA/src
//#include "fn_trans.hpp"
//#include "fn_norm.hpp"
//#include "fn_reshape.hpp"
//#include "op_times.hpp" //LNA/src
//#include "op_equal_equal.hpp" //LNA/src
//#include "cl_FEM_Interpolator.hpp" //FEM/INT/src
//#include "cl_FEM_Element.hpp" //FEM/INT/src
//
//
//
//namespace moris
//{
//    namespace fem
//    {
//
////------------------------------------------------------------------------------
//
//        Interpolator::Interpolator( const Matrix< DDRMat >      & aNodalCoords,
//                                    const uint                  & aNumberOfFields,
//                                    const Interpolation_Rule    & aFieldInterpolationRule,
//                                    const Interpolation_Rule    & aGeometryInterpolationRule,
//                                    const Integration_Rule      & aIntegrationRule ) : mNumberOfFields( aNumberOfFields )
//        {
//            // create interpolation functions
////            if( aFieldInterpolationRule.has_two_rules() )
////            {
//                // create space function
//                mSpaceInterpolation = aFieldInterpolationRule.create_space_interpolation_function();
//
//                // create time function
//                mTimeInterpolation = aFieldInterpolationRule.create_time_interpolation_function();
//
//                mMatrixCreator.resize( 2, nullptr );
//
//                // set matrix creator
//                mMatrixCreator( 0 ) = mSpaceInterpolation;
//                mMatrixCreator( 1 ) = mTimeInterpolation;
////            }
////            else
////            {
////                mSpaceInterpolation = aFieldInterpolationRule.create_space_time_interpolation_function();
////
////                mMatrixCreator.resize( 1, nullptr );
////
////                // set matrix creator
////                mMatrixCreator( 0 ) = mSpaceInterpolation;
////            }
//
//            // create interpolator
//            mGeometryInterpolator = new Geometry_Interpolator( aGeometryInterpolationRule );
//            // get node coordinates
//            mNodeCoords = aNodalCoords;
//
//            // number of nodes expected by geometry interpola//            //get number of dimensions in time to set tSpaceTimedNdt size
//            //            uint tNDimTime = mTimeInterpolation->get_number_of_dimensions();tor
//            auto tNumberOfGeometryBasis = mGeometryInterpolator->get_number_of_bases();
//
//            uint tNumberOfNodes = mNodeCoords.n_rows();
//            uint tNumberOfDimensions = mNodeCoords.n_cols();
//
//            // test if node coords vector needs to be chopped to linear
//            if( tNumberOfNodes != tNumberOfGeometryBasis )
//            {
//                if ( mGeometryInterpolator->get_interpolation_order() == mtk::Interpolation_Order::LINEAR )
//                {
//                    mNodeCoords.resize( tNumberOfGeometryBasis, tNumberOfDimensions );
//                    mIsoparametricFlag = false;
//                }
//                else
//                {
//                    MORIS_ERROR( false, "Geometry interpolation must be either linear or of identical order as Field interpolation.");
//                }
//            }
//            else if ( mGeometryInterpolator->get_interpolation_type() == aFieldInterpolationRule.get_type_in_space() )
//            {
//                mIsoparametricFlag = true;
//            }
//            else
//            {
//                mIsoparametricFlag = false;
//            }
//
////            if( ! mIsoparametricFlag )
////            {
////                mGN     = mGeometryInterpolator->create_matrix_pointer( 0, 0 );
////                mGdNdXi = mGeometryInterpolator->create_matrix_pointer( 1, 0 );
////            }
////            else
////            {
////                mGN     = mMatrixCreator( 0 )->create_matrix_pointer( 1, 0, 0 );
////                mGdNdXi = mMatrixCreator( 0 )->create_matrix_pointer( 1, 1, 0 );
////            }
//
//            // create integrator
//            mIntegrator = new Integrator( aIntegrationRule );
//
//            // store integration points in memory
//            mIntegrationPoints = mIntegrator->get_points();
//
//            // store integration weights in memory
//            mIntegrationWeights = mIntegrator->get_weights();
//
//
//            // set dimension of last point
//            mLastPointJt.set_size( tNumberOfDimensions,
//                                   1,
//                                   1e12); // <- put an value in there which is never used
//                                         //    by an integration point
//        }
////------------------------------------------------------------------------------
//
//        Interpolator::~Interpolator()
//        {
//            if( mSpaceInterpolation != NULL )
//            {
//                delete mSpaceInterpolation;
//            }
//            if( mTimeInterpolation != NULL )
//            {
//                delete mTimeInterpolation;
//            }
////            if( mSpaceTimeInterpolation != NULL )
////            {
////                delete mSpaceTimeInterpolation;
////            }
//
//            delete mGeometryInterpolator;
//
//            if ( mIntegrator != NULL )
//            {
//                delete mIntegrator;
//            }
//
////            if( mGN != NULL )
////            {
////                delete mGN;
////            }
////
////            if( mGdNdXi != NULL )
////            {
////                delete mGdNdXi;
////            }
//        }
//
////------------------------------------------------------------------------------
//
//        uint Interpolator::get_number_of_dofs()
//        {
//            // this needs to be changed if interpolation in space and time
//            return mMatrixCreator( 0 )->get_number_of_bases();
//        }
//
////------------------------------------------------------------------------------
//
//        uint Interpolator::get_number_of_integration_points()
//        {
//            // test if integrator was initialized
//            if ( mIntegrator != NULL )
//            {
//                return mIntegrator->get_number_of_points();
//            }
//            else
//            {
//                return 0;
//            }
//        }
//
////------------------------------------------------------------------------------
//
////        Interpolation_Matrix * Interpolator::create_matrix( const uint & aDerivativeInSpace,
////                                                            const uint & aDerivativeInTime )
////        {
////            // pass through to member function
////            Interpolation_Matrix * aMatrix = this->mMatrixCreator( 0 )
////                                                 ->create_matrix_pointer( mNumberOfFields,
////                                                                         aDerivativeInSpace,
////                                                                         aDerivativeInTime);
////
////            // set interpolator and function
////            aMatrix->assign_interpolator_and_function( this );
////
////            return aMatrix;
////        }
//
////------------------------------------------------------------------------------
//
//        real Interpolator::get_integration_weight( const uint & aPoint )
//        {
//            return mIntegrationWeights( aPoint );
//        }
//
////------------------------------------------------------------------------------
//
////        void Interpolator::eval_N(       Interpolation_Matrix & aMatrix,
////                                   const Matrix< DDRMat >     & aPoint)
////        {
////                mMatrixCreator( 0 )->eval_N( aMatrix, aPoint );
////        }
//
//        Matrix< DDRMat > Interpolator::eval_N( const Matrix< DDRMat > & aPoint)
//        {
//            Matrix< DDRMat > tN;
//            return tN = mMatrixCreator( 0 )->eval_N(aPoint );
//        }
//
////        void Interpolator::eval_N(       Interpolation_Matrix & aMatrix,
////                                   const Matrix< DDRMat >     & aPoint,
////                                   const Matrix< DDRMat >     & aTime )
////        {
////            // Get size of ...
////            moris::uint tNumInterpolationMatrices = mMatrixCreator.size();
////
////            if ( tNumInterpolationMatrices == 1 )
////            {
////                mMatrixCreator( 0 )->eval_N( aMatrix, aPoint );
////            }
////            else if ( tNumInterpolationMatrices == 2 )
////            {
////                moris::Cell < Interpolation_Matrix > ListInterpolationMatrices( tNumInterpolationMatrices );
////
////                for ( moris::uint Ik = 0; Ik < tNumInterpolationMatrices; Ik++ )
////                {
////                    if ( Ik == 0 )
////                    {
////                        ListInterpolationMatrices( Ik ).set_size( 1, 4 );
////                        mMatrixCreator( Ik )->eval_N( ListInterpolationMatrices( Ik ), aPoint );
////                    }
////                    if ( Ik == 1 )
////                    {
////                        ListInterpolationMatrices( Ik ).set_size( 1, 2 );
////                        mMatrixCreator( Ik )->eval_N( ListInterpolationMatrices( Ik ), aTime );
////                    }
////                    print( ListInterpolationMatrices( Ik ).matrix(), "tSpaceTimeN");
////                }
////
////                aMatrix =  trans( ListInterpolationMatrices( 0 ) ) * ListInterpolationMatrices( 1 ) ;
////
////                aMatrix.matrix() = reshape( aMatrix.matrix(), (moris::size_t)1,(moris::size_t) 8);
////            }
////            else
////            {
////                MORIS_ERROR( false, "Interpolator::eval_N: more than 2 interpolation matrices");
////            }
////
////        }
//
//        Matrix< DDRMat > Interpolator::eval_N( const Matrix< DDRMat > & aPoint,
//                                               const Matrix< DDRMat > & aTime )
//        {
//             // Get size of ...
//             moris::uint tNumInterpolationMatrices = mMatrixCreator.size();
//
//             if ( tNumInterpolationMatrices == 1 )
//             {
//                  Matrix< DDRMat > tN;
//                  return tN = mMatrixCreator( 0 )->eval_N( aPoint );
//             }
//             else if ( tNumInterpolationMatrices == 2 )
//             {
//                  moris::Cell < Matrix< DDRMat > > ListInterpolationMatrices( tNumInterpolationMatrices );
//
//                  for ( moris::uint Ik = 0; Ik < tNumInterpolationMatrices; Ik++ )
//                  {
//                       if ( Ik == 0 )
//                       {
//                           ListInterpolationMatrices( Ik ) = mMatrixCreator( Ik )->eval_N( aPoint );
//                        }
//                        if ( Ik == 1 )
//                        {
//                            ListInterpolationMatrices( Ik ) = mMatrixCreator( Ik )->eval_N( aTime );
//                        }
//                  }
//
//                  Matrix< DDRMat > aMatrix;
//                  aMatrix =  trans( ListInterpolationMatrices( 0 ) ) * ListInterpolationMatrices( 1 ) ;
//                  aMatrix = reshape( aMatrix, (moris::size_t)1,(moris::size_t) 8);
//                  return aMatrix;
//             }
//             else
//             {
//                 MORIS_ERROR( false, "Interpolator::eval_N: more than 2 interpolation matrices");
//                 Matrix< DDRMat > aEmpty;
//                 return aEmpty;
//              }
//        }
//
////------------------------------------------------------------------------------
//
////        void Interpolator::eval_dNdx(       Interpolation_Matrix & aMatrix,
////                                      const Matrix< DDRMat >     & aPoint )
////        {
////            mMatrixCreator( 0 )->eval_dNdXi( aMatrix, aPoint );
////
////            // test if element is isoparametric
////            if ( mIsoparametricFlag )
////            {
////                // calculate geometry jacobian
////                mGeometryInterpolator->eval_jacobian( mJt, aMatrix, mNodeCoords );
////            }
////            else
////            {
////                // calculate geometry derivative
////                mGeometryInterpolator->eval_dNdXi( *mGdNdXi, aPoint );
////
////                // calculate geometry Jacobian
////                mGeometryInterpolator->eval_jacobian( mJt, *mGdNdXi, mNodeCoords );
////            }
////
////            // transform output matrix
////            aMatrix.matrix_data() = inv( mJt ) * aMatrix.matrix_data();
////
////            // remember point
////            mLastPointJt = aPoint;
////        }
//        Matrix< DDRMat > Interpolator::eval_dNdx( const Matrix< DDRMat > & aPoint )
//        {
//            Matrix < DDRMat > tdNdx;
//            tdNdx = mMatrixCreator( 0 )->eval_dNdXi( aPoint );
//
//            // test if element is isoparametric
//            if ( mIsoparametricFlag )
//            {
//                // calculate geometry jacobian
//            	mJt = mGeometryInterpolator->eval_jacobian( tdNdx );
//            }
//            else
//            {
//                // calculate geometry derivative
//                Matrix < DDRMat > tGdNdXi = mGeometryInterpolator->eval_dNdXi( aPoint );
//
//                // calculate geometry Jacobian
//                mJt = mGeometryInterpolator->eval_jacobian( tGdNdXi );
//            }
//
//            // transform output matrix
//            tdNdx = inv( mJt ) * tdNdx;
//            return tdNdx;
//
//            // remember point
//            mLastPointJt = aPoint;
//        }
//
////------------------------------------------------------------------------------
//
//        real Interpolator::get_det_J( const uint & aPoint )
//        {
//            // pass to other function
//            return this->get_det_J( mIntegrationPoints.get_column( aPoint ) );
//        }
//
////------------------------------------------------------------------------------
//
//        real Interpolator::get_det_J( const Matrix< DDRMat > & aPoint )
//        {
//            // test if Jacobi matrix is up to date
//            //if ( norm( aPoint - mLastPointJt ) == 0.0 )
//            {
//                // calculate derivative
//            	Matrix < DDRMat > tGdNdXi;
//                if ( mIsoparametricFlag )
//                {
//                    tGdNdXi = mMatrixCreator( 0 )->eval_dNdXi( aPoint );
//                }
//                else
//                {
//                    tGdNdXi = mGeometryInterpolator->eval_dNdXi( aPoint );
//                }
//
//                // calculate geometry Jacobian
//                mJt = mGeometryInterpolator->eval_jacobian( tGdNdXi );
//
//                // remember point
//                mLastPointJt = aPoint;
//            }
//            return det( mJt );
//        }
//
////------------------------------------------------------------------------------
//
////        Matrix< DDRMat > Interpolator::eval_geometry_coords( const Matrix< DDRMat >    & aPoint )
////        {
////            //@fixme make sure that this is only calculated once
////            // calculate matrix
////            if ( mIsoparametricFlag )
////            {
////                mMatrixCreator( 0 )->eval_N( *mGN, aPoint );
////            }
////            else
////            {
////                mGeometryInterpolator->eval_N( *mGN, aPoint );
////            }
////
////            return moris::trans( ( * mGN ) * (mNodeCoords) );
////        }
//
//        Matrix< DDRMat > Interpolator::eval_geometry_coords( const Matrix< DDRMat >    & aPoint )
//        {
//            //@fixme make sure that this is only calculated once
//            // calculate matrix
//            Matrix< DDRMat > tGN;
//            if ( mIsoparametricFlag )
//            {
//                tGN = mMatrixCreator( 0 )->eval_N( aPoint );
//            }
//            else
//            {
//                tGN = mGeometryInterpolator->eval_N( aPoint );
//            }
//
//            return moris::trans( ( tGN ) * (mNodeCoords) );
//        }
//
////------------------------------------------------------------------------------
//
//        Matrix< DDRMat > Interpolator::eval_geometry_coords( const uint    & aPoint )
//        {
//            return this->eval_geometry_coords( mIntegrationPoints.get_column( aPoint ) );
//        }
//
////------------------------------------------------------------------------------
//    } /* namespace fem */
//} /* namespace moris */
//
