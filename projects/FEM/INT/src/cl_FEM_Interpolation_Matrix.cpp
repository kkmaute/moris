//#include "cl_FEM_Interpolation_Matrix.hpp"
//#include "cl_FEM_Interpolator.hpp"
//
//namespace moris
//{
//    namespace fem
//    {
////------------------------------------------------------------------------------
//
//        Interpolation_Matrix::Interpolation_Matrix(
//                const uint        & aSpaceFlag,
//                const uint        & aTimeFlag,
//                const uint        & aNumberOfRows,
//                const uint        & aNumberOfCols ) :
//                                        mSpaceFlag( aSpaceFlag ),
//                                        mTimeFlag( aTimeFlag ),
//                                        mData( aNumberOfRows, aNumberOfCols )
//        {
//        }
//
////------------------------------------------------------------------------------
//
//      Interpolation_Matrix::Interpolation_Matrix( const uint             & aSpaceFlag,
//                                                  const uint             & aTimeFlag,
//                                                  const Matrix< DDRMat > & aData ) : mSpaceFlag( aSpaceFlag ),
//                                                                                     mTimeFlag( aTimeFlag ),
//                                                                                     mData( aData )
//        {
//        }
//
////------------------------------------------------------------------------------
//
//        void Interpolation_Matrix::assign_interpolator_and_function( Interpolator * aInterpolator )
//        {
//            // set pointer to interpolator
//            mInterpolator = aInterpolator;
//
//            // set pointer to evaluation function
//            switch ( mTimeFlag )
//            {
//                case( 0 ) :
//                {
//                    switch( mSpaceFlag )
//                    {
////                        case ( 0 ) :
////                        {
////                            mEvaluate = &interpolator_eval_N;
////                            break;
////                        }
////                        case( 1 ) :
////                        {
////                            mEvaluate = &interpolator_eval_dNdx;
////                            break;
////                        }
//                        default :
//                        {
//                            MORIS_ERROR( false, "interpolation function not implemented");
//                        }
//                    }
//
//                    break;
//                }
//                default :
//                {
//                    MORIS_ERROR( false, "interpolation function not implemented");
//                }
//            }
//        }
//
////------------------------------------------------------------------------------
//
//        /**
//         * evaluates the matrix at given point
//         */
//        void Interpolation_Matrix::compute( const Matrix< DDRMat > & aPoint )
//        {
//            this->mEvaluate( mInterpolator, this, aPoint );
//        }
//
////------------------------------------------------------------------------------
//
//        /**
//         * evaluates the matrix at given integration point
//         */
//        void Interpolation_Matrix::compute( const uint & aPoint )
//        {
//            this->mEvaluate( mInterpolator, this, mInterpolator->get_point( aPoint ) );
//
//            auto tPoint = mInterpolator->get_point( aPoint );
//        }
//
////------------------------------------------------------------------------------
//
//    } /* namespace fem */
//} /* namespace moris */
