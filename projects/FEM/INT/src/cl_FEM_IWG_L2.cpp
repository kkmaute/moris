#include "cl_FEM_IWG_L2.hpp"

#include "op_times.hpp" //LINALG/src
#include "fn_norm.hpp"  //LINALG/src
#include "fn_trans.hpp" //LINALG/src
#include "fn_dot.hpp"   //LINALG/src
#include "fn_print.hpp" //LINALG/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

    IWG_L2::IWG_L2( const real aAlpha )
    {
        // set the residual dof type
        mResidualDofType = MSI::Dof_Type::L2;

        // set the active dof types
        mActiveDofTypes = { { MSI::Dof_Type::L2} };

        // set alpha
        this->set_alpha( aAlpha );
    }

//------------------------------------------------------------------------------

    void IWG_L2::set_alpha( const real aAlpha  )
    {

        mAlpha = aAlpha;

        if(  aAlpha == 0.0 )
        {
            mComputeFunction
                = & IWG_L2::compute_jacobian_and_residual_without_alpha;

            mComputeJacFunction
                = & IWG_L2::compute_jacobian_without_alpha;

            mComputeResFunction
                = & IWG_L2::compute_residual_without_alpha;
        }
        else
        {
            mComputeFunction
                = & IWG_L2::compute_jacobian_and_residual_with_alpha;

            mComputeJacFunction
                = & IWG_L2::compute_jacobian_with_alpha;

            mComputeResFunction
                = & IWG_L2::compute_residual_with_alpha;
        }
    }

//------------------------------------------------------------------------------

//        void IWG_L2::create_matrices( Field_Interpolator * aInterpolator )
//        {
//            // copy pointer to interpolator class
//            //mInterpolator = aInterpolator;
//
//            // create N-Matrix
//            //mN = aInterpolator->create_matrix( 0, 0 );
//
//            // create B-Matrix
//            //mB = aInterpolator->create_matrix( 1, 0 );
//        }

//------------------------------------------------------------------------------

//        void
//        IWG_L2::delete_matrices()
//        {
//            delete mN;
//            delete mB;
//        }

//------------------------------------------------------------------------------
//
//        void
//        IWG_L2::compute_jacobian_and_residual(
//                Matrix< DDRMat >       & aJacobian,
//                Matrix< DDRMat >       & aResidual,
//                const Matrix< DDRMat > & aNodalDOF,
//                const Matrix< DDRMat > & aNodalWeakBC,
//                const uint             & aPointIndex )
//        {
//            ( this->*mComputeFunction )(
//                    aJacobian,
//                    aResidual,
//                    aNodalDOF,
//                    aNodalWeakBC,
//                    aPointIndex );
//        }

        void IWG_L2::compute_jacobian_and_residual(       Matrix< DDRMat > & aJacobian,
                                                          Matrix< DDRMat > & aResidual,
                                                    const Matrix< DDRMat > & aNodalDOF,
                                                    const Matrix< DDRMat > & aNodalWeakBC )
        {
            ( this->*mComputeFunction )( aJacobian,
                                         aResidual,
                                         aNodalDOF,
                                         aNodalWeakBC );
        }

//------------------------------------------------------------------------------
//
//        real
//        IWG_L2::compute_integration_error(
//                const Matrix< DDRMat >                    & aNodalDOF,
//                real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
//                const uint                                & aPointIndex )
//        {
//            mN->compute( aPointIndex );
//
//            Matrix< DDRMat > tCoords = mN->matrix_data() * mInterpolator->get_node_coords();
//
//            // get shape function
//            Matrix< DDRMat > tPhiHat = mN->matrix_data() * aNodalDOF.matrix_data();
//
//            return std::pow( tPhiHat( 0 ) - aFunction( tCoords ), 2 );
//        }

//------------------------------------------------------------------------------
//
//        void
//        IWG_L2::compute_jacobian_and_residual_without_alpha(
//                Matrix< DDRMat >       & aJacobian,
//                Matrix< DDRMat >       & aResidual,
//                const Matrix< DDRMat > & aNodalDOF,
//                const Matrix< DDRMat > & aNodalWeakBC,
//                const uint             & aPointIndex )
//        {
//
//            // get shape function
//            mN->compute( aPointIndex );
//
//            // calculate Jacobian
//            aJacobian = trans( mN ) * mN;
//            aResidual = aJacobian * ( aNodalDOF - aNodalWeakBC );
//        }


        void IWG_L2::compute_jacobian_and_residual_without_alpha(       Matrix< DDRMat > & aJacobian,
                                                                        Matrix< DDRMat > & aResidual,
                                                                  const Matrix< DDRMat > & aNodalDOF,
                                                                  const Matrix< DDRMat > & aNodalWeakBC )
        {

            // get shape function
            Matrix< DDRMat > tN = mInterpolator->N();

            // compute Jacobian
            aJacobian = trans( tN ) * tN;

            // compute residual
            aResidual = aJacobian * ( aNodalDOF - aNodalWeakBC );
        }

//------------------------------------------------------------------------------
//
//        void
//        IWG_L2::compute_jacobian_and_residual_with_alpha(
//                Matrix< DDRMat >       & aJacobian,
//                Matrix< DDRMat >       & aResidual,
//                const Matrix< DDRMat > & aNodalDOF,
//                const Matrix< DDRMat > & aNodalWeakBC,
//                const uint             & aPointIndex )
//        {
//
//            // get shape function
//            mN->compute( aPointIndex );
//
//            // compute derivative
//            mB->compute( aPointIndex );
//
//            // calculate Jacobian
//            aJacobian = trans( mN ) * mN + mAlpha * ( trans( mB ) * mB );
//
//            aResidual = aJacobian * ( aNodalDOF - aNodalWeakBC );
//        }


        void IWG_L2::compute_jacobian_and_residual_with_alpha(       Matrix< DDRMat > & aJacobian,
                                                                     Matrix< DDRMat > & aResidual,
                                                               const Matrix< DDRMat > & aNodalDOF,
                                                               const Matrix< DDRMat > & aNodalWeakBC )
        {

            // get shape function
            Matrix< DDRMat > tN = mInterpolator->N();

            // compute derivative
            Matrix< DDRMat > tB = mInterpolator->Bx();

            // compute Jacobian
            aJacobian = trans( tN ) * tN + mAlpha * ( trans( tB ) * tB );

            // compute residual
            aResidual = aJacobian * ( aNodalDOF - aNodalWeakBC );
        }

//------------------------------------------------------------------------------
//
//        real
//        IWG_L2::interpolate_scalar_at_point(
//                                    const Matrix< DDRMat > & aNodalWeakBC,
//                                    const uint             & aPointIndex )
//        {
//            // get shape function
//            mN->compute( aPointIndex );
//
//            // return interpolation
//            return dot( mN->matrix() , aNodalWeakBC );
//        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                       Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            ( this->*mComputeJacFunction )( aJacobians,
                                            aFieldInterpolators );
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_without_alpha( Cell< Matrix< DDRMat > >    & aJacobians,
                                                     Cell< Field_Interpolator* > & aFieldInterpolators )
        {
                   // set field interpolator
                   Field_Interpolator* tFI = aFieldInterpolators( 0 );

                   // compute Jacobian
                   aJacobians( 0 ) = trans( tFI->N() ) * tFI->N();
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_with_alpha( Cell< Matrix< DDRMat > >    & aJacobians,
                                                  Cell< Field_Interpolator* > & aFieldInterpolators )
               {
                   // set field interpolator
                   Field_Interpolator* tInterpolator = aFieldInterpolators( 0 );

                   // get shape function
                   Matrix< DDRMat > tN = tInterpolator->N();

                   // compute derivative
                   Matrix< DDRMat > tB = tInterpolator->Bx();

                   // compute Jacobian
                   aJacobians( 0 ) = trans( tN ) * tN + mAlpha * ( trans( tB ) * tB );
               }

//------------------------------------------------------------------------------

        void IWG_L2::compute_residual( Matrix< DDRMat >            & aResidual,
                                       Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            ( this->*mComputeResFunction )( aResidual,
                                            aFieldInterpolators );
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_residual_without_alpha( Matrix< DDRMat >            & aResidual,
                                                     Cell< Field_Interpolator* > & aFieldInterpolators )
               {
                   // set field interpolator
                   Field_Interpolator* tFieldInterpolator = aFieldInterpolators( 0 );

                   // get shape function
                   Matrix< DDRMat > tN = tFieldInterpolator->N();

                   // get field value
                   Matrix< DDRMat > tField = tFieldInterpolator->val();

                   //FIXME: enforced weak BCs
                   Matrix< DDRMat > tNodalWeakBCs( tN.n_cols(), 1, 1.0);

                   // compute Jacobian
                   aResidual = trans( tN ) * ( tField - tN * tNodalWeakBCs );
               }

//------------------------------------------------------------------------------

        void IWG_L2::compute_residual_with_alpha( Matrix< DDRMat >            & aResidual,
                                                  Cell< Field_Interpolator* > & aFieldInterpolators )
               {
                   // set field interpolator
                   Field_Interpolator* mFieldInterpolator = aFieldInterpolators( 0 );

                   // get shape function
                   Matrix< DDRMat > tN = mFieldInterpolator->N();

                   // compute derivative
                   Matrix< DDRMat > tB = mFieldInterpolator->Bx();

                   // get shape function
                   Matrix< DDRMat > tField = mFieldInterpolator->val();

                   // compute derivative
                   Matrix< DDRMat > tFieldGradx = mFieldInterpolator->gradx( 1 );

                   //FIXME: enforced weak BCs
                   Matrix< DDRMat > tNodalWeakBCs( tN.n_cols(), 1, 1.0);
                   // compute Jacobian
                   aResidual = trans( tN ) * ( tField - tN * tNodalWeakBCs )
                             + mAlpha * trans( tB ) * ( tFieldGradx - tB * tNodalWeakBCs );
               }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
