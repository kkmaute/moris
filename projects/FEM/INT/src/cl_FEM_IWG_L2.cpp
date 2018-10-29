#include "cl_FEM_IWG_L2.hpp"
#include "fn_norm.hpp"
#include "op_times.hpp"                     //LINALG/src
#include "fn_trans.hpp"                     //LINALG/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

    IWG_L2::IWG_L2( const real aAlpha ) : mAlpha( aAlpha )
    {
        if(  aAlpha == 0.0 )
        {
            mComputeFunction
                = & IWG_L2::compute_jacobian_and_residual_without_alpha;
        }
        else
        {
            mComputeFunction
                = & IWG_L2::compute_jacobian_and_residual_with_alpha;
        }
    }

//------------------------------------------------------------------------------

        void
        IWG_L2::create_matrices( Interpolator * aInterpolator )
        {
            // copy pointer to interpolator class
            mInterpolator = aInterpolator;

            // create N-Matrix
            mN = aInterpolator->create_matrix( 0, 0 );

            // create B-Matrix
            mB = aInterpolator->create_matrix( 1, 0 );

        }

//------------------------------------------------------------------------------

        void
        IWG_L2::delete_matrices()
        {
            delete mN;
            delete mB;
        }

//------------------------------------------------------------------------------

        void
        IWG_L2::compute_jacobian_and_residual(
                Matrix< DDRMat >       & aJacobian,
                Matrix< DDRMat >       & aResidual,
                const Matrix< DDRMat > & aNodalDOF,
                const Matrix< DDRMat > & aNodalWeakBC,
                const uint             & aPointIndex )
        {
            ( this->*mComputeFunction )(
                    aJacobian,
                    aResidual,
                    aNodalDOF,
                    aNodalWeakBC,
                    aPointIndex );
        }
//------------------------------------------------------------------------------

        real
        IWG_L2::compute_integration_error(
                const Matrix< DDRMat >                    & aNodalDOF,
                real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
                const uint                                & aPointIndex )
        {
            mN->compute( aPointIndex );

            Matrix< DDRMat > tCoords = mN * mInterpolator->get_node_coords();

            // get shape function
            Matrix< DDRMat > tPhiHat = mN->matrix_data() * aNodalDOF.matrix_data();

            return std::pow(
                    tPhiHat( 0 )
                    - aFunction( tCoords ), 2 );
        }

//------------------------------------------------------------------------------

        void
        IWG_L2::compute_jacobian_and_residual_without_alpha(
                Matrix< DDRMat >       & aJacobian,
                Matrix< DDRMat >       & aResidual,
                const Matrix< DDRMat > & aNodalDOF,
                const Matrix< DDRMat > & aNodalWeakBC,
                const uint        & aPointIndex )
        {

            // get shape function
            mN->compute( aPointIndex );

            // calculate Jacobian
            aJacobian = trans( mN ) * mN;

            aResidual = aJacobian * ( aNodalDOF - aNodalWeakBC );
        }

//------------------------------------------------------------------------------

        void
        IWG_L2::compute_jacobian_and_residual_with_alpha(
                Matrix< DDRMat >       & aJacobian,
                Matrix< DDRMat >       & aResidual,
                const Matrix< DDRMat > & aNodalDOF,
                const Matrix< DDRMat > & aNodalWeakBC,
                const uint        & aPointIndex )
        {

            // get shape function
            mN->compute( aPointIndex );

            // compute derivative
            mB->compute( aPointIndex );

            // calculate Jacobian
            aJacobian = trans( mN ) * mN + mAlpha * ( trans( mB ) * mB );

            aResidual = aJacobian * ( aNodalDOF - aNodalWeakBC );
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
