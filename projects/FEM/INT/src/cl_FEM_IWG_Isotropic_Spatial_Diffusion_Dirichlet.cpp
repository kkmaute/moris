
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Isotropic_Spatial_Diffusion_Dirichlet::IWG_Isotropic_Spatial_Diffusion_Dirichlet()
        {
            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::TEMP };

            // set the active dof types
            mActiveDofTypes = { { MSI::Dof_Type::TEMP } };

            // set the active mp type
            mActivePropertyTypes = { fem::Property_Type::TEMP_DIRICHLET };

            // FIXME set a penalty
            mGamma = 1.0;

            //FIXME forced diffusion parameter
            //      forced dimensions for 3D
            eye( mSpaceDim, mSpaceDim, mKappa );
            mKappa = 1.0 * mKappa;
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_residual
            ( Matrix< DDRMat >                   & aResidual,
              moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            // interpolated TBar from nodal values
            Matrix< DDRMat > tTBar = tTemp->N() * mNodalWeakBCs;

            // compute the residual r_T
            aResidual = - trans( tTemp->N() ) * dot( mKappa * tTemp->gradx( 1 ), mNormal )
                        + trans( mKappa * tTemp->Bx() ) * mNormal * ( tTemp->val()( 0 ) - tTBar( 0 ) )
                        + mGamma * trans( tTemp->N() ) * ( tTemp->val()( 0 ) - tTBar( 0 ) );
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian
            ( moris::Cell< Matrix< DDRMat > >    & aJacobians,
              moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = - trans( tTemp->N() ) * trans( mNormal ) * mKappa * tTemp->Bx()
                              + trans( mKappa * tTemp->Bx() ) * mNormal * tTemp->N()
                              + mGamma * trans( tTemp->N() ) * tTemp->N();

        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual
            ( moris::Cell< Matrix< DDRMat > >    & aJacobians,
              Matrix< DDRMat >                   & aResidual,
              moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            // interpolated TBar from nodal values
            Matrix< DDRMat > tTBar = tTemp->N() * mNodalWeakBCs;

            // compute the residual r_T
            aResidual = - trans( tTemp->N() ) * dot( mKappa * tTemp->gradx( 1 ), mNormal )
                        + trans( mKappa * tTemp->Bx() ) * mNormal * ( tTemp->val()( 0 ) - tTBar( 0 ) )
                        + mGamma * trans( tTemp->N() ) * ( tTemp->val()( 0 ) - tTBar( 0 ) );

            // set the jacobian size
                       aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = - trans( tTemp->N() ) * trans( mNormal ) * mKappa * tTemp->Bx()
                              + trans( mKappa * tTemp->Bx() ) * mNormal * tTemp->N()
                              + mGamma * trans( tTemp->N() ) * tTemp->N();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
