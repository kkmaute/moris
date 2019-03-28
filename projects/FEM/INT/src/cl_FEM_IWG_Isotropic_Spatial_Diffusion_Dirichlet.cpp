
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

            // set the active dof type
            mActiveDofTypes = { { MSI::Dof_Type::TEMP } };
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_residual
            ( Matrix< DDRMat >            & aResidual,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            // compute the residual r_T
            aResidual = trans( tTemp->N() ) * ( tTemp->val() - tTemp->N() * mNodalWeakBCs );
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian
            ( Cell< Matrix< DDRMat > >    & aJacobians,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = trans( tTemp->N() ) * tTemp->N();
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual
            ( Cell< Matrix< DDRMat > >    & aJacobians,
              Matrix< DDRMat >            & aResidual,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            // compute the residual r_T
            aResidual = trans( tTemp->N() ) * ( tTemp->val() - tTemp->N() * mNodalWeakBCs );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = trans( tTemp->N() ) * tTemp->N();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
