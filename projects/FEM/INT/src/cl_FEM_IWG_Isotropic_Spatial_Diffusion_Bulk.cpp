
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Isotropic_Spatial_Diffusion_Bulk::IWG_Isotropic_Spatial_Diffusion_Bulk()
        {
            //FIXME forced diffusion parameter
            //      forced dimensions for 3D
            eye( mSpaceDim, mSpaceDim, mKappa );
            mKappa = 1.0 * mKappa;

            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::TEMP };

            // set the active dof type
            mActiveDofTypes = { { MSI::Dof_Type::TEMP } };
        }

//------------------------------------------------------------------------------

        void
        IWG_Isotropic_Spatial_Diffusion_Bulk::compute_residual( Matrix< DDRMat >            & aResidual,
                                                                Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            //fixme heat load enforced
            Matrix< DDRMat > tQ( 1, 1, 0.0 );

            // compute the residual r_T
            aResidual = trans( tTemp->Bx() ) * mKappa * tTemp->gradx( 1 )
                      + trans( tTemp->N() ) * tQ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = trans( tTemp->Bx() ) * mKappa * tTemp->Bx();
        }

//------------------------------------------------------------------------------

        void
        IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                             Matrix< DDRMat >            & aResidual,
                                                                             Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            //fixme heat load enforced
            Matrix< DDRMat > tQ( 1, 1, 0.0 );

            // compute the residual r_T
            aResidual = trans( tTemp->Bx() ) * mKappa * tTemp->gradx( 1 )
                      + trans( tTemp->N() ) * tQ;

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = trans( tTemp->Bx() ) * mKappa * tTemp->Bx();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
