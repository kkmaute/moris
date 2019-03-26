
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Sideset.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Isotropic_Spatial_Diffusion_Sideset::IWG_Isotropic_Spatial_Diffusion_Sideset()
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
        IWG_Isotropic_Spatial_Diffusion_Sideset::compute_residual( Matrix< DDRMat >            & aResidual,
                                                                   Cell< Field_Interpolator* > & aFieldInterpolators )
        {

            // set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            // compute the residual r_T
            aResidual = - trans( tTemp->N() ) * dot( mKappa * tTemp->gradx( 1 ), mNormal );
        }

//------------------------------------------------------------------------------

        void
        IWG_Isotropic_Spatial_Diffusion_Sideset::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                   Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = - trans( tTemp->N() ) *  trans( trans( tTemp->Bx() ) * mKappa * mNormal ) ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Isotropic_Spatial_Diffusion_Sideset::compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                                Matrix< DDRMat >            & aResidual,
                                                                                Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            // compute the residual r_T
            aResidual = trans( tTemp->N() ) * dot( mKappa * tTemp->gradx( 1 ), mNormal );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = - trans( tTemp->N() ) *  trans( trans( tTemp->Bx() ) * mKappa * mNormal ) ;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
