
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
            eye( 3, 3, mKappa );
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
            Field_Interpolator* T  = aFieldInterpolators( 0 );

            // FIXME: set normal value for 3D
            Matrix < DDRMat > aNormal( 3, 1, 1.0 );

            // compute the residual r_T
            aResidual = trans( T->N() ) * dot( mKappa * T->gradx( 1 ), aNormal );
        }

//------------------------------------------------------------------------------

        void
		IWG_Isotropic_Spatial_Diffusion_Sideset::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* T  = aFieldInterpolators( 0 );

            // FIXME: set normal value for 3D
            Matrix < DDRMat > aNormal( 3, 1, 1.0 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = trans( T->N() ) *  trans( trans( T->Bx() ) * mKappa * aNormal ) ;
        }

//------------------------------------------------------------------------------

        void
		IWG_Isotropic_Spatial_Diffusion_Sideset::compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                                 Matrix< DDRMat >            & aResidual,
                                                                                 Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* T  = aFieldInterpolators( 0 );

            // FIXME: set normal value for 3D
            Matrix < DDRMat > aNormal( 3, 1, 1.0 );

            // compute the residual r_T
            aResidual = trans( T->N() ) * dot( mKappa * T->gradx( 1 ), aNormal );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = trans( T->N() ) *  trans( trans( T->Bx() ) * mKappa * aNormal ) ;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
