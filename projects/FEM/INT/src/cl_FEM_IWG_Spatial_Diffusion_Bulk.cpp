
#include "cl_FEM_IWG_Spatial_Diffusion_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Spatial_Diffusion_Bulk::IWG_Spatial_Diffusion_Bulk()
        {
            //FIXME set diffusion parameter
            mKappa = 1.0;

            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::TEMP };

            // set the active dof type
            mActiveDofTypes = { { MSI::Dof_Type::TEMP } };
        }

//------------------------------------------------------------------------------

        void IWG_Spatial_Diffusion_Bulk::compute_residual( Matrix< DDRMat >            & aResidual,
                                                           Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* T  = aFieldInterpolators( 0 );

            // compute the residual r_T
            aResidual = mKappa * trans( T->Bx() ) * T->gradx( 1 );
        }

//------------------------------------------------------------------------------

        void IWG_Spatial_Diffusion_Bulk::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                                           Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* T  = aFieldInterpolators( 0 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = mKappa * trans( T->Bx() ) * T->Bx();
        }

//------------------------------------------------------------------------------

        void IWG_Spatial_Diffusion_Bulk::compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                        Matrix< DDRMat >            & aResidual,
                                                                        Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* T  = aFieldInterpolators( 0 );

            // compute the residual r_T
            aResidual = mKappa * trans( T->Bx() ) * T->gradx( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = mKappa * trans( T->Bx() ) * T->Bx();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
