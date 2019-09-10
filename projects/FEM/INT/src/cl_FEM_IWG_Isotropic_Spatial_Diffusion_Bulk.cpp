
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Isotropic_Spatial_Diffusion_Bulk::IWG_Isotropic_Spatial_Diffusion_Bulk(){}

//------------------------------------------------------------------------------

        void
        IWG_Isotropic_Spatial_Diffusion_Bulk::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            //FIXME heat load enforced
            Matrix< DDRMat > tQ( 1, 1, 0.0 );

            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // compute conductivity matrix
            Matrix< DDRMat > K;
            eye( mSpaceDim, mSpaceDim, K );
            K = mMasterProp( 0 )->val()( 0 ) * K;

            // set residual size
            this->set_residual( aResidual );

            // compute the residual r_T
            aResidual( 0 ) = trans( mMasterFI( 0 )->Bx() ) * K * mMasterFI( 0 )->gradx( 1 )
                           - trans( mMasterFI( 0 )->N() ) * tQ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // compute conductivity matrix
            Matrix< DDRMat > K;
            eye( mSpaceDim, mSpaceDim, K );
            K = mMasterProp( 0 )->val()( 0 ) * K;

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian j_T_T
            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->Bx() ) * K * mMasterFI( 0 )->Bx();
        }

//------------------------------------------------------------------------------

        void
        IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                             moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            //FIXME heat load enforced
            Matrix< DDRMat > tQ( 1, 1, 0.0 );

            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // compute conductivity matrix
            Matrix< DDRMat > K;
            eye( mSpaceDim, mSpaceDim, K );
            K = mMasterProp( 0 )->val()( 0 ) * K;

            // set the jacobian size
            this->set_residual( aResidual );

            // compute the residual r_T
            aResidual( 0 ) = trans( mMasterFI( 0 )->Bx() ) * K * mMasterFI( 0 )->gradx( 1 )
                           + trans( mMasterFI( 0 )->N() ) * tQ;

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian j_T_T
            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->Bx() ) * K * mMasterFI( 0 )->Bx();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
