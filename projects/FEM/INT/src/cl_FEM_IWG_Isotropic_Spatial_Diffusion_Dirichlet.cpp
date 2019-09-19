
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
            // FIXME set a penalty
            mGamma = 1.0;
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
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

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * dot( K * mMasterFI( 0 )->gradx( 1 ), mNormal )
                           + trans( K * mMasterFI( 0 )->Bx() ) * mNormal * ( mMasterFI( 0 )->val()( 0 ) - mMasterProp( 1 )->val()( 0 ) )
                           + mGamma * trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->val()( 0 ) - mMasterProp( 1 )->val()( 0 ) );
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
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
            aJacobians( 0 )( 0 ) = - trans( mMasterFI( 0 )->N() ) * trans( mNormal ) * K * mMasterFI( 0 )->Bx()
                                 + trans( K * mMasterFI( 0 )->Bx() ) * mNormal * mMasterFI( 0 )->N()
                                 + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();



        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                       moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // compute conductivity matrix
            Matrix< DDRMat > K;
            eye( mSpaceDim, mSpaceDim, K );
            K = mMasterProp( 0 )->val()( 0 ) * K;

            // set the residual size
            this->set_residual( aResidual );

            // compute the residual r_T
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * dot( K * mMasterFI( 0 )->gradx( 1 ), mNormal )
                           + trans( K * mMasterFI( 0 )->Bx() ) * mNormal * ( mMasterFI( 0 )->val()( 0 ) - mMasterProp( 1 )->val()( 0 ) )
                           + mGamma * trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->val()( 0 ) - mMasterProp( 1 )->val()( 0 ) );

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian j_T_T
            aJacobians( 0 )( 0 ) = - trans( mMasterFI( 0 )->N() ) * trans( mNormal ) * K * mMasterFI( 0 )->Bx()
                                 + trans( K * mMasterFI( 0 )->Bx() ) * mNormal * mMasterFI( 0 )->N()
                                 + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
