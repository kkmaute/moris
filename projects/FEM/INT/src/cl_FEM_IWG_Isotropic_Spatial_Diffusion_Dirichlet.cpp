
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
            mMasterDofTypes = { { MSI::Dof_Type::TEMP } };

            // set the active property types
            mMasterPropTypes = { fem::Property_Type::CONDUCTIVITY,
                                 fem::Property_Type::TEMP_DIRICHLET };

            // FIXME set a penalty
            mGamma = 1.0;

            //FIXME forced diffusion parameter
            //      forced dimensions for 3D
            eye( mSpaceDim, mSpaceDim, mKappa );
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // compute conductivity
            mKappa = mMasterProp( 0 )->val()( 0 ) * mKappa;

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * dot( mKappa * mMasterFI( 0 )->gradx( 1 ), mNormal )
                           + trans( mKappa * mMasterFI( 0 )->Bx() ) * mNormal * ( mMasterFI( 0 )->val()( 0 ) - mMasterProp( 1 )->val()( 0 ) )
                           + mGamma * trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->val()( 0 ) - mMasterProp( 1 )->val()( 0 ) );
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // compute conductivity
            mKappa = mMasterProp( 0 )->val()( 0 ) * mKappa;

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian j_T_T
            aJacobians( 0 )( 0 ) = - trans( mMasterFI( 0 )->N() ) * trans( mNormal ) * mKappa * mMasterFI( 0 )->Bx()
                                 + trans( mKappa * mMasterFI( 0 )->Bx() ) * mNormal * mMasterFI( 0 )->N()
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

            // compute conductivity
            mKappa = mMasterProp( 0 )->val()( 0 ) * mKappa;

            // set the residual size
            this->set_residual( aResidual );

            // compute the residual r_T
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * dot( mKappa * mMasterFI( 0 )->gradx( 1 ), mNormal )
                           + trans( mKappa * mMasterFI( 0 )->Bx() ) * mNormal * ( mMasterFI( 0 )->val()( 0 ) - mMasterProp( 1 )->val()( 0 ) )
                           + mGamma * trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->val()( 0 ) - mMasterProp( 1 )->val()( 0 ) );

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian j_T_T
            aJacobians( 0 )( 0 ) = - trans( mMasterFI( 0 )->N() ) * trans( mNormal ) * mKappa * mMasterFI( 0 )->Bx()
                                 + trans( mKappa * mMasterFI( 0 )->Bx() ) * mNormal * mMasterFI( 0 )->N()
                                 + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
