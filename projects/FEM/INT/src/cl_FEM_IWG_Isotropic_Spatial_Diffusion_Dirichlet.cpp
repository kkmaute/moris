
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
            // check field interpolators, properties, constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set residual size
            this->set_residual( aResidual );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

//            print( mMasterFI( 0 )->N(),"mMasterCM( 0 )->traction( mNormal )");
//            print( mMasterFI( 0 )->get_space(),"mMasterCM( 0 )->traction( mNormal )");

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * mMasterCM( 0 )->traction( mNormal )
                             + mMasterCM( 0 )->testTraction( mNormal ) * tJump
                             + mGamma * trans( mMasterFI( 0 )->N() ) * tJump;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check field interpolators, properties, constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = mMasterCM( 0 )->testTraction( mNormal ) * mMasterFI( 0 )->N()
                                 + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if dependency on the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += -1.0 * mMasterCM( 0 )->testTraction( mNormal ) * mMasterProp( 0 )->dPropdDOF( tDofType )
                       - mGamma * trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( tDofType );
                }

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * trans( mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal ) )
                       + mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal, tJump ) * tJump( 0 );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                       moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // check field interpolators, properties, constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set residual size
            this->set_residual( aResidual );

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * mMasterCM( 0 )->traction( mNormal )
                             + mMasterCM( 0 )->testTraction( mNormal ) * tJump
                             + mGamma * trans( mMasterFI( 0 )->N() ) * tJump;

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = mMasterCM( 0 )->testTraction( mNormal ) * mMasterFI( 0 )->N()
                                 + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if dependency on the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += -1.0 * mMasterCM( 0 )->testTraction( mNormal ) * mMasterProp( 0 )->dPropdDOF( tDofType )
                       - mGamma * trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( tDofType );
                }

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal )
                       + mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal, tJump ) * tJump( 0 ) ;
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
