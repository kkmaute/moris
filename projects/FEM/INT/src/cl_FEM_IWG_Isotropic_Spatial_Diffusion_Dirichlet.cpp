
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
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // compute flux
            Matrix< DDRMat > tFlux;
            mMasterCM( 0 )->eval_flux( tFlux );

            // compute conductivity matrix
            Matrix< DDRMat > tK;
            mMasterCM( 0 )->eval_const( tK );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * dot( tFlux, mNormal )
                             + trans( mMasterFI( 0 )->dnNdxn( 1 ) ) * tK * mNormal * tJump
                             + mGamma * trans( mMasterFI( 0 )->N() ) * tJump;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check field interpolators, properties, constitutive models
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

//            // compute flux
//            Matrix< DDRMat > tFlux;
//            mMasterCM( 0 )->eval_flux( tFlux );

            // compute conductivity matrix
            Matrix< DDRMat > tK;
            mMasterCM( 0 )->eval_const( tK );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->dnNdxn( 1 ) ) * tK * mNormal * mMasterFI( 0 )->N()
                                 + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // if dependency on the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->dnNdxn( 1 ) ) * tK * mNormal * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( iDOF ) )
                       - mGamma * trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( iDOF ) );
                }

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // evaluate stress derivative
                    Matrix< DDRMat > tdFlux;
                    mMasterCM( 0 )->eval_dFluxdDOF( mMasterGlobalDofTypes( iDOF ), tdFlux );

                    // evaluate constitutive matrix derivative
                    Matrix< DDRMat > tdK;
                    mMasterCM( 0 )->eval_dConstdDOF( mMasterGlobalDofTypes( iDOF ), tdK );

                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * trans( mNormal ) * tdFlux
                       + trans( mMasterFI( 0 )->dnNdxn( 1 ) ) * mNormal * tdK * tJump( 0 ) ;
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                       moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // check field interpolators, properties, constitutive models
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // compute flux
            Matrix< DDRMat > tFlux;
            mMasterCM( 0 )->eval_flux( tFlux );

            // compute conductivity matrix
            Matrix< DDRMat > tK;
            mMasterCM( 0 )->eval_const( tK );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * dot( tFlux, mNormal )
                             + trans( mMasterFI( 0 )->dnNdxn( 1 ) ) * tK * mNormal * tJump
                             + mGamma * trans( mMasterFI( 0 )->N() ) * tJump;

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->dnNdxn( 1 ) ) * tK * mNormal * mMasterFI( 0 )->N()
                                 + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // if dependency on the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->dnNdxn( 1 ) ) * tK * mNormal * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( iDOF ) )
                       - mGamma * trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( iDOF ) );
                }

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // evaluate stress derivatives
                    Matrix< DDRMat > tdFlux;
                    mMasterCM( 0 )->eval_dFluxdDOF( mMasterGlobalDofTypes( iDOF ), tdFlux );

                    // evaluate constitutive matrix derivatives
                    Matrix< DDRMat > tdK;
                    mMasterCM( 0 )->eval_dConstdDOF( mMasterGlobalDofTypes( iDOF ), tdK );

                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * trans( mNormal ) * tdFlux
                       + trans( mMasterFI( 0 )->dnNdxn( 1 ) ) * mNormal * tdK * tJump( 0 ) ;
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
