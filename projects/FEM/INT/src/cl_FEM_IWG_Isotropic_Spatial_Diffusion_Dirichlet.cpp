
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check field interpolators, properties, constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // set residual size
            this->set_residual( aResidual );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->val();

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * mMasterCM( 0 )->traction( mNormal )
                             + mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->testTraction( mNormal ) * tJump
                             + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * tJump;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check field interpolators, properties, constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->val();

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->testTraction( mNormal ) * mMasterFI( 0 )->N()
                                 + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if dependency on the dof type
                if ( mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += -1.0 * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->testTraction( mNormal ) * mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->dPropdDOF( tDofType )
                       - mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->dPropdDOF( tDofType );
                }

                // if dependency on the dof type
                if ( mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->dTractiondDOF( tDofType, mNormal )
                       + mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->dTestTractiondDOF( tDofType, mNormal ) * tJump( 0 );
                }

                // if dependency on the dof type
                if ( mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += trans( mMasterFI( 0 )->N() ) * tJump( 0 ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->dSPdMasterDOF( tDofType );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                       moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual - Not implemeted." );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
