
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Interface.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        IWG_Isotropic_Spatial_Diffusion_Interface::IWG_Isotropic_Spatial_Diffusion_Interface()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "DiffLinIso" ] = IWG_Constitutive_Type::DIFF_LIN_ISO;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ]      = IWG_Stabilization_Type::NITSCHE_INTERFACE;
            mStabilizationMap[ "MasterWeightInterface" ] = IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE;
            mStabilizationMap[ "SlaveWeightInterface" ]  = IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Interface::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master and slave field interpolators
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

            // set residual cell size
            this->set_residual_double( aResidual );

            // evaluate average traction
            Matrix< DDRMat > tTraction = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->val()( 0 ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->traction( mNormal )
                                       + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->val()( 0 ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

            // compute master residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * tTraction
                             + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->val()( 0 ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * tJump
                             + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * tJump;
            // compute slave residual
            aResidual( 1 ) =   trans( mSlaveFI( 0 )->N() ) * tTraction
                             + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->val()( 0 ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * tJump
                             - mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->val()( 0 ) * trans( mSlaveFI( 0 )->N() ) * tJump;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Interface::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master and slave field interpolators
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

            // set the jacobian cell size
            this->set_jacobian_double( aJacobians );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mMasterGlobalDofTypes.size();

            // evaluate temperature jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 )
            =   mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->val()( 0 ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * mMasterFI( 0 )->N()
              + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            aJacobians( 0 )( tMasterNumDofDependencies )
            = - mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->val()( 0 ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * mSlaveFI( 0 )->N()
              - mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * mSlaveFI( 0 )->N();

            aJacobians( 1 )( 0 )
            =   mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->val()( 0 ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * mMasterFI( 0 )->N()
              - mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->val()( 0 ) * trans( mSlaveFI( 0 )->N() ) * mMasterFI( 0 )->N();

            aJacobians( 1 )( tMasterNumDofDependencies )
            = - mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->val()( 0 ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * mSlaveFI( 0 )->N()
              + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->val()( 0 ) * trans( mSlaveFI( 0 )->N() ) * mSlaveFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through master constitutive models
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if dependency of constitutive models on the dof type
                if ( mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->val()( 0 ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->dTractiondDOF( tDofType, mNormal )
                       + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->val()( 0 ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->dTestTractiondDOF( tDofType, mNormal ) * tJump( 0 );

                    aJacobians( 1 )( iDOF ).matrix_data()
                    += trans( mSlaveFI( 0 )->N() ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->val()( 0 ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->dTractiondDOF( tDofType, mNormal );
                }

                // if dependency of stabilization parameters on the dof type
                if ( mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += trans( mMasterFI( 0 )->N() ) * tJump * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->dSPdMasterDOF( tDofType );

                    aJacobians( 1 )( iDOF ).matrix_data()
                    += - trans( mSlaveFI( 0 )->N() ) * tJump * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->dSPdMasterDOF( tDofType );
                }

                if ( mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->traction( mNormal ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->dSPdMasterDOF( tDofType )
                       + mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * tJump * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->dSPdMasterDOF( tDofType );

                    aJacobians( 1 )( iDOF ).matrix_data()
                    += trans( mSlaveFI( 0 )->N() ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->traction( mNormal ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->dSPdMasterDOF( tDofType );
                }

                if ( mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->traction( mNormal ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->dSPdMasterDOF( tDofType );

                    aJacobians( 1 )( iDOF ).matrix_data()
                    +=   trans( mSlaveFI( 0 )->N() ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->traction( mNormal ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->dSPdMasterDOF( tDofType )
                       + mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * tJump * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->dSPdMasterDOF( tDofType );
                }
            }

            // compute the jacobian for indirect dof dependencies through slave constitutive models
            uint tSlaveNumDofDependencies = mSlaveGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mSlaveGlobalDofTypes( iDOF );

                // if dependency on the dof type
                if ( mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->val()( 0 ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->dTractiondDOF( tDofType, mNormal );

                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    +=   trans( mSlaveFI( 0 )->N() ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->val()( 0 ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->dTractiondDOF( tDofType, mNormal )
                       + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->val()( 0 ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->dTestTractiondDOF( tDofType, mNormal ) * tJump( 0 );
                }

                // if dependency of stabilization parameters on the dof type
                if ( mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += trans( mMasterFI( 0 )->N() ) * tJump * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->dSPdSlaveDOF( tDofType );

                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += - trans( mSlaveFI( 0 )->N() ) * tJump * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) )->dSPdSlaveDOF( tDofType );
                }

                if ( mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->traction( mNormal ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->dSPdSlaveDOF( tDofType )
                       + mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * tJump * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->dSPdSlaveDOF( tDofType );

                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += trans( mSlaveFI( 0 )->N() ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->traction( mNormal ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) )->dSPdSlaveDOF( tDofType );
                }

                if ( mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->traction( mNormal ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->dSPdSlaveDOF( tDofType );

                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    +=   trans( mSlaveFI( 0 )->N() ) * mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->traction( mNormal ) * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->dSPdSlaveDOF( tDofType )
                       + mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testTraction( mNormal ) * tJump * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) )->dSPdSlaveDOF( tDofType );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Interface::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                       moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
