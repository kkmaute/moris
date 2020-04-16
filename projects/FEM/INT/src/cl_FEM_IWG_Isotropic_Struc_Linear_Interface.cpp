
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        IWG_Isotropic_Struc_Linear_Interface::IWG_Isotropic_Struc_Linear_Interface()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IWG_Constitutive_Type::ELAST_LIN_ISO;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ]      = IWG_Stabilization_Type::NITSCHE_INTERFACE;
            mStabilizationMap[ "MasterWeightInterface" ] = IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE;
            mStabilizationMap[ "SlaveWeightInterface" ]  = IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the elasticity constitutive model
            std::shared_ptr< Constitutive_Model > tCMMasterElasticity
            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveElasticity
            = mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPNitsche
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // get the master weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPMasterWeight
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );

            // get the master weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPSlaveWeight
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );

            // evaluate average traction
            Matrix< DDRMat > tTraction = tSPMasterWeight->val()( 0 ) * tCMMasterElasticity->traction( mNormal )
                                       + tSPSlaveWeight->val()( 0 )  * tCMSlaveElasticity->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJump = tFIMaster->val() - tFISlave->val();

            // compute master residual
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
            += aWStar * ( - trans( tFIMaster->N() ) * tTraction
                          + tSPMasterWeight->val()( 0 ) * tCMMasterElasticity->testTraction( mNormal, mResidualDofType ) * tJump
                          + tSPNitsche->val()( 0 ) * trans( tFIMaster->N() ) * tJump ) ;

            // compute slave residual
            mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } )
            += aWStar * (   trans( tFISlave->N() ) * tTraction
                          + tSPSlaveWeight->val()( 0 ) * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType ) * tJump
                          - tSPNitsche->val()( 0 ) * trans( tFISlave->N() ) * tJump );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the elasticity constitutive model
            std::shared_ptr< Constitutive_Model > tCMMasterElasticity
            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveElasticity
            = mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPNitsche
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // get the master weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPMasterWeight
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );

            // get the master weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPSlaveWeight
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // evaluate displacement jump
            Matrix< DDRMat > tJump = tFIMaster->val() - tFISlave->val();

            // compute the jacobian for indirect dof dependencies through master constitutive models
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * (   tSPMasterWeight->val()( 0 ) * tCMMasterElasticity->testTraction( mNormal, mResidualDofType ) * tFIMaster->N()
                                  + tSPNitsche->val()( 0 ) * trans( tFIMaster->N() ) * tFIMaster->N() );

                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * (   tSPSlaveWeight->val()( 0 ) * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType ) * tFIMaster->N()
                                  - tSPNitsche->val()( 0 ) * trans( tFISlave->N() ) * tFIMaster->N() );
                }

                // if dependency on the dof type
                if ( tCMMasterElasticity->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * ( - trans( tFIMaster->N() ) * tSPMasterWeight->val()( 0 ) * tCMMasterElasticity->dTractiondDOF( tDofType, mNormal ) );
                       //+ tSPMasterWeight->val()( 0 ) * tCMMasterElasticity->dTestTractiondDOF( tDofType, mNormal ) * tJump;

                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * ( trans( tFISlave->N() ) * tSPMasterWeight->val()( 0 ) * tCMMasterElasticity->dTractiondDOF( tDofType, mNormal ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * ( trans( tFIMaster->N() ) * tJump * tSPNitsche->dSPdMasterDOF( tDofType ) );

                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    -= aWStar * ( trans( tFISlave->N() ) * tJump * tSPNitsche->dSPdMasterDOF( tDofType ) );
                }

                if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * ( - trans( tFIMaster->N() ) * tCMMasterElasticity->traction( mNormal ) * tSPMasterWeight->dSPdMasterDOF( tDofType )
                                  + tCMMasterElasticity->testTraction( mNormal, mResidualDofType ) * tJump * tSPMasterWeight->dSPdMasterDOF( tDofType ) );

                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * ( trans( tFISlave->N() ) * tCMMasterElasticity->traction( mNormal ) * tSPMasterWeight->dSPdMasterDOF( tDofType ) );
                }

                if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    -= aWStar * ( trans( tFIMaster->N() ) * tCMSlaveElasticity->traction( mNormal ) * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );

                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * (   trans( tFISlave->N() ) * tCMSlaveElasticity->traction( mNormal ) * tSPSlaveWeight->dSPdMasterDOF( tDofType )
                                  + tCMSlaveElasticity->testTraction( mNormal, mResidualDofType ) * tJump * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );
                }
            }

            // compute the jacobian for indirect dof dependencies through slave constitutive models
            uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                    += aWStar * ( - tSPMasterWeight->val()( 0 ) * tCMMasterElasticity->testTraction( mNormal, mResidualDofType ) * tFISlave->N()
                                  - tSPNitsche->val()( 0 ) * trans( tFIMaster->N() ) * tFISlave->N() );

                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                          { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                    += aWStar * ( - tSPSlaveWeight->val()( 0 ) * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType ) * tFISlave->N()
                                  + tSPNitsche->val()( 0 ) * trans( tFISlave->N() ) * tFISlave->N() );
                }

                // if dependency on the dof type
                if ( tCMSlaveElasticity->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                    -= aWStar * ( trans( tFIMaster->N() ) * tSPSlaveWeight->val()( 0 ) * tCMSlaveElasticity->dTractiondDOF( tDofType, mNormal ) );

                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                          { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                    += aWStar * ( trans( tFISlave->N() ) * tSPSlaveWeight->val()( 0 ) * tCMSlaveElasticity->dTractiondDOF( tDofType, mNormal ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                    += aWStar * ( trans( tFIMaster->N() ) * tJump * tSPNitsche->dSPdSlaveDOF( tDofType ) );

                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                          { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                    -= aWStar * ( trans( tFISlave->N() ) * tJump * tSPNitsche->dSPdSlaveDOF( tDofType ) );
                }

                if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                    += aWStar * ( - trans( tFIMaster->N() ) * tCMMasterElasticity->traction( mNormal ) * tSPMasterWeight->dSPdSlaveDOF( tDofType )
                                  + tCMMasterElasticity->testTraction( mNormal, mResidualDofType ) * tJump * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );

                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                          { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                    += aWStar * ( trans( tFISlave->N() ) * tCMMasterElasticity->traction( mNormal ) * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );
                }

                if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                    -= aWStar * ( trans( tFIMaster->N() ) * tCMSlaveElasticity->traction( mNormal ) * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );

                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                          { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                    += aWStar * (   trans( tFISlave->N() ) * tCMSlaveElasticity->traction( mNormal ) * tSPSlaveWeight->dSPdSlaveDOF( tDofType )
                                  + tCMSlaveElasticity->testTraction( mNormal, mResidualDofType ) * tJump * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Interface::compute_dRdp - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
