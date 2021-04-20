/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.cpp
 *
 *  Created on: Feb 18, 2020
 *      Author: ritzert
 */
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Contact_Nitsche::IWG_Isotropic_Struc_Linear_Contact_Nitsche( sint aBeta )
        {
            // sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            // .resize: gives aValue:(The value to initialize the new elements with) and aCount:(new size of the Cell)
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "PenaltyContact" ]  = static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT );
            mStabilizationMap[ "MasterWeightInterface" ] = static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE );
            mStabilizationMap[ "SlaveWeightInterface" ]  = static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMMasterElasticity
            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveElasticity
            = mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPPenalty
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT ) );

            // get the master weight stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPMasterWeight
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );

            // get the master weight stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPSlaveWeight
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            //             // get master index for residual dof type
            //             uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            //
            //             // get slave index for residual dof type
            //             uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            //
            //             // get master field interpolator for the residual dof type
            //             Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            //
            //             // get slave field interpolator for the residual dof type
            //             Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            //
            //             // get indices for SP, CM and properties
            //             uint tElastLinIsoIndex  = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
            //             uint tPenIndex      = static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT);
            //             uint tMasterWeightIndex      = static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE);
            //             uint tSlaveWeightIndex      = static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE);



            //            uint tNitschetIndex = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_CONTACT );

            //             // evaluate average traction
            //             Matrix< DDRMat > tAverTraction = trans(mNormal) * (mMasterCM( tElastLinIsoIndex )->traction( mNormal ) + mSlaveCM( tElastLinIsoIndex )->traction( mNormal ));
            //
            // //moris::print( tJumpTraction, "tJumpTraction" );
            // moris::print( mMasterCM( tElastLinIsoIndex )->traction( mNormal ), "normal_Traction_masterCM" );
            // moris::print( mSlaveCM( tElastLinIsoIndex )->traction( mNormal ), "normal_Traction_slaveCM" );
            // moris::print( mMasterCM( tElastLinIsoIndex )->testTraction( mNormal ), "normal_test_Traction_masterCM" );
            // moris::print( mSlaveCM( tElastLinIsoIndex )->testTraction( mNormal ), "normal_test_Traction_slaveCM" );
            //
            //
            //
            //             // evaluate gap
            //             Matrix< DDRMat > tGap = trans(mNormal) * ( tFISlave->val() - tFIMaster->val()); // mNormal is normal on Masterside
            //
            //             // compute contact residual on slave side
            //
            //             mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) }, { 0, 0 } )
            //                 +=   + mStabilizationParam( tSlaveWeightIndex )->val()( 0 ) * (mSlaveCM( tElastLinIsoIndex )->testTraction( mNormal )) * mNormal * tGap * tWStar
            //                      + mStabilizationParam( tSlaveWeightIndex )->val()( 0 ) * trans( tFISlave -> N()) * mNormal * tAverTraction * tWStar
            //                      + mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFISlave -> N()) * mNormal * tGap * tWStar;
            //
            //// moris::print(mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) }, { 0, 0 } ),"residual_slave");
            //
            //             // compute contact residual on master side
            //             mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } )
            //                 +=   - mStabilizationParam( tMasterWeightIndex )->val()( 0 ) * (mMasterCM( tElastLinIsoIndex )->testTraction( mNormal )) * mNormal * tGap * tWStar
            //                      - mStabilizationParam( tMasterWeightIndex )->val()( 0 ) * trans( tFIMaster -> N()) * mNormal * tAverTraction * tWStar
            //                      - mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFIMaster -> N()) * mNormal * tGap * tWStar;




            // NEW NITSCHE: unsymmetric Nitsche


            // evaluate gap
            Matrix< DDRMat > tGap = trans(mNormal) * ( tFISlave->val() - tFIMaster->val()); // mNormal is normal on Masterside

            // evaluate average pressure
            Matrix< DDRMat > tAverConPressure = tSPMasterWeight->val()( 0 ) * trans(mNormal) * tCMMasterElasticity->traction( mNormal )
                                                                - tSPSlaveWeight->val()( 0 ) * trans(mNormal) * tCMSlaveElasticity->traction( mNormal );
            //moris::print(tAverConPressure,"tAverConPressure");

            // evaluate discrete test pressure on slave side
            Matrix< DDRMat > tTestPressureSlave = (-1) * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * mNormal;

            // evaluate discrete test pressure on master side
            Matrix< DDRMat > tTestPressureMaster = tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * mNormal;

            // compute contact residual on slave side

            mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } )
                                         +=   (-1) * trans( tFISlave -> N()) * mNormal * tAverConPressure * aWStar
                                         +  mBeta * tSPSlaveWeight->val()( 0 ) * tTestPressureSlave * tGap * aWStar
                                         + tSPPenalty->val()( 0 ) * trans( tFISlave -> N()) * mNormal * tGap * aWStar;

            // moris::print(mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ),"residual_slave");

            // compute contact residual on master side
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                                         +=   trans( tFIMaster -> N()) * mNormal * tAverConPressure * aWStar
                                         + mBeta * tSPMasterWeight->val()( 0 ) * tTestPressureMaster * tGap * aWStar
                                         - tSPPenalty->val()( 0 ) * trans( tFIMaster -> N()) * mNormal * tGap * aWStar;

            //moris::print(mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ),"residual_master");

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            //            // get master index for residual dof type
            //            uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            //
            //            // get slave index for residual dof type
            //            uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            //
            //            // get master field interpolator for the residual dof type
            //            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            //
            //            // get slave field interpolator for the residual dof type
            //            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            //
            //            // get indices for SP, CM and properties
            //            uint tElastLinIsoIndex  = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
            //            uint tPenIndex      = static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT );
            //            uint tMasterWeightIndex      = static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE);
            //            uint tSlaveWeightIndex      = static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE);
            //            //uint tNitscheIndex  = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_CONTACT );

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMMasterElasticity
            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveElasticity
            = mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPPenalty
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT ) );

            // get the master weight stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPMasterWeight
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );

            // get the master weight stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPSlaveWeight
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();



            //            // evaluate average traction
            //            Matrix< DDRMat > tJumpTraction = (mMasterCM( tElastLinIsoIndex )->traction( mNormal ) + mSlaveCM( tElastLinIsoIndex )->traction( mNormal ));
            //
            //
            //            // get number of master dof dependencies
            //            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            //
            //            // create quadrature of normal
            //            Matrix< DDRMat > tNormQ = mNormal * trans (mNormal);
            //
            //            // evaluate gap
            //            Matrix< DDRMat > tGap = trans(mNormal) * ( tFISlave->val() - tFIMaster->val()); // mNormal is normal on Master side
            //
            //            // compute the jacobian for direct dof dependencies
            //            if ( mResidualDofTypeRequested )
            //            {
            //            // compute jacobian fon slave side
            //            // slave slave
            //            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
            //                    { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 1 ) } )
            //                    += ( mStabilizationParam( tSlaveWeightIndex )->val()( 0 ) * (mSlaveCM( tElastLinIsoIndex )->testTraction( mNormal )) * tNormQ * tFISlave -> N()
            //                            + mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFISlave -> N()) * tNormQ * tFISlave -> N() )* tWStar;
            //
            //            // slave master
            //            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
            //                    { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 1 ) } )
            //                    += ( -mStabilizationParam( tSlaveWeightIndex )->val()( 0 ) * (mSlaveCM( tElastLinIsoIndex )->testTraction( mNormal )) * tNormQ * tFIMaster -> N()
            //                            - mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFISlave -> N()) * tNormQ * tFIMaster -> N() )* tWStar;
            //
            //            // compute the jacobian on master side
            //            //master master
            //            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
            //                    { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) } )
            //                    +=( mStabilizationParam( tMasterWeightIndex )->val()( 0 ) * (mMasterCM( tElastLinIsoIndex )->testTraction( mNormal )) * tNormQ * tFIMaster -> N()
            //                            + mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFIMaster -> N()) * tNormQ * tFIMaster -> N() )* tWStar;
            //
            //            //master slave
            //            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
            //                    { mSet->get_jac_dof_assembly_map()( tDofIndexMaster  )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 1 ) } )
            //                    += ( -mStabilizationParam( tMasterWeightIndex )->val()( 0 ) * (mMasterCM( tElastLinIsoIndex )->testTraction( mNormal )) * tNormQ * tFISlave -> N()
            //                            - mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFIMaster -> N()) * tNormQ * tFISlave -> N() )* tWStar;
            //            }
            //
            //
            //            // compute the jacobian for indirect dof dependencies through master constitutive models
            //            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            //            {
            //                // get the dof type
            //                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );
            //
            //                // get index for the dof type
            //                sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
            //
            //                // if dependency on the dof type
            //
            //                if ( mMasterCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
            //                {
            //                    // add contribution to jacobian
            //                    // master master
            //                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
            //                                          { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 ) } )
            //                    +=  -mStabilizationParam( tMasterWeightIndex )->val()( 0 ) * trans(tFIMaster -> N()) * tNormQ * mMasterCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
            //
            //                    // master slave
            //                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
            //                                          { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 ) } )
            //                    += mStabilizationParam( tMasterWeightIndex )->val()( 0 ) * trans(tFIMaster -> N()) * tNormQ * mSlaveCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
            //                }
            //            }
            //
            //            // compute the jacobian for indirect dof dependencies through slave constitutive models
            //            uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
            //            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            //            {
            //                // get dof type
            //                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );
            //
            //                // get index for dof type
            //                sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
            //
            //                // if dependency on the dof type
            //                if ( mSlaveCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
            //                {
            //                    // add contribution to jacobian
            //                    // slave master
            //                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
            //                                          { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 ) } )
            //                    += -mStabilizationParam( tSlaveWeightIndex )->val()( 0 ) * trans(tFISlave -> N()) * tNormQ * mMasterCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
            //
            //                    //slave slave
            //                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
            //                                          { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 ) } )
            //                    += mStabilizationParam( tSlaveWeightIndex )->val()( 0 ) * trans(tFISlave -> N()) * tNormQ * mSlaveCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
            //                }
            //            }

            // NEW NITSHE: unsymmetric Nitsche


            // create quadrature of normal
            Matrix< DDRMat > tNormQ = mNormal * trans (mNormal);

            // evaluate discrete test pressure on slave side
            Matrix< DDRMat > tTestPressureSlave = (-1) * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * mNormal;

            // evaluate discrete test pressure on master side
            Matrix< DDRMat > tTestPressureMaster = tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * mNormal;

            // evaluate derivative of gap with respect to slave DOFs
            Matrix< DDRMat > tdGapdDOFslave = trans(mNormal) *  tFISlave->N(); // mNormal is normal on Master side

            // evaluate derivative of gap with respect to slave DOFs
            Matrix< DDRMat > tdGapdDOFmaster =  (-1) * trans(mNormal) * tFIMaster->N(); // mNormal is normal on Master side

            // compute the jacobian for indirect dof dependencies through master constitutive models
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // evaluate derivative of test pressure with respect to slave DOFs
                Matrix< DDRMat > tdAverConPressuredDOFslave = (-1) * tSPSlaveWeight->val()( 0 ) * trans(mNormal) * tCMSlaveElasticity->dTractiondDOF( tDofType, mNormal ) ; // mNormal is normal on Master side

                // evaluate derivative of test pressure with respect to master DOFs
                Matrix< DDRMat > tdAverConPressuredDOFmaster = tSPMasterWeight->val()( 0 ) * trans(mNormal) * tCMMasterElasticity->dTractiondDOF( tDofType, mNormal ) ; // mNormal is normal on Master side

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } )
                            +=  ( mBeta * tSPMasterWeight->val()( 0 ) *  tTestPressureMaster * tdGapdDOFmaster
                                    + tSPPenalty->val()( 0 ) * trans( tFIMaster -> N()) * tNormQ * tFIMaster -> N() ) * aWStar;

                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } )
                            +=  ( mBeta * tSPMasterWeight->val()( 0 ) * tTestPressureMaster * tdGapdDOFslave
                                    - tSPPenalty->val()( 0 ) * trans( tFIMaster -> N()) * tNormQ * tFISlave -> N() ) * aWStar;
                }

                // if dependency on the dof type
                // dTractiondDOF
                if ( tCMMasterElasticity->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } )
                            += aWStar * ( trans( tFIMaster->N() ) * mNormal * tdAverConPressuredDOFmaster );
                    //+ tSPMasterWeight->val()( 0 ) * tCMMasterElasticity->dTestTractiondDOF( tDofType, mNormal ) * tJump;

                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } )
                            += aWStar * ( trans( tFIMaster->N() ) * mNormal * tdAverConPressuredDOFslave );
                }

                //                // if dependency of stabilization parameters on the dof type
                //                if ( tSPPenalty->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                //                {
                //                    // add contribution to jacobian
                //                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                //                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                //                    += aWStar * ( trans( tFIMaster->N() ) * tJump * tSPNitsche->dSPdMasterDOF( tDofType ) );
                //
                //                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                //                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                //                    -= aWStar * ( trans( tFISlave->N() ) * tJump * tSPNitsche->dSPdMasterDOF( tDofType ) );
                //                }
                //
                //                if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                //                {
                //                    // add contribution to jacobian
                //                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                //                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                //                    += aWStar * ( - trans( tFIMaster->N() ) * tCMMasterElasticity->traction( mNormal ) * tSPMasterWeight->dSPdMasterDOF( tDofType )
                //                                  + tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tJump * tSPMasterWeight->dSPdMasterDOF( tDofType ) );
                //
                //                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                //                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                //                    += aWStar * ( trans( tFISlave->N() ) * tCMMasterElasticity->traction( mNormal ) * tSPMasterWeight->dSPdMasterDOF( tDofType ) );
                //                }
                //
                //                if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                //                {
                //                    // add contribution to jacobian
                //                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                //                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                //                    -= aWStar * ( trans( tFIMaster->N() ) * tCMSlaveElasticity->traction( mNormal ) * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );
                //
                //                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                //                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                //                    += aWStar * (   trans( tFISlave->N() ) * tCMSlaveElasticity->traction( mNormal ) * tSPSlaveWeight->dSPdMasterDOF( tDofType )
                //                                  + tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tJump * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );
                //                }

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

                // evaluate derivative of test pressure with respect to slave DOFs
                Matrix< DDRMat > tdAverConPressuredDOFslave = (-1) * tSPSlaveWeight->val()( 0 ) * trans(mNormal) * tCMSlaveElasticity->dTractiondDOF( tDofType, mNormal ) ; // mNormal is normal on Master side

                // evaluate derivative of test pressure with respect to master DOFs
                Matrix< DDRMat > tdAverConPressuredDOFmaster = tSPMasterWeight->val()( 0 ) * trans(mNormal) * tCMMasterElasticity->dTractiondDOF( tDofType, mNormal ) ; // mNormal is normal on Master side


                //                // get the index for the dof type
                //                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                //                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                //                uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                            +=  ( tSPSlaveWeight->val()( 0 ) * tTestPressureSlave * tdGapdDOFmaster
                                    - tSPPenalty->val()( 0 ) * trans( tFISlave -> N()) * tNormQ * tFIMaster -> N() )* aWStar;
                    //std::cout << "hier_1" << std::endl;
                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                            +=  ( tSPSlaveWeight->val()( 0 ) * tTestPressureSlave * tdGapdDOFslave
                                    + tSPPenalty->val()( 0 ) * trans( tFISlave -> N()) * tNormQ * tFISlave -> N() ) * aWStar;
                }

                // if dependency on the dof type
                if ( tCMSlaveElasticity->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                            += (-1) * aWStar * ( trans( tFISlave->N() ) * mNormal * tdAverConPressuredDOFmaster );

                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                            += (-1) * aWStar * ( trans( tFISlave->N() ) * mNormal * tdAverConPressuredDOFslave );
                    //std::cout << "hier_2" << std::endl;

                }

                //                // if dependency of stabilization parameters on the dof type
                //                if ( tSPPenalty->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                //                {
                //                    // add contribution to jacobian
                //                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                //                                          { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                //                    += aWStar * ( trans( tFIMaster->N() ) * tJump * tSPPenalty->dSPdSlaveDOF( tDofType ) );
                //
                //                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                //                                          { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                //                    -= aWStar * ( trans( tFISlave->N() ) * tJump * tSPPenalty->dSPdSlaveDOF( tDofType ) );
                //                }
                //
                //                if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                //                {
                //                    // add contribution to jacobian
                //                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                //                                          { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                //                    += aWStar * ( - trans( tFIMaster->N() ) * tCMMasterElasticity->traction( mNormal ) * tSPMasterWeight->dSPdSlaveDOF( tDofType )
                //                                  + tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tJump * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );
                //
                //                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                //                                          { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                //                    += aWStar * ( trans( tFISlave->N() ) * tCMMasterElasticity->traction( mNormal ) * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );
                //                }
                //
                //                if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                //                {
                //                    // add contribution to jacobian
                //                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                //                                          { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                //                    -= aWStar * ( trans( tFIMaster->N() ) * tCMSlaveElasticity->traction( mNormal ) * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );
                //
                //                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                //                                          { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                //                    += aWStar * (   trans( tFISlave->N() ) * tCMSlaveElasticity->traction( mNormal ) * tSPSlaveWeight->dSPdSlaveDOF( tDofType )
                //                                  + tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tJump * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );
                //                }

            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


