/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.hpp"
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

        IWG_Isotropic_Struc_Linear_Contact_Penalty::IWG_Isotropic_Struc_Linear_Contact_Penalty()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            // .resize: gives aValue:(The value to initialize the new elements with) and aCount:(new size of the Cell)
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] =  static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "PenaltyContact" ]      =  static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT );
            mStabilizationMap[ "StabPenaltyContact" ]  =  static_cast< uint >( IWG_Stabilization_Type::STAB_PENALTY_CONTACT );

            //            mStabilizationMap[ "MasterWeightInterface" ] = IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE;
            //            mStabilizationMap[ "SlaveWeightInterface" ]  = IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE;
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_residual( real aWStar )
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
            std::shared_ptr< Constitutive_Model > & tCMMasterElasticity
            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            std::shared_ptr< Constitutive_Model > & tCMSlaveElasticity
            = mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPPenalty
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPStabPen
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::STAB_PENALTY_CONTACT ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // evaluate average traction
            Matrix< DDRMat > tJumpTraction = tCMMasterElasticity->traction( mNormal ) - tCMSlaveElasticity->traction( mNormal );
            Matrix< DDRMat > tJumpPressure = trans(mNormal) * tJumpTraction;

            // moris::print( tJumpTraction, "tJumpTraction" );

            // evaluate gap
            Matrix< DDRMat > tGap = trans(mNormal) * ( tFISlave->val() - tFIMaster->val()); // mNormal is normal on Masterside

            if ( tGap(0,0) > 0 )
            {
                tGap(0,0) = 0;
            }

            // compute contact residual on slave side
            mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } )
                                         +=   tSPPenalty->val()( 0 ) * trans( tFISlave -> N()) * mNormal * tGap(0,0) * aWStar;

            // compute contact residual on master side
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                                         += (-1) * tSPPenalty->val()( 0 ) * trans( tFIMaster -> N()) * mNormal * tGap(0,0) * aWStar;

            // if stabilized Penalty
            if ( tSPStabPen != nullptr )
            {
                // compute contact residual on slave side
                mSet->get_residual()( 0 )(
                        { tSlaveResStartIndex, tSlaveResStopIndex },
                        { 0, 0 } ) += aWStar * (
                                tSPStabPen->val()( 0 ) * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) )* mNormal * tJumpPressure );

                // compute contact residual on master side
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { 0, 0 } ) += aWStar * (
                                tSPStabPen->val()( 0 ) * tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * mNormal * tJumpPressure );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

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
            //             uint tPenIndex          = static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT );
            //             uint tStabPenIndex      = static_cast< uint >( IWG_Stabilization_Type::STAB_PENALTY_CONTACT );

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
            std::shared_ptr< Constitutive_Model > & tCMMasterElasticity
            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            std::shared_ptr< Constitutive_Model > & tCMSlaveElasticity
            = mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPPenalty
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPStabPen
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::STAB_PENALTY_CONTACT ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // evaluate average traction
            Matrix< DDRMat > tJumpTraction = tCMMasterElasticity->traction( mNormal ) - tCMSlaveElasticity->traction( mNormal );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // create quadrature of normal
            Matrix< DDRMat > tNormQ = mNormal * trans (mNormal);

            // evaluate gap
            Matrix< DDRMat > tGap = trans(mNormal) * ( tFISlave->val() - tFIMaster->val()); // mNormal is normal on Master side

            if ( tGap(0,0) <= 0 )
            {
                // hier
                // compute the jacobian through master constitutive models
                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                    uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                    uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                    // compute jacobian direct dependencies
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                += tSPPenalty->val()( 0 ) * trans( tFIMaster -> N()) * tNormQ * tFIMaster -> N() * aWStar;

                        mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                += (-1) * tSPPenalty->val()( 0 ) * trans( tFIMaster -> N()) * tNormQ * tFISlave -> N() * aWStar;
                    }
                }

                // compute the jacobian through slave constitutive models
                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                    uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                    uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                    // compute jacobian direct dependencies
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                                += (-1) * tSPPenalty->val()( 0 ) * trans( tFISlave -> N()) * tNormQ * tFIMaster -> N() * aWStar;

                        mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                                += tSPPenalty->val()( 0 ) * trans( tFISlave -> N()) * tNormQ * tFISlave -> N() * aWStar;
                    }
                }

                //hier

                //                 // compute the jacobian for direct dof dependencies
                //                 if ( mResidualDofTypeRequested )
                //                 {
                //                 // compute jacobian fon slave side
                //                 // slave slave
                //                 mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                //                         { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 1 ) } )
                //                         += tSPPenalty->val()( 0 ) * trans( tFISlave -> N()) * tNormQ * tFISlave -> N() * aWStar;
                //
                //                 // slave master
                //                 mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                //                         { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 1 ) } )
                //                         += (-1) * tSPPenalty->val()( 0 ) * trans( tFISlave -> N()) * tNormQ * tFIMaster -> N() * aWStar;
                //
                //
                //                 // compute the jacobian on master side
                //                 //master master
                //                 mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                //                         { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) } )
                //                         += mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFIMaster -> N()) * tNormQ * tFIMaster -> N() * tWStar;
                //
                //                 //master slave
                //                 mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                //                         { mSet->get_jac_dof_assembly_map()( tDofIndexMaster  )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 1 ) } )
                //                         += (-1) * mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFIMaster -> N()) * tNormQ * tFISlave -> N() * tWStar;
                //                 }
            }

            // if stabilized Penalty
            if ( tSPStabPen != nullptr )
            {
                //hier

                // compute the jacobian for indirect dof dependencies through master constitutive models
                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                    uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                    uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                    // if dependency on the dof type
                    if ( tCMMasterElasticity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                += tSPStabPen->val()( 0 ) * tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) )
                                * tNormQ * tCMMasterElasticity->dTractiondDOF( tDofType, mNormal ) * aWStar;

                        mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                += tSPStabPen->val()( 0 ) * tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) )
                                * tNormQ * tCMSlaveElasticity->dTractiondDOF( tDofType, mNormal ) * aWStar;
                    }
                }

                // compute the jacobian for indirect dof dependencies through slave constitutive models
                uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
                for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                    uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                    uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                    // if dependency on the dof type
                    if ( tCMSlaveElasticity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                                -= tSPStabPen->val()( 0 ) * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) )
                                * tNormQ * tCMMasterElasticity->dTractiondDOF( tDofType, mNormal ) * aWStar;

                        mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                                -= tSPStabPen->val()( 0 ) * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) )
                                * tNormQ * tCMSlaveElasticity->dTractiondDOF( tDofType, mNormal ) * aWStar;
                    }

                }
                //hier ende

                //                 // compute the jacobian for indirect dof dependencies through master constitutive models
                //                 for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                //                 {
                //                     // get the dof type
                //                     Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );
                //
                //                     // get index for the dof type
                //                     sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                //
                //                     // if dependency on the dof type
                //
                //                     if ( mMasterCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                //                     {
                //                         // add contribution to jacobian
                //
                //                         // master master
                //                         mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                //                                               { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 ) } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mMasterCM( tElastLinIsoIndex )->testTraction( mNormal )
                //                             * mMasterCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
                //
                //                         // master slave
                //                         mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                //                                               { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 ) } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mMasterCM( tElastLinIsoIndex )->testTraction( mNormal )
                //                             * mSlaveCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
                //                     }
                //                 }
                //
                //                 // compute the jacobian for indirect dof dependencies through slave constitutive models
                //                 uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
                //                 for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
                //                 {
                //                     // get dof type
                //                     Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );
                //
                //                     // get index for dof type
                //                     sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                //
                //                     // if dependency on the dof type
                //                     if ( mSlaveCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                //                     {
                //                         // add contribution to jacobian
                //
                //                         // slave master
                //                         mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                //                                               { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mSlaveCM( tElastLinIsoIndex )->testTraction( mNormal )
                //                             * mMasterCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
                //
                //                         //slave slave
                //                         mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                //                                               { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mSlaveCM( tElastLinIsoIndex )->testTraction( mNormal )
                //                             * mSlaveCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
                //                     }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

