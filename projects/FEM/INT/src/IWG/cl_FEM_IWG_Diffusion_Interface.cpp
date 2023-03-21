/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Interface.cpp
 *
 */

#include "cl_FEM_IWG_Diffusion_Interface.hpp"
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

        IWG_Diffusion_Interface::IWG_Diffusion_Interface( sint aBeta )
        {
            // set sint for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the constitutive model pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );
            mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Interface::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // Check that master and slave CMs are different
            MORIS_ASSERT( &mMasterCM( 0 ) != &mSlaveCM( 0 ), "Master and Slave constitutive model are the same. This will cause problems " );

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the master field interpolator for the residual dof type
            Field_Interpolator* tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the slave field interpolator for the residual dof type
            Field_Interpolator* tFISlave =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMMasterDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            const std::shared_ptr< Constitutive_Model >& tCMSlaveDiffusion =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            real tNitsche      = tSPNitsche->val()( 0 );
            real tMasterWeight = tSPNitsche->val()( 1 );
            real tSlaveWeight  = tSPNitsche->val()( 2 );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelectMaster =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            const std::shared_ptr< Property >& tPropSelectSlave =
                    mSlaveProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // get selection property from master and slave
            real tMMaster = 1.0;
            if ( tPropSelectMaster != nullptr )
            {
                tMMaster = tPropSelectMaster->val()( 0 );
            }

            real tMSlave = 1.0;
            if ( tPropSelectSlave != nullptr )
            {
                tMSlave = tPropSelectSlave->val()( 0 );
            }

            // compute effective selection property
            real tM = tMMaster * tMSlave;

            // skip computing residual if selection property is zero
            if ( std::abs( tM ) < MORIS_REAL_EPS ) return;

            // evaluate average traction
            real tTraction =
                    tMasterWeight * tCMMasterDiffusion->traction( mNormal )( 0 ) +    //
                    tSlaveWeight * tCMSlaveDiffusion->traction( mNormal )( 0 );

            // evaluate temperature jump
            real tJump = tFIMaster->val()( 0 ) - tFISlave->val()( 0 );

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex } ) +=
                    aWStar * (                                                                                                           //
                            -tFIMaster->N_trans() * tM * tTraction                                                                       //
                            + mBeta * tMasterWeight * tCMMasterDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump    //
                            + tNitsche * tFIMaster->N_trans() * tM * tJump );

            // compute slave residual
            mSet->get_residual()( 0 )(
                    { tSlaveResStartIndex, tSlaveResStopIndex } ) +=
                    aWStar * (                                                                                                         //
                            +tFISlave->N_trans() * tM * tTraction                                                                      //
                            + mBeta * tSlaveWeight * tCMSlaveDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump    //
                            - tNitsche * tFISlave->N_trans() * tM * tJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Interface::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Interface::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
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
            Field_Interpolator* tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get slave field interpolator for the residual dof type
            Field_Interpolator* tFISlave =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMMasterDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model >& tCMSlaveDiffusion =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            real tNitsche      = tSPNitsche->val()( 0 );
            real tMasterWeight = tSPNitsche->val()( 1 );
            real tSlaveWeight  = tSPNitsche->val()( 2 );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            real tM = 1.0;
            if ( tPropSelect != nullptr )
            {
                tM = tPropSelect->val()( 0 );

                // skip computing residual if projection matrix is zero
                if ( std::abs( tM ) < MORIS_REAL_EPS ) return;
            }

            // evaluate temperature jump
            real tJump = tFIMaster->val()( 0 ) - tFISlave->val()( 0 );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof dependencies
            for ( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // get sub-matrices
                auto tJacMaster = mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } );

                auto tJacSlave = mSet->get_jacobian()(
                        { tSlaveResStartIndex, tSlaveResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } );

                // compute Jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJacMaster +=                                                                                                                        //
                            aWStar * (                                                                                                                   //
                                    +mBeta * tMasterWeight * tCMMasterDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tFIMaster->N()    //
                                    + tNitsche * tFIMaster->N_trans() * tM * tFIMaster->N() );

                    tJacSlave +=                                                                                                                       //
                            aWStar * (                                                                                                                 //
                                    +mBeta * tSlaveWeight * tCMSlaveDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tFIMaster->N()    //
                                    - tNitsche * tFISlave->N_trans() * tM * tFIMaster->N() );
                }

                // if dependency of constitutive models on the dof type
                if ( tCMMasterDiffusion->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    tJacMaster +=                                                                                                          //
                            aWStar * (                                                                                                     //
                                    -tFIMaster->N_trans() * tMasterWeight * tM * tCMMasterDiffusion->dTractiondDOF( tDofType, mNormal )    //
                                    + mBeta * tMasterWeight * tCMMasterDiffusion->dTestTractiondDOF( tDofType, mNormal, mResidualDofType( 0 ) ) * tM * tJump );

                    tJacSlave +=          //
                            aWStar * (    //
                                    +tFISlave->N_trans() * tMasterWeight * tM * tCMMasterDiffusion->dTractiondDOF( tDofType, mNormal ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    Matrix< DDRMat > tTractionDer =
                            tCMMasterDiffusion->traction( mNormal )( 0 ) * tMasterWeightDer    //
                            + tCMSlaveDiffusion->traction( mNormal )( 0 ) * tSlaveWeightDer;

                    // add contribution to Jacobian
                    tJacMaster +=                                                                                                                   //
                            aWStar * (                                                                                                              //
                                    -tFIMaster->N_trans() * tTractionDer                                                                            //
                                    + mBeta * tCMMasterDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump * tMasterWeightDer    //
                                    + tFIMaster->N_trans() * tM * tJump * tNitscheDer );

                    tJacSlave +=                                                                                                                  //
                            aWStar * (                                                                                                            //
                                    +tFISlave->N_trans() * tTractionDer                                                                           //
                                    + mBeta * tCMSlaveDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump * tSlaveWeightDer    //
                                    - tFISlave->N_trans() * tM * tJump * tNitscheDer );
                }
            }

            // get number of slave dof dependencies
            uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();

            // loop over slave dof dependencies
            for ( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // get dof type
                const Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                // get sub-matrices
                auto tJacMaster = mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tSlaveDepStartIndex, tSlaveDepStopIndex } );

                auto tJacSlave = mSet->get_jacobian()(
                        { tSlaveResStartIndex, tSlaveResStopIndex },
                        { tSlaveDepStartIndex, tSlaveDepStopIndex } );

                // compute Jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJacMaster +=                                                                                                                       //
                            aWStar * (                                                                                                                  //
                                    -mBeta * tMasterWeight * tCMMasterDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tFISlave->N()    //
                                    - tNitsche * tFIMaster->N_trans() * tM * tFISlave->N() );

                    tJacSlave +=                                                                                                                      //
                            aWStar * (                                                                                                                //
                                    -mBeta * tSlaveWeight * tCMSlaveDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tFISlave->N()    //
                                    + tNitsche * tFISlave->N_trans() * tM * tFISlave->N() );
                }

                // if dependency on the dof type
                if ( tCMSlaveDiffusion->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    tJacMaster +=         //
                            aWStar * (    //
                                    -tFIMaster->N_trans() * tSlaveWeight * tM * tCMSlaveDiffusion->dTractiondDOF( tDofType, mNormal ) );

                    tJacSlave +=                                                                                                        //
                            aWStar * (                                                                                                  //
                                    +tFISlave->N_trans() * tSlaveWeight * tM * tCMSlaveDiffusion->dTractiondDOF( tDofType, mNormal )    //
                                    + mBeta * tSlaveWeight * tCMSlaveDiffusion->dTestTractiondDOF( tDofType, mNormal, mResidualDofType( 0 ) ) * tM * tJump );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    Matrix< DDRMat > tTractionDer =                                            //
                            tCMMasterDiffusion->traction( mNormal )( 0 ) * tMasterWeightDer    //
                            + tCMSlaveDiffusion->traction( mNormal )( 0 ) * tSlaveWeightDer;

                    // add contribution to Jacobian
                    tJacMaster +=                                                                                                                   //
                            aWStar * (                                                                                                              //
                                    -tFIMaster->N_trans() * tTractionDer                                                                            //
                                    + mBeta * tCMMasterDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump * tMasterWeightDer    //
                                    + tFIMaster->N_trans() * tM * tJump * tNitscheDer );

                    tJacSlave +=                                                                                                                  //
                            aWStar * (                                                                                                            //
                                    +tFISlave->N_trans() * tTractionDer                                                                           //
                                    + mBeta * tCMSlaveDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump * tSlaveWeightDer    //
                                    - tFISlave->N_trans() * tM * tJump * tNitscheDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Diffusion_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Interface::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Interface::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
