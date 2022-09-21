/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Pressure_Interface.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Pressure_Interface.hpp"
//LINALG/src
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Incompressible_NS_Pressure_Interface::IWG_Incompressible_NS_Pressure_Interface( sint aBeta )
        {
            // set mBeta for symmetric/skew symmetric Nitsche
            mBeta = aBeta;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Interface::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the velocity dof type FI
            // FIXME protect dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );
            Field_Interpolator * tFISlave =
                    mSlaveFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get master/slave fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMMasterFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tMasterWeight = tSPNitsche->val()( 1 );
            real tSlaveWeight  = tSPNitsche->val()( 2 );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIMaster->val() - tFISlave->val();

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            mBeta * tMasterWeight * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump );

            // compute slave residual
            mSet->get_residual()( 0 )(
                    { tSlaveResStartIndex, tSlaveResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            mBeta * tSlaveWeight * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Pressure_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Interface::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the velocity dof type FI
            // FIXME protect dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );
            Field_Interpolator * tFISlave =
                    mSlaveFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get master/slave fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMMasterFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tMasterWeight = tSPNitsche->val()( 1 );
            real tSlaveWeight  = tSPNitsche->val()( 2 );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIMaster->val() - tFISlave->val();

            // get number of master dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if dof type is velocity
                // FIXME protect dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    mBeta * tMasterWeight * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFIMaster->N() );

                    // add contribution to slave jacobian
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    mBeta * tSlaveWeight * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFIMaster->N() );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMMasterFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    mBeta * tMasterWeight * tCMMasterFluid->dTestTractiondDOF(
                                            tDofType, mNormal, tVelocityJump, mResidualDofType( 0 ) ) );
                }

                // if master SP depends on dof type
                if( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // get derivatives of SP
                    Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 2 );

                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    mBeta * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump * tMasterWeightDer
                                    + mBeta * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump * tSlaveWeightDer );
                }
            }

            // get number of slave dependencies
            uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master
            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                // if dof type is velocity
                // FIXME protect dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    mBeta * tMasterWeight * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFISlave->N() );

                    // add contribution to slave jacobian
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    mBeta * tSlaveWeight * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFISlave->N() );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMSlaveFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to slave jacobian
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) -= aWStar * (
                                    mBeta * tSlaveWeight * tCMSlaveFluid->dTestTractiondDOF(
                                            tDofType, mNormal, tVelocityJump, mResidualDofType( 0 ) ) );
                }

                // if master SP depends on dof type
                if( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // get derivatives of SP
                    Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 2 );

                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) -= aWStar * (
                                    mBeta * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump * tMasterWeightDer
                                    + mBeta * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump * tSlaveWeightDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Incompressible_NS_Pressure_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Interface::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Interface::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Interface::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Interface::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

