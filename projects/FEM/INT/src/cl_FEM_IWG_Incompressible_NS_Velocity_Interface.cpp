
#include "cl_FEM_IWG_Incompressible_NS_Velocity_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
        IWG_Incompressible_NS_Velocity_Interface::IWG_Incompressible_NS_Velocity_Interface()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE;
            mConstitutiveMap[ "TurbulenceFluid" ]     = IWG_Constitutive_Type::FLUID_TURBULENCE;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ]      = IWG_Stabilization_Type::NITSCHE_INTERFACE;
            mStabilizationMap[ "MasterWeightInterface" ] = IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE;
            mStabilizationMap[ "SlaveWeightInterface" ]  = IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Velocity_Interface::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aConstitutiveString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Velocity_Interface::set_constitutive_model - Unknown aConstitutiveString: " ) +
                    aConstitutiveString;
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(), tErrMsg.c_str() );

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Velocity_Interface::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Velocity_Interface::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                    aStabilizationString;
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Velocity_Interface::compute_residual( real aWStar )
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

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMMasterFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // get the master weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPMasterWeight =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );

            // get the master weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPSlaveWeight =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );

            // evaluate average traction
            Matrix< DDRMat > tTractionFluid =
                    tSPMasterWeight->val()( 0 ) * tCMMasterFluid->traction( mNormal ) +
                    tSPSlaveWeight->val()( 0 )  * tCMSlaveFluid->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpVelocity = tFIMaster->val() - tFISlave->val();

            // compute master residual
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) += aWStar * (
                    - trans( tFIMaster->N() ) * tTractionFluid
                    + tSPMasterWeight->val()( 0 ) * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity
                    + tSPNitsche->val()( 0 ) * trans( tFIMaster->N() ) * tJumpVelocity ) ;

            // compute slave residual
            mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ) += aWStar * (
                    trans( tFISlave->N() ) * tTractionFluid
                    + tSPSlaveWeight->val()( 0 ) * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity
                    - tSPNitsche->val()( 0 ) * trans( tFISlave->N() ) * tJumpVelocity );

            // get the turbulence constitutive model
            std::shared_ptr< Constitutive_Model > tCMMasterTurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_TURBULENCE ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveTurbulence =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_TURBULENCE ) );

            // if turbulence
            if( tCMMasterTurbulence != nullptr && tCMSlaveTurbulence != nullptr )
            {
                // evaluate average traction
                Matrix< DDRMat > tTractionTurbulence =
                        tSPMasterWeight->val()( 0 ) * tCMMasterTurbulence->traction( mNormal ) +
                        tSPSlaveWeight->val()( 0 )  * tCMSlaveTurbulence->traction( mNormal );

                // compute master residual
                mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) += aWStar * (
                        - trans( tFIMaster->N() ) * tTractionTurbulence
                        + tSPMasterWeight->val()( 0 ) * trans( tCMMasterTurbulence->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity ) ;

                // compute slave residual
                mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ) += aWStar * (
                        trans( tFISlave->N() ) * tTractionTurbulence
                        + tSPSlaveWeight->val()( 0 ) * trans( tCMSlaveTurbulence->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity );
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Velocity_Interface::compute_jacobian( real aWStar )
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

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMMasterFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the turbulence constitutive model
            std::shared_ptr< Constitutive_Model > tCMMasterTurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_TURBULENCE ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveTurbulence =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_TURBULENCE ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // get the master weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPMasterWeight =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );

            // get the master weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPSlaveWeight =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );

            // evaluate average traction
            Matrix< DDRMat > tTractionFluid =
                    tSPMasterWeight->val()( 0 ) * tCMMasterFluid->traction( mNormal ) +
                    tSPSlaveWeight->val()( 0 )  * tCMSlaveFluid->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpVelocity = tFIMaster->val() - tFISlave->val();

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

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
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    tSPMasterWeight->val()( 0 ) * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tFIMaster->N()
                                    + tSPNitsche->val()( 0 ) * trans( tFIMaster->N() ) * tFIMaster->N() );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    tSPSlaveWeight->val()( 0 ) * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tFIMaster->N()
                                    - tSPNitsche->val()( 0 ) * trans( tFISlave->N() ) * tFIMaster->N() );
                }

                // if dependency on the dof type
                if ( tCMMasterFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - trans( tFIMaster->N() ) * tSPMasterWeight->val()( 0 ) * tCMMasterFluid->dTractiondDOF( tDofType, mNormal )
                                    + tSPMasterWeight->val()( 0 ) * tCMMasterFluid->dTestTractiondDOF( tDofType, mNormal, tJumpVelocity, mResidualDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tSPMasterWeight->val()( 0 ) * tCMMasterFluid->dTractiondDOF( tDofType, mNormal ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFIMaster->N() ) * tJumpVelocity * tSPNitsche->dSPdMasterDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    trans( tFISlave->N() ) * tJumpVelocity * tSPNitsche->dSPdMasterDOF( tDofType ) );
                }

                if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - trans( tFIMaster->N() ) * tCMMasterFluid->traction( mNormal ) * tSPMasterWeight->dSPdMasterDOF( tDofType )
                                    + trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity * tSPMasterWeight->dSPdMasterDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tCMMasterFluid->traction( mNormal ) * tSPMasterWeight->dSPdMasterDOF( tDofType ) );
                }

                if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    trans( tFIMaster->N() ) * tCMSlaveFluid->traction( mNormal ) * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tCMSlaveFluid->traction( mNormal ) * tSPSlaveWeight->dSPdMasterDOF( tDofType )
                                    + trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );
                }

                // if turbulence
                if( tCMMasterTurbulence != nullptr && tCMSlaveTurbulence != nullptr )
                {
                    // compute jacobian direct dependencies
                    if ( tDofType( 0 ) == mResidualDofType( 0 ) )
                    {
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        tSPMasterWeight->val()( 0 ) * trans( tCMMasterTurbulence->testTraction( mNormal, mResidualDofType ) ) * tFIMaster->N() );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex,  tSlaveResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        tSPSlaveWeight->val()( 0 ) * trans( tCMSlaveTurbulence->testTraction( mNormal, mResidualDofType ) ) * tFIMaster->N() );
                    }

                    // if dependency on the dof type
                    if ( tCMMasterTurbulence->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        - trans( tFIMaster->N() ) * tSPMasterWeight->val()( 0 ) * tCMMasterTurbulence->dTractiondDOF( tDofType, mNormal )
                                        + tSPMasterWeight->val()( 0 ) * tCMMasterTurbulence->dTestTractiondDOF( tDofType, mNormal, tJumpVelocity, mResidualDofType ) );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex,  tSlaveResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        trans( tFISlave->N() ) * tSPMasterWeight->val()( 0 ) * tCMMasterTurbulence->dTractiondDOF( tDofType, mNormal ) );
                    }

                    if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        - trans( tFIMaster->N() ) * tCMMasterTurbulence->traction( mNormal ) * tSPMasterWeight->dSPdMasterDOF( tDofType )
                                        + trans( tCMMasterTurbulence->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity * tSPMasterWeight->dSPdMasterDOF( tDofType ) );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex,  tSlaveResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        trans( tFISlave->N() ) * tCMMasterTurbulence->traction( mNormal ) * tSPMasterWeight->dSPdMasterDOF( tDofType ) );
                    }

                    if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        trans( tFIMaster->N() ) * tCMSlaveTurbulence->traction( mNormal ) * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex,  tSlaveResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        trans( tFISlave->N() ) * tCMSlaveTurbulence->traction( mNormal ) * tSPSlaveWeight->dSPdMasterDOF( tDofType )
                                        + trans( tCMSlaveTurbulence->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );
                    }
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
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    - tSPMasterWeight->val()( 0 ) * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tFISlave->N()
                                    - tSPNitsche->val()( 0 ) * trans( tFIMaster->N() ) * tFISlave->N() );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    - tSPSlaveWeight->val()( 0 ) * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tFISlave->N()
                                    + tSPNitsche->val()( 0 ) * trans( tFISlave->N() ) * tFISlave->N() );
                }

                // if dependency on the dof type
                if ( tCMSlaveFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) -= aWStar * (
                                    trans( tFIMaster->N() ) * tSPSlaveWeight->val()( 0 ) * tCMSlaveFluid->dTractiondDOF( tDofType, mNormal ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tSPSlaveWeight->val()( 0 ) * tCMSlaveFluid->dTractiondDOF( tDofType, mNormal )
                                    + tSPSlaveWeight->val()( 0 ) * tCMSlaveFluid->dTestTractiondDOF( tDofType, mNormal, tJumpVelocity, mResidualDofType ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    trans( tFIMaster->N() ) * tJumpVelocity * tSPNitsche->dSPdSlaveDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) -= aWStar * (
                                    trans( tFISlave->N() ) * tJumpVelocity * tSPNitsche->dSPdSlaveDOF( tDofType ) );
                }

                if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    - trans( tFIMaster->N() ) * tCMMasterFluid->traction( mNormal ) * tSPMasterWeight->dSPdSlaveDOF( tDofType )
                                    + trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tCMMasterFluid->traction( mNormal ) * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );
                }

                if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) -= aWStar * (
                                    trans( tFIMaster->N() ) * tCMSlaveFluid->traction( mNormal ) * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tCMSlaveFluid->traction( mNormal ) * tSPSlaveWeight->dSPdSlaveDOF( tDofType )
                                    + trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );
                }

                // if turbulence
                if( tCMMasterTurbulence != nullptr && tCMSlaveTurbulence != nullptr )
                {
                    // if dof type is residual dof type
                    if( tDofType( 0 ) == mResidualDofType( 0 ) )
                    {
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                        - tSPMasterWeight->val()( 0 ) * trans( tCMMasterTurbulence->testTraction( mNormal, mResidualDofType ) ) * tFISlave->N() );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex, tSlaveResStopIndex },
                                { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                        - tSPSlaveWeight->val()( 0 ) * trans( tCMSlaveTurbulence->testTraction( mNormal, mResidualDofType ) ) * tFISlave->N() );
                    }

                    // if dependency on the dof type
                    if ( tCMSlaveTurbulence->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) -= aWStar * (
                                        trans( tFIMaster->N() ) * tSPSlaveWeight->val()( 0 ) * tCMSlaveTurbulence->dTractiondDOF( tDofType, mNormal ) );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex, tSlaveResStopIndex },
                                { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                        trans( tFISlave->N() ) * tSPSlaveWeight->val()( 0 ) * tCMSlaveTurbulence->dTractiondDOF( tDofType, mNormal )
                                        + tSPSlaveWeight->val()( 0 ) * tCMSlaveTurbulence->dTestTractiondDOF( tDofType, mNormal, tJumpVelocity, mResidualDofType ) );
                    }

                    if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                        - trans( tFIMaster->N() ) * tCMMasterTurbulence->traction( mNormal ) * tSPMasterWeight->dSPdSlaveDOF( tDofType )
                                        + trans( tCMMasterTurbulence->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex, tSlaveResStopIndex },
                                { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                        trans( tFISlave->N() ) * tCMMasterTurbulence->traction( mNormal ) * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );
                    }

                    if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) -= aWStar * (
                                        trans( tFIMaster->N() ) * tCMSlaveTurbulence->traction( mNormal ) * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex, tSlaveResStopIndex },
                                { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                        trans( tFISlave->N() ) * tCMSlaveTurbulence->traction( mNormal ) * tSPSlaveWeight->dSPdSlaveDOF( tDofType )
                                        + trans( tCMSlaveTurbulence->testTraction( mNormal, mResidualDofType ) ) * tJumpVelocity * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Velocity_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Velocity_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Interface::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
