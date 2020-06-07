
//FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Pressure_Interface.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        IWG_Incompressible_NS_Pressure_Interface::IWG_Incompressible_NS_Pressure_Interface()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "MasterWeightInterface" ] = IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE;
            mStabilizationMap[ "SlaveWeightInterface" ]  = IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Interface::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aConstitutiveString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Pressure_Interface::set_constitutive_model - Unknown aConstitutiveString: " ) +
                    aConstitutiveString;
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(), tErrMsg.c_str() );

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------
         void IWG_Incompressible_NS_Pressure_Interface::set_stabilization_parameter(
                 std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                 std::string                                aStabilizationString )
         {
             // check that aStabilizationString makes sense
             std::string tErrMsg =
                     std::string( "IWG_Incompressible_NS_Pressure_Interface::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                     aStabilizationString;
             MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );

             // set the stabilization parameter in the stabilization parameter cell
             this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
         }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Interface::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the velocity dof type FI
            // FIXME protect dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );
            Field_Interpolator * tFISlave =
                    mSlaveFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get master/slave fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMMasterFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the master/slave weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPMasterWeight =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
            std::shared_ptr< Stabilization_Parameter > tSPSlaveWeight =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIMaster->val() - tFISlave->val();

            // compute master residual
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) += aWStar * (
                    - tSPMasterWeight->val()( 0 ) * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tVelocityJump );

            // compute slave residual
            mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ) += aWStar * (
                    - tSPSlaveWeight->val()( 0 ) * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tVelocityJump );
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Interface::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the velocity dof type FI
            // FIXME protect dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );
            Field_Interpolator * tFISlave =
                    mSlaveFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get master/slave fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMMasterFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the master/slave weight stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPMasterWeight =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
            std::shared_ptr< Stabilization_Parameter > tSPSlaveWeight =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIMaster->val() - tFISlave->val();

            // get number of master dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

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
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - tSPMasterWeight->val()( 0 ) * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tFIMaster->N() );

                    // add contribution to slave jacobian
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - tSPSlaveWeight->val()( 0 ) * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tFIMaster->N() );
               }

                // if fluid constitutive model depends on dof type
                if ( tCMMasterFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - tSPMasterWeight->val()( 0 ) * tCMMasterFluid->dTestTractiondDOF( tDofType, mNormal, tVelocityJump, mResidualDofType ) );
                }

                // if master SP depends on dof type
                if( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tVelocityJump * tSPMasterWeight->dSPdMasterDOF( tDofType ) );
                }

                // if slave SP depends on dof type
                if( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tVelocityJump * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );
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
                                    + tSPMasterWeight->val()( 0 ) * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tFISlave->N() );

                    // add contribution to slave jacobian
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    + tSPSlaveWeight->val()( 0 ) * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tFISlave->N() );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMSlaveFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to slave jacobian
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    - tSPSlaveWeight->val()( 0 ) * tCMSlaveFluid->dTestTractiondDOF( tDofType, mNormal, tVelocityJump, mResidualDofType ) );
                }

                // if master SP depends on dof type
                if( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    - trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType ) ) * tVelocityJump * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );
                }

                // if master SP depends on dof type
                if( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to master jacobian
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    - trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType ) ) * tVelocityJump * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Interface::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Interface::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Interface::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Interface::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
