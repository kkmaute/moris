//FEM/INT/src
#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        IWG_Spalart_Allmaras_Turbulence_Interface::IWG_Spalart_Allmaras_Turbulence_Interface()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Viscosity" ] = IWG_Property_Type::VISCOSITY;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ]      = IWG_Stabilization_Type::NITSCHE_INTERFACE;
            mStabilizationMap[ "MasterWeightInterface" ] = IWG_Stabilization_Type::MASTER_WEIGHT_INTERFACE;
            mStabilizationMap[ "SlaveWeightInterface" ]  = IWG_Stabilization_Type::SLAVE_WEIGHT_INTERFACE;
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Spalart_Allmaras_Turbulence_Interface::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Spalart_Allmaras_Turbulence_Interface::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                    aStabilizationString;
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_residual( real aWStar )
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

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

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
            Matrix< DDRMat > tMasterTraction;
            this->compute_traction( tMasterTraction, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlaveTraction;
            this->compute_traction( tSlaveTraction, mtk::Master_Slave::SLAVE );
            Matrix< DDRMat > tTraction =
                    tSPMasterWeight->val()( 0 ) * tMasterTraction +
                    tSPSlaveWeight->val()( 0 )  * tSlaveTraction;

            // evaluate test traction
            Matrix< DDRMat > tMasterTestTraction;
            this->compute_testtraction( mResidualDofType, tMasterTestTraction, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlaveTestTraction;
            this->compute_testtraction( mResidualDofType, tSlaveTestTraction, mtk::Master_Slave::SLAVE );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpViscosity = tFIMaster->val() - tFISlave->val();

            // compute master residual
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) += aWStar * (
                    - trans( tFIMaster->N() ) * tTraction
                    + tSPMasterWeight->val()( 0 ) * trans( tMasterTestTraction ) * tJumpViscosity
                    + tSPNitsche->val()( 0 ) * trans( tFIMaster->N() ) * tJumpViscosity ) ;

            // compute slave residual
            mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ) += aWStar * (
                    trans( tFISlave->N() ) * tTraction
                    + tSPSlaveWeight->val()( 0 ) * trans( tSlaveTestTraction ) * tJumpViscosity
                    - tSPNitsche->val()( 0 ) * trans( tFISlave->N() ) * tJumpViscosity );
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_jacobian( real aWStar )
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

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

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
            Matrix< DDRMat > tMasterTraction;
            this->compute_traction( tMasterTraction, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlaveTraction;
            this->compute_traction( tSlaveTraction, mtk::Master_Slave::SLAVE );
            Matrix< DDRMat > tTraction =
                    tSPMasterWeight->val()( 0 ) * tMasterTraction +
                    tSPSlaveWeight->val()( 0 )  * tSlaveTraction;

            // evaluate test traction
            Matrix< DDRMat > tMasterTestTraction;
            this->compute_testtraction( mResidualDofType, tMasterTestTraction, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlaveTestTraction;
            this->compute_testtraction( mResidualDofType, tSlaveTestTraction, mtk::Master_Slave::SLAVE );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpViscosity = tFIMaster->val() - tFISlave->val();

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
                                    tSPMasterWeight->val()( 0 ) * trans( tMasterTestTraction ) * tFIMaster->N()
                                    + tSPNitsche->val()( 0 ) * trans( tFIMaster->N() ) * tFIMaster->N() );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    tSPSlaveWeight->val()( 0 ) * trans( tSlaveTestTraction ) * tFIMaster->N()
                                    - tSPNitsche->val()( 0 ) * trans( tFISlave->N() ) * tFIMaster->N() );
                }

                // evaluate derivative of traction
                Matrix< DDRMat > tMasterdTractiondu;
                this->compute_dtractiondu( tDofType, tMasterdTractiondu, mtk::Master_Slave::MASTER );

                // evaluate derivative of testtraction
                Matrix< DDRMat > tMasterdTestTractiondu;
                this->compute_dtesttractiondu( mResidualDofType, tDofType, tMasterdTestTractiondu, mtk::Master_Slave::MASTER );

                // add contribution to jacobian
                mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                - trans( tFIMaster->N() ) * tSPMasterWeight->val()( 0 ) * tMasterdTractiondu )
                                + tSPMasterWeight->val()( 0 ) * tMasterdTestTractiondu * tJumpViscosity( 0 );

                mSet->get_jacobian()(
                        { tSlaveResStartIndex,  tSlaveResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                trans( tFISlave->N() ) * tSPMasterWeight->val()( 0 ) * tMasterdTractiondu );

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFIMaster->N() ) * tJumpViscosity * tSPNitsche->dSPdMasterDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    trans( tFISlave->N() ) * tJumpViscosity * tSPNitsche->dSPdMasterDOF( tDofType ) );
                }

                if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - trans( tFIMaster->N() ) * tMasterTraction * tSPMasterWeight->dSPdMasterDOF( tDofType )
                                    + trans( tMasterTestTraction ) * tJumpViscosity * tSPMasterWeight->dSPdMasterDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tMasterTraction * tSPMasterWeight->dSPdMasterDOF( tDofType ) );
                }

                if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    trans( tFIMaster->N() ) * tSlaveTraction * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tSlaveTraction * tSPSlaveWeight->dSPdMasterDOF( tDofType )
                                    + trans( tSlaveTestTraction ) * tJumpViscosity * tSPSlaveWeight->dSPdMasterDOF( tDofType ) );
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
                                    - tSPMasterWeight->val()( 0 ) * trans( tMasterTestTraction ) * tFISlave->N()
                                    - tSPNitsche->val()( 0 ) * trans( tFIMaster->N() ) * tFISlave->N() );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    - tSPSlaveWeight->val()( 0 ) * trans( tSlaveTestTraction ) * tFISlave->N()
                                    + tSPNitsche->val()( 0 ) * trans( tFISlave->N() ) * tFISlave->N() );
                }

                // evaluate the derivative of the traction
                Matrix< DDRMat > tSlavedTractiondu;
                this->compute_dtractiondu( tDofType, tSlavedTractiondu, mtk::Master_Slave::SLAVE );

                // evaluate derivative of testtraction
                Matrix< DDRMat > tSlavedTestTractiondu;
                this->compute_dtesttractiondu( mResidualDofType, tDofType, tSlavedTestTractiondu, mtk::Master_Slave::SLAVE );

                // add contribution to jacobian
                mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) -= aWStar * (
                                trans( tFIMaster->N() ) * tSPSlaveWeight->val()( 0 ) * tSlavedTractiondu );

                mSet->get_jacobian()(
                        { tSlaveResStartIndex, tSlaveResStopIndex },
                        { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                trans( tFISlave->N() ) * tSPSlaveWeight->val()( 0 ) * tSlavedTractiondu )
                                + tSPSlaveWeight->val()( 0 ) * tSlavedTestTractiondu * tJumpViscosity( 0 );

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    trans( tFIMaster->N() ) * tJumpViscosity * tSPNitsche->dSPdSlaveDOF( tDofType ) );

                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) -= aWStar * (
                                    trans( tFISlave->N() ) * tJumpViscosity * tSPNitsche->dSPdSlaveDOF( tDofType ) );
                }

                if ( tSPMasterWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    - trans( tFIMaster->N() ) * tMasterTraction * tSPMasterWeight->dSPdSlaveDOF( tDofType )
                                    + trans( tMasterTestTraction ) * tJumpViscosity * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tMasterTraction * tSPMasterWeight->dSPdSlaveDOF( tDofType ) );
                }

                if ( tSPSlaveWeight->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) -= aWStar * (
                                    trans( tFIMaster->N() ) * tSlaveTraction * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    trans( tFISlave->N() ) * tSlaveTraction * tSPSlaveWeight->dSPdSlaveDOF( tDofType )
                                    + trans( tSlaveTestTraction ) * tJumpViscosity * tSPSlaveWeight->dSPdSlaveDOF( tDofType ) );
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Interface::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Interface::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_traction(
                Matrix< DDRMat >  & aTraction,
                mtk::Master_Slave   aIsMaster )
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity;

            // get the viscosity property
            std::shared_ptr< Property > tPropViscosity;

            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // get the residual dof FI (here viscosity)
                    tFIViscosity = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

                    // get the viscosity property
                    tPropViscosity = mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    // get the residual dof FI (here viscosity)
                    tFIViscosity = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

                    // get the viscosity property
                    tPropViscosity = mSlaveProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
                    break;
                }
                default:
                    MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Interface::compute_traction - can only be master or slave. ");
                    break;
            }

            // compute the traction
            aTraction =
                    ( tFIViscosity->val()( 0 ) + tPropViscosity->val()( 0 ) ) *
                    trans ( tFIViscosity->gradx( 1 ) ) * mNormal;
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_dtractiondu(
                Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >      & adtractiondu,
                mtk::Master_Slave       aIsMaster )
        {
            // get the der FI
            Field_Interpolator * tFIDer;

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity;

            // get the viscosity property
            std::shared_ptr< Property > tPropViscosity;

            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // get the der FI
                    tFIDer = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

                    // get the residual dof FI (here viscosity)
                    tFIViscosity = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

                    // get the viscosity property
                    tPropViscosity = mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    // get the der FI
                    tFIDer = mSlaveFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

                    // get the residual dof FI (here viscosity)
                    tFIViscosity = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

                    // get the viscosity property
                    tPropViscosity = mSlaveProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
                    break;
                }
                default:
                    MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Interface::compute_traction - can only be master or slave. ");
                    break;
            }

            // init the derivative of the divergence of the flux
            adtractiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // if derivative dof type is residual dof type
            if( aDofTypes( 0 ) == mResidualDofType( 0 ) )
            {
                // add contribution to dtractiondu
                adtractiondu.matrix_data() +=
                         trans( mNormal ) * tFIViscosity->gradx( 1 ) * tFIViscosity->N() +
                         ( tFIViscosity->val()( 0 ) + tPropViscosity->val()( 0 ) ) * trans( mNormal ) * tFIViscosity->dnNdxn( 1 );
            }

            // if viscosity property depends on derivative dof type
            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                // add contribution to dtractiondu
                adtractiondu.matrix_data() +=
                        trans( mNormal ) * tFIViscosity->gradx( 1 ) * tPropViscosity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_testtraction(
                moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                Matrix< DDRMat >             & aTestTraction,
                mtk::Master_Slave              aIsMaster )
        {
            // get the der FI
            Field_Interpolator * tFITest;

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity;

            // get the viscosity property
            std::shared_ptr< Property > tPropViscosity;

            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // get the der FI
                    tFITest = mMasterFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

                    // get the residual dof FI (here viscosity)
                    tFIViscosity = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

                    // get the viscosity property
                    tPropViscosity = mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    // get the der FI
                    tFITest = mSlaveFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

                    // get the residual dof FI (here viscosity)
                    tFIViscosity = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

                    // get the viscosity property
                    tPropViscosity = mSlaveProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
                    break;
                }
                default:
                    MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Interface::compute_traction - can only be master or slave. ");
                    break;
            }

            // init the derivative of the divergence of the flux
            aTestTraction.set_size( 1, tFITest->get_number_of_space_time_coefficients(), 0.0 );

            // if derivative dof type is residual dof type
            if( aTestDofTypes( 0 ) == mResidualDofType( 0 ) )
            {
                // add contribution to dtractiondu
                aTestTraction.matrix_data() +=
                        trans( mNormal ) * tFIViscosity->gradx( 1 ) * tFIViscosity->N() +
                        ( tFIViscosity->val()( 0 ) + tPropViscosity->val()( 0 ) ) * trans( mNormal ) * tFIViscosity->dnNdxn( 1 );
            }

            // if viscosity property depends on derivative dof type
            if( tPropViscosity->check_dof_dependency( aTestDofTypes ) )
            {
                // add contribution to dtractiondu
                aTestTraction.matrix_data() +=
                        trans( mNormal ) * tFIViscosity->gradx( 1 ) * tPropViscosity->dPropdDOF( aTestDofTypes );
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_dtesttractiondu(
                moris::Cell< MSI::Dof_Type> & aTestDofTypes,
                moris::Cell< MSI::Dof_Type> & aDofTypes,
                Matrix< DDRMat >            & adtesttractiondu,
                mtk::Master_Slave             aIsMaster )
        {
            // get the derivative dof type FI
            Field_Interpolator * tFIDer;

            // get the test dof type FI
            Field_Interpolator * tFITest;

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity;

            // get the viscosity property
            std::shared_ptr< Property > tPropViscosity;

            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // get the der FI
                    tFIDer = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

                    // get the test dof type FI
                    tFITest = mMasterFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

                    // get the residual dof FI (here viscosity)
                    tFIViscosity = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

                    // get the viscosity property
                    tPropViscosity = mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    // get the der FI
                    tFIDer = mSlaveFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

                    // get the test dof type FI
                    tFITest = mMasterFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

                    // get the residual dof FI (here viscosity)
                    tFIViscosity = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

                    // get the viscosity property
                    tPropViscosity = mSlaveProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
                    break;
                }
                default:
                    MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Interface::compute_traction - can only be master or slave. ");
                    break;
            }

            // init the derivative of the divergence of the flux
            adtesttractiondu.set_size(
                    tFITest->get_number_of_space_time_coefficients(),
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // if derivative dof type is residual dof type
            if( ( aTestDofTypes( 0 ) == mResidualDofType( 0 ) ) && ( aDofTypes( 0 ) == mResidualDofType( 0 ) ) )
            {
                adtesttractiondu.matrix_data() +=
                        trans( tFIViscosity->N() ) * trans( mNormal ) * tFIViscosity->dnNdxn( 1 )
                        + trans( tFIViscosity->dnNdxn( 1 ) ) * mNormal * tFIViscosity->N() ;

                if( tPropViscosity->check_dof_dependency( aTestDofTypes ) )
                {
                    adtesttractiondu.matrix_data() +=
                            trans( tPropViscosity->dPropdDOF( aTestDofTypes ) ) * trans( mNormal ) * tFIViscosity->dnNdxn( 1 );
                }

                if( tPropViscosity->check_dof_dependency( aDofTypes ) )
                {
                    adtesttractiondu.matrix_data() +=
                            trans( tFIViscosity->dnNdxn( 1 ) ) * mNormal * tPropViscosity->dPropdDOF( aDofTypes );
                }

            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
