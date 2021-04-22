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

        IWG_Spalart_Allmaras_Turbulence_Interface::IWG_Spalart_Allmaras_Turbulence_Interface( sint aBeta )
        {
            // set beta for symmetric/unsymmetric Nitsche formulation
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Viscosity" ] = static_cast< uint >( IWG_Property_Type::VISCOSITY );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
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

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tNitsche      = tSPNitsche->val()( 0 );
            real tMasterWeight = tSPNitsche->val()( 1 );
            real tSlaveWeight  = tSPNitsche->val()( 2 );

            // evaluate average traction
            Matrix< DDRMat > tMasterTraction;
            this->compute_traction( tMasterTraction, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlaveTraction;
            this->compute_traction( tSlaveTraction, mtk::Master_Slave::SLAVE );
            Matrix< DDRMat > tTraction =
                    tMasterWeight * tMasterTraction +
                    tSlaveWeight  * tSlaveTraction;

            // evaluate test traction
            Matrix< DDRMat > tMasterTestTraction;
            this->compute_testtraction( mResidualDofType( 0 ), tMasterTestTraction, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlaveTestTraction;
            this->compute_testtraction( mResidualDofType( 0 ), tSlaveTestTraction,  mtk::Master_Slave::SLAVE );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpViscosity = tFIMaster->val() - tFISlave->val();

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            - trans( tFIMaster->N() ) * tTraction
                            - mBeta * tMasterWeight * trans( tMasterTestTraction ) * tJumpViscosity
                            + tNitsche * trans( tFIMaster->N() ) * tJumpViscosity ) ;

            // compute slave residual
            mSet->get_residual()( 0 )(
                    { tSlaveResStartIndex, tSlaveResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            + trans( tFISlave->N() ) * tTraction
                            - mBeta * tSlaveWeight * trans( tSlaveTestTraction ) * tJumpViscosity
                            - tNitsche * trans( tFISlave->N() ) * tJumpViscosity );

            // check for nan, infinity
            MORIS_ASSERT( isfinite(  mSet->get_residual()( 0 ) ),
                    "IWG_Spalart_Allmaras_Turbulence_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
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

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tNitsche      = tSPNitsche->val()( 0 );
            real tMasterWeight = tSPNitsche->val()( 1 );
            real tSlaveWeight  = tSPNitsche->val()( 2 );

            // evaluate average traction
            Matrix< DDRMat > tMasterTraction;
            this->compute_traction( tMasterTraction, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlaveTraction;
            this->compute_traction( tSlaveTraction, mtk::Master_Slave::SLAVE );
            Matrix< DDRMat > tTraction =
                    tMasterWeight * tMasterTraction +
                    tSlaveWeight  * tSlaveTraction;

            // evaluate test traction
            Matrix< DDRMat > tMasterTestTraction;
            this->compute_testtraction( mResidualDofType( 0 ), tMasterTestTraction, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlaveTestTraction;
            this->compute_testtraction( mResidualDofType( 0 ), tSlaveTestTraction,  mtk::Master_Slave::SLAVE );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpViscosity = tFIMaster->val() - tFISlave->val();

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master constitutive models
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
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - mBeta * tMasterWeight * trans( tMasterTestTraction ) * tFIMaster->N()
                                    + tNitsche * trans( tFIMaster->N() ) * tFIMaster->N() );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - mBeta * tSlaveWeight * trans( tSlaveTestTraction ) * tFIMaster->N()
                                    - tNitsche * trans( tFISlave->N() ) * tFIMaster->N() );
                }

                // evaluate derivative of traction
                Matrix< DDRMat > tMasterdTractiondu;
                this->compute_dtractiondu( tDofType, tMasterdTractiondu, mtk::Master_Slave::MASTER );

                // evaluate derivative of testtraction
                Matrix< DDRMat > tMasterdTestTractiondu;
                this->compute_dtesttractiondu( tDofType, mResidualDofType( 0 ), tMasterdTestTractiondu, mtk::Master_Slave::MASTER );

                // add contribution to jacobian
                mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                - trans( tFIMaster->N() ) * tMasterWeight * tMasterdTractiondu )
                                - mBeta * tMasterWeight * tMasterdTestTractiondu * tJumpViscosity( 0 );

                mSet->get_jacobian()(
                        { tSlaveResStartIndex,  tSlaveResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                + trans( tFISlave->N() ) * tMasterWeight * tMasterdTractiondu );

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 2 );

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFIMaster->N() ) * tJumpViscosity * tNitscheDer
                                    - trans( tFIMaster->N() ) * tMasterTraction * tMasterWeightDer
                                    - mBeta * trans( tMasterTestTraction ) * tJumpViscosity * tMasterWeightDer
                                    - trans( tFIMaster->N() ) * tSlaveTraction * tSlaveWeightDer );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - trans( tFISlave->N() ) * tJumpViscosity * tNitscheDer
                                    + trans( tFISlave->N() ) * tMasterTraction * tMasterWeightDer
                                    + trans( tFISlave->N() ) * tSlaveTraction * tSlaveWeightDer
                                    - mBeta * trans( tSlaveTestTraction ) * tJumpViscosity * tSlaveWeightDer );
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
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    + mBeta * tMasterWeight * trans( tMasterTestTraction ) * tFISlave->N()
                                    - tNitsche * trans( tFIMaster->N() ) * tFISlave->N() );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    + mBeta * tSlaveWeight * trans( tSlaveTestTraction ) * tFISlave->N()
                                    + tNitsche * trans( tFISlave->N() ) * tFISlave->N() );
                }

                // evaluate the derivative of the traction
                Matrix< DDRMat > tSlavedTractiondu;
                this->compute_dtractiondu( tDofType, tSlavedTractiondu, mtk::Master_Slave::SLAVE );

                // evaluate derivative of testtraction
                Matrix< DDRMat > tSlavedTestTractiondu;
                this->compute_dtesttractiondu( tDofType, mResidualDofType( 0 ), tSlavedTestTractiondu, mtk::Master_Slave::SLAVE );

                // add contribution to jacobian
                mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                - trans( tFIMaster->N() ) * tSlaveWeight * tSlavedTractiondu );

                mSet->get_jacobian()(
                        { tSlaveResStartIndex, tSlaveResStopIndex },
                        { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                + trans( tFISlave->N() ) * tSlaveWeight * tSlavedTractiondu )
                                - mBeta * tSlaveWeight * tSlavedTestTractiondu * tJumpViscosity( 0 );

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // get derivatives of SP
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 2 );

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    + trans( tFIMaster->N() ) * tJumpViscosity * tNitscheDer
                                    - trans( tFIMaster->N() ) * tMasterTraction * tMasterWeightDer
                                    - mBeta * trans( tMasterTestTraction ) * tJumpViscosity * tMasterWeightDer
                                    - trans( tFIMaster->N() ) * tSlaveTraction * tSlaveWeightDer );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    - trans( tFISlave->N() ) * tJumpViscosity * tNitscheDer
                                    + trans( tFISlave->N() ) * tMasterTraction * tMasterWeightDer
                                    + trans( tFISlave->N() ) * tSlaveTraction * tSlaveWeightDer
                                    - mBeta * trans( tSlaveTestTraction ) * tJumpViscosity * tSlaveWeightDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Spalart_Allmaras_Turbulence_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
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

        real IWG_Spalart_Allmaras_Turbulence_Interface::compute_diffusion_coefficient(
                mtk::Master_Slave aIsMaster )
        {
            // init diffusion coeff
            real tDiffusionTerm = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the wall distance property
            std::shared_ptr< Property > tPropKinViscosity =
                    this->get_properties( aIsMaster )( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the fluid kinematic viscosity value
            real tKinViscosity = tPropKinViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // compute diffusion term
                tDiffusionTerm = ( tKinViscosity + tModViscosity ) / mSigma;
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = this->compute_fn( aIsMaster );

                // compute diffusion term
                tDiffusionTerm = ( tKinViscosity + tModViscosity * tFn ) / mSigma;
            }

            return tDiffusionTerm;
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_ddiffusiondu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & addiffusiondu,
                mtk::Master_Slave                    aIsMaster )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( tDerDofType );

            // init ddiffusiondu
            addiffusiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the kinematic viscosity property
            std::shared_ptr< Property > tPropKinViscosity =
                    this->get_properties( aIsMaster )( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // if derivative dof type is viscosity
                if( tDerDofType == mResidualDofType( 0 )( 0 ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tFIModViscosity->N() / mSigma;
                }

                // if kinematic viscosity depends on derivative dof type
                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = this->compute_fn( aIsMaster );

                // compute dfndu
                Matrix< DDRMat > tdfndu;
                this->compute_dfndu( aDofTypes, tdfndu, aIsMaster );

                // if derivative dof type is viscosity
                if( tDerDofType == mResidualDofType( 0 )( 0 ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tFn * tFIModViscosity->N() / mSigma;
                }

                // if kinematic viscosity depends on derivative dof type
                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }

                // add contribution from fn to ddiffusiondu
                addiffusiondu += tModViscosity * tdfndu / mSigma;
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_traction(
                Matrix< DDRMat >  & aTraction,
                mtk::Master_Slave   aIsMaster )
        {
            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // compute the diffusion coefficient
            real tDiffusionCoeff = this->compute_diffusion_coefficient( aIsMaster );

            // compute the traction
            aTraction = tDiffusionCoeff * trans ( tFIModViscosity->gradx( 1 ) ) * mNormal;
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_dtractiondu(
                const Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >            & adtractiondu,
                mtk::Master_Slave             aIsMaster )
        {
            // get the der FI
            Field_Interpolator * tFIDer =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the derivative of the divergence of the flux
            adtractiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // compute the diffusion coefficient
            real tDiffusionCoeff = this->compute_diffusion_coefficient( aIsMaster );

            // if derivative dof type is residual dof type
            if( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                // add contribution to dtractiondu
                adtractiondu +=
                        tDiffusionCoeff * trans( mNormal ) * tFIModViscosity->dnNdxn( 1 );
            }

            // compute ddiffusiondu
            Matrix< DDRMat > tddiffusiondu;
            this->compute_ddiffusiondu( aDofTypes, tddiffusiondu, aIsMaster );

            // add contribution to dtractiondu
            adtractiondu +=
                    trans( mNormal ) * tFIModViscosity->gradx( 1 ) * tddiffusiondu;
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_testtraction(
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                Matrix< DDRMat >                   & aTestTraction,
                mtk::Master_Slave                    aIsMaster )
        {
            // get the der FI
            Field_Interpolator * tFITest =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // init the derivative of the divergence of the flux
            aTestTraction.set_size( 1, tFITest->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // compute the diffusion coefficient
            real tDiffusionCoeff = this->compute_diffusion_coefficient( aIsMaster );

            // if derivative dof type is residual dof type
            if( aTestDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                // add contribution to dtractiondu
                aTestTraction +=
                        tDiffusionCoeff * trans( mNormal ) * tFIModViscosity->dnNdxn( 1 );
            }

            // compute ddiffusiondu
            Matrix< DDRMat > tddiffusiondu;
            this->compute_ddiffusiondu( aTestDofTypes, tddiffusiondu, aIsMaster );

            // add contribution of diffusion to dtractiondu
            aTestTraction +=
                    trans( mNormal ) * tFIModViscosity->gradx( 1 ) * tddiffusiondu;
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_dtesttractiondu(
                const moris::Cell< MSI::Dof_Type> & aDofTypes,
                const moris::Cell< MSI::Dof_Type> & aTestDofTypes,
                Matrix< DDRMat >                  & adtesttractiondu,
                mtk::Master_Slave                   aIsMaster )
        {
            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the test dof type FI
            Field_Interpolator * tFITest =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // init the derivative of the test traction
            adtesttractiondu.set_size(
                    tFITest->get_number_of_space_time_coefficients(),
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            if( ( tModViscosity >= 0.0 ) && ( aTestDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) ) )
            {
                // if derivative dof type is residual dof type
                if( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute ddiffusiondutest
                    Matrix< DDRMat > tddiffusiondutest;
                    this->compute_ddiffusiondu( aTestDofTypes, tddiffusiondutest, aIsMaster );

                    // add contribution
                    adtesttractiondu +=
                            trans( tddiffusiondutest ) * trans( mNormal ) * tFIModViscosity->dnNdxn( 1 );
                }

                // compute ddiffusiondu
                Matrix< DDRMat > tddiffusiondu;
                this->compute_ddiffusiondu( aDofTypes, tddiffusiondu, aIsMaster );

                // add contribution
                adtesttractiondu +=
                        trans( tFIModViscosity->dnNdxn( 1 ) ) * mNormal * tddiffusiondu;

                // FIXME assumed that second order derivative of diffusion coeff is zero
            }
            else
            {
                // FIXME compute the test traction derivative by FD
                this->compute_dtesttractiondu_FD( aDofTypes, aTestDofTypes, adtesttractiondu, aIsMaster );
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_dtesttractiondu_FD(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                Matrix< DDRMat >                   & adtesttractiondu_FD,
                mtk::Master_Slave                    aIsMaster,
                real                                 aPerturbation,
                fem::FDScheme_Type                   aFDSchemeType )
        {
            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the field interpolator for type
            Field_Interpolator* tFIDerivative =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of coefficients, fields and bases for the considered FI
            uint tDerNumDof    = tFIDerivative->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFIDerivative->get_number_of_space_time_bases();
            uint tDerNumFields = tFIDerivative->get_number_of_fields();

            // set size for derivative
            Matrix< DDRMat > tTestTractionForSize;
            this->compute_testtraction( aTestDofTypes, tTestTractionForSize, aIsMaster );
            uint tNumRow = tTestTractionForSize.n_cols();
            adtesttractiondu_FD.set_size( tNumRow, tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFIDerivative->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( ( tDeltaH < 1e-12 ) && ( tDeltaH > - 1e-12 ) )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // loop over the points for FD
                    for( uint iPoint = 0; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficient
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFIDerivative->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // compute test traction
                        Matrix< DDRMat > tTestTraction;
                        this->compute_testtraction( aTestDofTypes, tTestTraction, aIsMaster );

                        // assemble the jacobian
                        adtesttractiondu_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( iPoint ) *
                                trans( tTestTraction ) /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFIDerivative->set_coeff( tCoeff );
        }

        //------------------------------------------------------------------------------

        real IWG_Spalart_Allmaras_Turbulence_Interface::compute_fn(
                mtk::Master_Slave aIsMaster )
        {
            // compute chi, chi³
            real tChi  = this->compute_chi( aIsMaster );
            real tChi3 = std::pow( tChi, 3 );

            // compute fn
            return ( mCn1 + tChi3 ) / ( mCn1 - tChi3 );
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_dfndu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfndu,
                mtk::Master_Slave                    aIsMaster )
        {
            // compute chi, chi², chi³
            real tChi = this->compute_chi( aIsMaster );
            real tChi2 = std::pow( tChi, 2 );
            real tChi3 = std::pow( tChi, 3 );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu, aIsMaster );

            // compute adfndu
            adfndu = 6.0 * mCn1 * tChi2 * tdchidu / std::pow( mCn1 - tChi3, 2 );
        }

        //------------------------------------------------------------------------------

        real IWG_Spalart_Allmaras_Turbulence_Interface::compute_chi(
                mtk::Master_Slave aIsMaster )
        {
            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the kinematic viscosity property
            std::shared_ptr< Property > tPropKinViscosity =
                    this->get_properties( aIsMaster )( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // compute chi
            return tFIModViscosity->val()( 0 ) / tPropKinViscosity->val()( 0 );
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_dchidu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adchidu,
                mtk::Master_Slave                    aIsMaster )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidu
            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the kinematic viscosity property
            std::shared_ptr< Property > tPropKinViscosity =
                    this->get_properties( aIsMaster )( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // compute chi
            real tChi = tFIModViscosity->val()( 0 ) / tPropKinViscosity->val()( 0 );

            // if residual dof type (here viscosity)
            if( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                adchidu += tDerFI->N() / tPropKinViscosity->val()( 0 );
            }

            if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
            {
                adchidu -=
                        tChi * tPropKinViscosity->dPropdDOF( aDofTypes ) /
                        tPropKinViscosity->val()( 0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
