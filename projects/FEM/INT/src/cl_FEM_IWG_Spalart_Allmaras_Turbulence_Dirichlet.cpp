//FEM/INT/src
#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        IWG_Spalart_Allmaras_Turbulence_Dirichlet::IWG_Spalart_Allmaras_Turbulence_Dirichlet()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = IWG_Property_Type::DIRICHLET;
            mPropertyMap[ "Viscosity" ] = IWG_Property_Type::VISCOSITY;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "Nitsche" ] = IWG_Stabilization_Type::NITSCHE;
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Spalart_Allmaras_Turbulence_Dirichlet::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Spalart_Allmaras_Turbulence_Dirichlet::set_property - No slave allowed." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Spalart_Allmaras_Turbulence_Dirichlet::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                    aStabilizationString;
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the imposed viscosity property
            std::shared_ptr< Property > tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE ) );

            // compute the jump
            Matrix< DDRMat > tJump = tFIViscosity->val() - tPropDirichlet->val();

            // compute the traction
            Matrix< DDRMat > tTraction;
            this->compute_traction( tTraction );

            // compute the test traction
            Matrix< DDRMat > tTestTraction;
            this->compute_testtraction( mResidualDofType, tTestTraction );

            // compute the residual weak form
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) +=
                    aWStar * (
                            - trans( tFIViscosity->N() ) * tTraction / mSigma
                            + mBeta * trans( tTestTraction ) * tJump / mSigma
                            + trans( tFIViscosity->N() ) * tSPNitsche->val()( 0 ) * tJump );
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the viscosity property
            std::shared_ptr< Property > tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE ) );

            // compute the jump
            Matrix< DDRMat > tJump = tFIViscosity->val() - tPropDirichlet->val();

            // compute the test traction
            Matrix< DDRMat > tTestTraction;
            this->compute_testtraction( mResidualDofType, tTestTraction );

            // get number of dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if residual dof type (here viscosity)
                if( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    mBeta * trans( tTestTraction  ) * tFIViscosity->N() / mSigma
                                    + tSPNitsche->val()( 0 ) * trans( tFIViscosity->N() ) * tFIViscosity->N() );
                }

                // if imposed viscosity depends on dof type
                if( tPropDirichlet->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    mBeta * trans( tTestTraction ) * tPropDirichlet->dPropdDOF( tDofType ) / mSigma
                                    + tSPNitsche->val()( 0 ) * trans( tFIViscosity->N() ) * tPropDirichlet->dPropdDOF( tDofType ) );
                }

                // if Nitsche SP depends on dof type
                if( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFIViscosity->N() ) * tJump * tSPNitsche->dSPdMasterDOF( tDofType ) );
                }

                // compute dtractiondu
                Matrix< DDRMat > tdtractiondu;
                this->compute_dtractiondu( tDofType, tdtractiondu );

                // compute the test traction
                Matrix< DDRMat > tdtesttractiondu;
                this->compute_dtesttractiondu( mResidualDofType, tDofType, tdtesttractiondu );

                // add contribution to jacobian
                mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                - trans( tFIViscosity->N() ) * tdtractiondu / mSigma
                                + mBeta * tJump( 0 ) * tdtesttractiondu / mSigma );
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_traction( Matrix< DDRMat > & aTraction )
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the viscosity property
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // compute the traction
            aTraction =
                    ( tFIViscosity->val()( 0 ) + tPropViscosity->val()( 0 ) ) *
                    trans ( tFIViscosity->gradx( 1 ) ) * mNormal;
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_dtractiondu(
                Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >      & adtractiondu )
        {
            // get the der FI
            Field_Interpolator * tFIDer = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the derivative of the divergence of the flux
            adtractiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the viscosity property
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

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
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_testtraction(
                moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                Matrix< DDRMat >             & aTestTraction )
        {
            // get the der FI
            Field_Interpolator * tFIDer = mMasterFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // init the derivative of the divergence of the flux
            aTestTraction.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the viscosity property
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

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
        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_dtesttractiondu(
                moris::Cell< MSI::Dof_Type> & aTestDofTypes,
                moris::Cell< MSI::Dof_Type> & aDofTypes,
                Matrix< DDRMat >            & adtesttractiondu )
        {
            // get the derivative dof type FI
            Field_Interpolator * tFIDer = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the test dof type FI
            Field_Interpolator * tFITest = mMasterFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the viscosity property
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

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
