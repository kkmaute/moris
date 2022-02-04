//FEM/INT/src
#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Spalart_Allmaras_Turbulence_Dirichlet::IWG_Spalart_Allmaras_Turbulence_Dirichlet( sint aBeta )
        {
            // set mBeta for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            mPropertyMap[ "Viscosity" ] = static_cast< uint >( IWG_Property_Type::VISCOSITY );
            mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );
            mPropertyMap[ "Upwind" ]    = static_cast< uint >( IWG_Property_Type::UPWIND );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "Nitsche" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE );
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual viscosity FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ));

            // get the velocity dof FI
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the imposed viscosity property
            const std::shared_ptr< Property > & tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE ) );

            // get the selection matrix property
            const std::shared_ptr< Property > & tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFIViscosity->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // compute the jump
            Matrix< DDRMat > tJump = tFIViscosity->val() - tPropDirichlet->val();

            // compute the traction
            Matrix< DDRMat > tTraction;
            this->compute_traction( tTraction );

            // compute the test traction
            Matrix< DDRMat > tTestTraction;
            this->compute_testtraction( mResidualDofType( 0 ), tTestTraction );

            // compute the residual weak form
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            - tFIViscosity->N_trans() * tM * tTraction -
                            mBeta * trans( tTestTraction ) * tM * tJump +
                            tSPNitsche->val()( 0 ) * tFIViscosity->N_trans() * tM * tJump );

            // get the upwind property
            const std::shared_ptr< Property > & tPropUpwind =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::UPWIND ) );

            // upwind term
            if ( tPropUpwind )
            {
                // compute modified velocity
                Matrix< DDRMat > tModVelocity =
                        tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

                // add upwind contribution to residual
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { 0, 0 } ) -= aWStar * (
                                tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() *
                                dot( tModVelocity, mNormal ) * tM * tJump );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite(  mSet->get_residual()( 0 ) ),
                    "IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the velocity dof FI
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the dirichlet property
            const std::shared_ptr< Property > & tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE ) );

            // get the selection matrix property
            const std::shared_ptr< Property > & tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFIViscosity->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // compute the jump
            Matrix< DDRMat > tJump = tFIViscosity->val() - tPropDirichlet->val();

            // compute the test traction
            Matrix< DDRMat > tTestTraction;
            this->compute_testtraction( mResidualDofType( 0 ), tTestTraction );

            // get number of dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if residual dof type (here viscosity)
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - mBeta * trans( tTestTraction ) * tM * tFIViscosity->N() +
                                    tSPNitsche->val()( 0 ) * tFIViscosity->N_trans() * tM * tFIViscosity->N() );
                }

                // if imposed viscosity depends on dof type
                if( tPropDirichlet->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    + mBeta * trans( tTestTraction ) * tM * tPropDirichlet->dPropdDOF( tDofType )
                                    - tSPNitsche->val()( 0 ) * tFIViscosity->N_trans() * tM * tPropDirichlet->dPropdDOF( tDofType ) );
                }

                // if Nitsche SP depends on dof type
                if( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    tFIViscosity->N_trans() * tM * tJump * tSPNitsche->dSPdMasterDOF( tDofType ) );
                }

                // compute dtractiondu
                Matrix< DDRMat > tdtractiondu;
                this->compute_dtractiondu( tDofType, tdtractiondu );

                // compute the test traction derivative
                Matrix< DDRMat > tdtesttractiondu;
                this->compute_dtesttractiondu( tDofType, mResidualDofType( 0 ), tdtesttractiondu );

                // add contribution to jacobian
                mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                - tFIViscosity->N_trans() * tM * tdtractiondu
                                - mBeta * tM( 0 ) * tJump( 0 ) * tdtesttractiondu );

                // get the upwind property
                const std::shared_ptr< Property > & tPropUpwind =
                        mMasterProp( static_cast< uint >( IWG_Property_Type::UPWIND ) );

                // upwind term
                if ( tPropUpwind )
                {
                    // compute modified velocity
                    Matrix< DDRMat > tModVelocity =
                            tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

                    // if dof type is residual dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // compute dModVelocitydModViscosity
                        Matrix< DDRMat > tModVelocityDer = - mCb2 * tFIViscosity->dnNdxn( 1 ) / mSigma;

                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() * dot( tModVelocity, mNormal ) * tM * tFIViscosity->N() +
                                        tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() * tM * tJump * trans( mNormal ) * tModVelocityDer );
                    }

                    // if dof type is residual dof type
                    if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                    {
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() *
                                        tM * tJump * trans( mNormal ) * tFIVelocity->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                    {
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() *
                                        dot( tModVelocity, mNormal ) * tM * tPropDirichlet->dPropdDOF( tDofType ) );
                    }

                    // if upwind parameter depends on the dof type
                    if ( tPropUpwind->check_dof_dependency( tDofType ) )
                    {
                        // add contribution of SP to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        tFIViscosity->N_trans() * dot( tModVelocity, mNormal ) * tM * tJump *
                                        tPropUpwind->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_jacobian - Jacobian contains NAN or INF, exiting!");
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

        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_traction(
                Matrix< DDRMat > & aTraction )
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ));

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // compute the diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropKinViscosity );

            // compute the traction
            aTraction = tDiffusionCoeff * trans ( tFIViscosity->gradx( 1 ) ) * mNormal;
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_dtractiondu(
                const Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >      & adtractiondu )
        {
            // get the der FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the derivative of the divergence of the flux
            adtractiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ));

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // compute the diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropKinViscosity );

            // if derivative dof type is residual dof type
            if( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                // add contribution to dtractiondu
                adtractiondu +=
                        tDiffusionCoeff * trans( mNormal ) * tFIViscosity->dnNdxn( 1 );
            }

            // compute ddiffusiondu
            Matrix< DDRMat > tddiffusiondu;
            compute_ddiffusiondu(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropKinViscosity,
                    aDofTypes, tddiffusiondu );

            // add contribution to dtractiondu
            adtractiondu +=
                    trans( mNormal ) * tFIViscosity->gradx( 1 ) * tddiffusiondu;
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_testtraction(
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                Matrix< DDRMat >                   & aTestTraction )
        {
            // get the der FI
            Field_Interpolator * tFITest =
                    mMasterFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // init the derivative of the divergence of the flux
            aTestTraction.set_size( 1, tFITest->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // compute the diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropKinViscosity );

            // if derivative dof type is residual dof type
            if( aTestDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                // add contribution to dtractiondu
                aTestTraction +=
                        tDiffusionCoeff * trans( mNormal ) * tFIViscosity->dnNdxn( 1 );
            }

            // compute ddiffusiondu
            Matrix< DDRMat > tddiffusiondu;
            compute_ddiffusiondu(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropKinViscosity,
                    aTestDofTypes,
                    tddiffusiondu );

            // add contribution of diffusion to dtractiondu
            aTestTraction +=
                    trans( mNormal ) * tFIViscosity->gradx( 1 ) * tddiffusiondu;
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_dtesttractiondu(
                const moris::Cell< MSI::Dof_Type> & aDofTypes,
                const moris::Cell< MSI::Dof_Type> & aTestDofTypes,
                Matrix< DDRMat >                  & adtesttractiondu )
        {
            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the test dof type FI
            Field_Interpolator * tFITest =
                    mMasterFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // init the derivative of the test traction
            adtesttractiondu.set_size(
                    tFITest->get_number_of_space_time_coefficients(),
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the modified viscosity dof FI
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // get modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            if( ( tModViscosity >= 0.0 ) && ( aTestDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) ) )
            {
                // if derivative dof type is residual dof type
                if( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute ddiffusiondutest
                    Matrix< DDRMat > tddiffusiondutest;
                    compute_ddiffusiondu(
                            mResidualDofType( 0 ),
                            mMasterFIManager,
                            tPropKinViscosity,
                            aTestDofTypes,
                            tddiffusiondutest );

                    // add contribution
                    adtesttractiondu +=
                            trans( tddiffusiondutest ) * trans( mNormal ) * tFIModViscosity->dnNdxn( 1 );
                }

                // compute ddiffusiondu
                Matrix< DDRMat > tddiffusiondu;
                compute_ddiffusiondu(
                        mResidualDofType( 0 ),
                        mMasterFIManager,
                        tPropKinViscosity,
                        aDofTypes,
                        tddiffusiondu );

                // add contribution
                adtesttractiondu +=
                        trans( tFIModViscosity->dnNdxn( 1 ) ) * mNormal * tddiffusiondu;

                // FIXME assumed that second order derivative of diffusion coeff is zero
            }
            else
            {
                // FIXME compute the test traction derivative by FD
                this->compute_dtesttractiondu_FD( aDofTypes, aTestDofTypes, adtesttractiondu );
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_dtesttractiondu_FD(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                Matrix< DDRMat >                   & adtesttractiondu_FD,
                real                                 aPerturbation,
                fem::FDScheme_Type                   aFDSchemeType )
        {
            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the field interpolator for type
            Field_Interpolator* tFIDerivative =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of coefficients, fields and bases for the considered FI
            uint tDerNumDof    = tFIDerivative->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFIDerivative->get_number_of_space_time_bases();
            uint tDerNumFields = tFIDerivative->get_number_of_fields();

            // set size for derivative
            Matrix< DDRMat > tTestTractionForSize;
            this->compute_testtraction( aTestDofTypes, tTestTractionForSize );
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
                        this->compute_testtraction( aTestDofTypes, tTestTraction );

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
    } /* namespace fem */
} /* namespace moris */
