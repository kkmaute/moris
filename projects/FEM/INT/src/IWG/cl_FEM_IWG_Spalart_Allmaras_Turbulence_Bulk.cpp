//FEM/INT/src
#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Spalart_Allmaras_Turbulence_Bulk::IWG_Spalart_Allmaras_Turbulence_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "WallDistance" ] = static_cast< uint >( IWG_Property_Type::WALL_DISTANCE );
            mPropertyMap[ "Viscosity" ]    = static_cast< uint >( IWG_Property_Type::VISCOSITY );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "SUPG" ] = static_cast< uint >( IWG_Stabilization_Type::SUPG );
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual( real aWStar )
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

            // get the velocity FI
            // FIXME protect dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // get the wall distance property
            const std::shared_ptr< Property > & tPropWallDistance =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // compute modified velocity
            Matrix< DDRMat > tModVelocity =
                    tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

            // compute production term
            real tProductionTerm = compute_production_term(
                    mResidualDofType( 0 ),
                    { MSI::Dof_Type::VX },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance );

            // compute wall destruction term
            real tWallDestructionTerm = compute_wall_destruction_term(
                    mResidualDofType( 0 ),
                    { MSI::Dof_Type::VX },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance );

            // compute diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropKinViscosity );

            // compute residual of the strong form
            Matrix< DDRMat > tR;
            this->compute_residual_strong_form( tR );

            // get sub-matrix of residual
            auto tRes = mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex });

            // compute the residual weak form
            tRes += aWStar * (
                    tFIViscosity->N_trans() * (
                            tFIViscosity->gradt( 1 ) +
                            trans( tModVelocity ) * tFIViscosity->gradx( 1 ) -
                            tProductionTerm +
                            tWallDestructionTerm ) +
                            trans( tFIViscosity->dnNdxn( 1 ) ) * (
                                    tDiffusionCoeff * tFIViscosity->gradx( 1 ) +
                                    tModVelocity * tSPSUPG->val()( 0 ) * tR( 0 ) ) );

            // show that wall distance is negative
            if( tPropWallDistance->val()( 0 ) < 0.0 )
            {
                //MORIS_LOG( "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual - Negative or zero wall distance");
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual viscosity FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the velocity FI
            // FIXME protect dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // get the wall distance property
            const std::shared_ptr< Property > & tPropWallDistance =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // compute modified velocity
            Matrix< DDRMat > tModVelocity =
                    tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

            // compute diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropKinViscosity );

            // compute residual of the strng form
            Matrix< DDRMat > tR;
            this->compute_residual_strong_form( tR );

            // get number of dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                const sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                const uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                const uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // extract sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } );

                // if residual dof type (here viscosity)
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute dModVelocitydModViscosity
                    Matrix< DDRMat > tModVelocityDer = - mCb2 * tFIViscosity->dnNdxn( 1 ) / mSigma;

                    // compute the jacobian
                    tJac += aWStar * (
                            tFIViscosity->N_trans() * (
                                    tFIViscosity->dnNdtn( 1 ) +
                                    trans( tFIVelocity->val() - 2.0 * mCb2 * tFIViscosity->gradx( 1 ) / mSigma ) * tFIViscosity->dnNdxn( 1 ) ) +
                                    trans( tFIViscosity->dnNdxn( 1 ) ) * (
                                            tDiffusionCoeff * tFIViscosity->dnNdxn( 1 ) +
                                            tModVelocityDer * tSPSUPG->val()( 0 ) * tR( 0 ) ) );
                }
                // if velocity dof type
                // FIXME protect dof type
                else if( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // compute dModVelocitydVelocity
                    const Matrix< DDRMat > & tModVelocityDer = tFIVelocity->N();

                    // compute the jacobian
                    tJac += aWStar * (
                            tFIViscosity->N_trans() * trans( tFIViscosity->gradx( 1 ) ) * tModVelocityDer +
                            trans( tFIViscosity->dnNdxn( 1 ) ) * tModVelocityDer * tSPSUPG->val()( 0 ) * tR( 0 ) );
                }

                if( tSPSUPG->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJac += aWStar * (
                            trans( tFIViscosity->dnNdxn( 1 ) ) * tModVelocity *
                            tR( 0 ) * tSPSUPG->dSPdMasterDOF( tDofType ) );
                }

                // compute jacobian of the strong form
                Matrix< DDRMat > tJ;
                this->compute_jacobian_strong_form( tDofType, tJ );

                // compute derivative of production term
                Matrix< DDRMat > tdProductiondu;
                compute_dproductiontermdu(
                        mResidualDofType( 0 ),
                        { MSI::Dof_Type::VX },
                        mMasterFIManager,
                        tPropKinViscosity,
                        tPropWallDistance,
                        tDofType,
                        tdProductiondu );

                // compute derivative of wall destruction term
                Matrix< DDRMat > tdWallDestructiondu;
                compute_dwalldestructiontermdu(
                        mResidualDofType( 0 ),
                        { MSI::Dof_Type::VX },
                        mMasterFIManager,
                        tPropKinViscosity,
                        tPropWallDistance,
                        tDofType,
                        tdWallDestructiondu );

                // compute derivative of diffusion coefficient
                Matrix< DDRMat > tdDiffdu;
                compute_ddiffusiondu(
                        mResidualDofType( 0 ),
                        mMasterFIManager,
                        tPropKinViscosity,
                        tDofType,
                        tdDiffdu );

                // compute the jacobian
                tJac += aWStar * (
                        tFIViscosity->N_trans() *
                        ( -1.0 * tdProductiondu + tdWallDestructiondu ) +
                        trans( tFIViscosity->dnNdxn( 1 ) ) *
                        ( tFIViscosity->gradx( 1 ) * tdDiffdu + tModVelocity * tSPSUPG->val()( 0 ) * tJ ) );
            }

            // show that wall distance is negative
            if( tPropWallDistance->val()( 0 ) < 0.0 )
            {
                //MORIS_LOG( "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual - Negative or zero wall distance");
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual_strong_form(
                Matrix< DDRMat > & aR )
        {
            // get the residual viscosity FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ));

            // get the velocity FI
            // FIXME protect dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // get the wall distance property
            const std::shared_ptr< Property > & tPropWallDistance =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // compute modified velocity
            Matrix< DDRMat > tModVelocity =
                    tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

            // compute production term
            real tProductionTerm = compute_production_term(
                    mResidualDofType( 0 ),
                    { MSI::Dof_Type::VX },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance );

            // compute wall destruction term
            real tWallDestructionTerm = compute_wall_destruction_term(
                    mResidualDofType( 0 ),
                    { MSI::Dof_Type::VX },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance );

            // compute divergence of flux
            Matrix< DDRMat > tDivFlux;
            this->compute_divflux( tDivFlux );

            // compute strong form of residual
            aR = tFIViscosity->gradt( 1 ) +
                    trans( tModVelocity ) * tFIViscosity->gradx( 1 ) -
                    tProductionTerm +
                    tWallDestructionTerm -
                    tDivFlux;
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_strong_form(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & aJ )
        {
            // get the residual viscosity FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the velocity FI
            // FIXME protect dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // get the wall distance property
            const std::shared_ptr< Property > & tPropWallDistance =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // get the der FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            //init aJ
            aJ.set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // compute derivative of the divergence of flux
            Matrix< DDRMat > tddivfluxdu;
            this->compute_ddivfluxdu( aDofTypes, tddivfluxdu );

            // compute derivative of production term
            Matrix< DDRMat > tdProductiondu;
            compute_dproductiontermdu(
                    mResidualDofType( 0 ),
                    { MSI::Dof_Type::VX },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance,
                    aDofTypes,
                    tdProductiondu );

            // compute derivative of wall destruction term
            Matrix< DDRMat > tdWallDestructiondu;
            compute_dwalldestructiontermdu(
                    mResidualDofType( 0 ),
                    { MSI::Dof_Type::VX },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance,
                    aDofTypes,
                    tdWallDestructiondu );

            // add contribution to jacobian
            aJ = -1.0 * tdProductiondu + tdWallDestructiondu - tddivfluxdu;

            // if dof type is residual df type (here viscosity)
            if( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                aJ += tFIViscosity->dnNdtn( 1 ) +
                        trans( tFIVelocity->val() - 2.0 * mCb2 * tFIViscosity->gradx( 1 ) / mSigma ) * tFIViscosity->dnNdxn( 1 );// -
            }
            // if dof type is velocity dof type
            // FIXME protect dof type
            else if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
            {
                aJ += trans( tFIViscosity->gradx( 1 ) ) * tFIVelocity->N();
            }
        }

//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_production_term()
//        {
//            // init production term
//            real tProductionTerm = 0.0;
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute ft2
//                real tFt2 = this->compute_ft2();
//
//                // compute Stilde
//                real tSTilde = this->compute_stilde();
//
//                // compute production term
//                tProductionTerm = mCb1 * ( 1.0 - tFt2 ) * tSTilde * tModViscosity;
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute S
//                real tS = this->compute_s();
//
//                // compute production term
//                tProductionTerm = mCb1 * ( 1.0 - mCt3 ) * tS * tModViscosity;
//            }
//
//            return tProductionTerm;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dproductiondu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adproductiondu )
//        {
//            //
//            MSI::Dof_Type tDerDofType = aDofTypes( 0 );
//
//            // get the derivative dof FI
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( tDerDofType );
//
//            // init adproductiondu
//            adproductiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute ft2
//                real tFt2 = this->compute_ft2();
//
//                // compute dft2du
//                Matrix< DDRMat > tdft2du;
//                this->compute_dft2du( aDofTypes, tdft2du );
//
//                // compute STilde
//                real tSTilde = this->compute_stilde();
//
//                // compute dSTildedu
//                Matrix< DDRMat > tdstildedu;
//                this->compute_dstildedu( aDofTypes, tdstildedu );
//
//                // if derivative dof type is viscosity
//                if( tDerDofType == mResidualDofType( 0 ) )
//                {
//                    // add contribution to dproductiondu
//                    adproductiondu +=
//                            mCb1 * ( 1.0 - tFt2 ) * tSTilde * tFIModViscosity->N();
//                }
//
//                // add contribution to dproductiondu
//                adproductiondu +=
//                        - mCb1 * tSTilde * tModViscosity * tdft2du +
//                        mCb1 * ( 1 - tFt2 ) * tModViscosity * tdstildedu;
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute S
//                real tS = this->compute_s();
//
//                // compute dsdu
//                Matrix< DDRMat > tdsdu;
//                this->compute_dsdu( aDofTypes, tdsdu );
//
//                // if derivative dof type is viscosity
//                if( tDerDofType == mResidualDofType( 0 ) )
//                {
//                    // add contribution to dproductiondu
//                    adproductiondu +=
//                            mCb1 * ( 1.0 - mCt3 ) * tS * tFIModViscosity->N();
//                }
//
//                // add contribution to dproductiondu
//                adproductiondu +=
//                        mCb1 * ( 1.0 - mCt3 ) * tModViscosity * tdsdu;
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_wall_destruction_term()
//        {
//            // init wall destruction term
//            real tWallDestructionTerm = 0.0;
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the wall distance property
//            std::shared_ptr< Property > tPropWallDistance =
//                    mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // get the wall distance value
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // FIXME check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_wall_destruction_term - Negative or zero wall distance, exiting!");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute fw
//                real tFw = this->compute_fw();
//
//                // compute ft2
//                real tFt2 = this->compute_ft2();
//
//                // compute wall destruction term
//                tWallDestructionTerm =
//                        ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
//                        std::pow( tModViscosity / tWallDistance, 2.0 );
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute wall destruction term
//                tWallDestructionTerm = - mCw1 * std::pow( tModViscosity / tWallDistance, 2.0 );
//            }
//
//            return tWallDestructionTerm;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dwalldestructiondu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adwalldestructiondu )
//        {
//            // get derivative dof type
//            MSI::Dof_Type tDerDofType = aDofTypes( 0 );
//
//            // get the derivative dof FI
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( tDerDofType );
//
//            // init adwalldestructiondu
//            adwalldestructiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the wall distance property
//            std::shared_ptr< Property > tPropWallDistance =
//                    mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // get the wall distance value
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dwalldestructiondu - Negative or zero wall distance, exiting!");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute fw
//                real tFw = this->compute_fw();
//
//                // compute dfwfu
//                Matrix< DDRMat > tdfwdu;
//                this->compute_dfwdu( aDofTypes, tdfwdu );
//
//                // compute ft2
//                real tFt2 = this->compute_ft2();
//
//                // compute dft2fu
//                Matrix< DDRMat > tdft2du;
//                this->compute_dft2du( aDofTypes, tdft2du );
//
//                // if derivative dof type is viscosity
//                if( tDerDofType == mResidualDofType( 0 ) )
//                {
//                    // add contribution to dwalldestructiondu
//                    adwalldestructiondu +=
//                            ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
//                            2.0 * tModViscosity * tFIModViscosity->N() /
//                            std::pow( tWallDistance, 2.0 );
//                }
//
//                // add contribution from ft2 and fw to dwalldestructiondu
//                adwalldestructiondu +=
//                        ( mCw1 * tdfwdu - mCb1 * tdft2du / std::pow( mKappa, 2.0 ) ) *
//                        std::pow( tModViscosity / tWallDistance, 2.0 );
//
//                // if wall distance depends on derivative dof type
//                if( ( tPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
//                        ( tWallDistance > mEpsilon ) )
//                {
//                    //MORIS_LOG( "compute_dwalldestructiondu - Wall distance depends on dof type - not tested yet" );
//
//                    // add contribution to dwalldestructiondu
//                    adwalldestructiondu -=
//                            2.0 * ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
//                            std::pow( tModViscosity, 2.0 ) *
//                            tPropWallDistance->dPropdDOF( aDofTypes ) /
//                            std::pow( tWallDistance, 3.0 );
//                }
//            }
//            // if viscosity is negative
//            else
//            {
//                // if derivative dof type is viscosity
//                if( tDerDofType == mResidualDofType( 0 ) )
//                {
//                    // add contribution to dwalldestructiondu
//                    adwalldestructiondu -=
//                            2.0 * mCw1 * tModViscosity * tFIModViscosity->N() /
//                            std::pow( tWallDistance, 2.0 );
//                }
//
//                // if wall distance depends on derivative dof type
//                if( ( tPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
//                        ( tWallDistance > mEpsilon ) )
//                {
//                    //MORIS_LOG( "compute_dwalldestructiondu - Wall distance depends on dof type - not tested yet" );
//
//                    // add contribution to dwalldestructiondu
//                    adwalldestructiondu +=
//                            2.0 * mCw1 * std::pow( tModViscosity, 2.0 ) *
//                            tPropWallDistance->dPropdDOF( aDofTypes ) /
//                            std::pow( tWallDistance, 3.0 );
//                }
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_diffusion_coefficient()
//        {
//            // init diffusion coeff
//            real tDiffusionTerm = 0.0;
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the wall distance property
//            std::shared_ptr< Property > tPropKinViscosity =
//                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // get the fluid kinematic viscosity value
//            real tKinViscosity = tPropKinViscosity->val()( 0 );
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute diffusion term
//                tDiffusionTerm = ( tKinViscosity + tModViscosity ) / mSigma;
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute fn
//                real tFn = this->compute_fn();
//
//                // compute diffusion term
//                tDiffusionTerm = ( tKinViscosity + tModViscosity * tFn ) / mSigma;
//            }
//
//            return tDiffusionTerm;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_ddiffusiondu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & addiffusiondu )
//        {
//            // get derivative dof type
//            MSI::Dof_Type tDerDofType = aDofTypes( 0 );
//
//            // get the derivative dof FI
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( tDerDofType );
//
//            // init ddiffusiondu
//            addiffusiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the fluid  kinematic viscosity property
//            std::shared_ptr< Property > tPropKinViscosity =
//                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // if derivative dof type is viscosity
//                if( tDerDofType == mResidualDofType( 0 ) )
//                {
//                    // add contribution to ddiffusiondu
//                    addiffusiondu += tFIModViscosity->N() / mSigma;
//                }
//
//                // if kinematic viscosity depends on derivative dof type
//                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
//                {
//                    // add contribution to ddiffusiondu
//                    addiffusiondu += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
//                }
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute fn
//                real tFn = this->compute_fn();
//
//                // compute dfndu
//                Matrix< DDRMat > tdfndu;
//                this->compute_dfndu( aDofTypes, tdfndu );
//
//                // if derivative dof type is viscosity
//                if( tDerDofType == mResidualDofType( 0 ) )
//                {
//                    // add contribution to ddiffusiondu
//                    addiffusiondu += tFn * tFIModViscosity->N() / mSigma;
//                }
//
//                // if kinematic viscosity depends on derivative dof type
//                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
//                {
//                    // add contribution to ddiffusiondu
//                    addiffusiondu += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
//                }
//
//                // add contribution from fn to ddiffusiondu
//                addiffusiondu += tModViscosity * tdfndu / mSigma;
//            }
//        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_divflux(
                Matrix< DDRMat > & aDivFlux )
        {
            // get the residual viscosity FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // compute diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropKinViscosity );

            // compute diffusion coefficient space derivative
            Matrix<DDRMat > tdDiffusionCoeffdx;
            compute_ddiffusiondx(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropKinViscosity,
                    tdDiffusionCoeffdx );

            // FIXME get spatial dimension
            uint tSpaceDim = tFIViscosity->dnNdxn( 1 ).n_rows();

            // compute the divergence of the flux
            aDivFlux = trans ( tdDiffusionCoeffdx ) * tFIViscosity->gradx( 1 ) +
                    tDiffusionCoeff * sum( tFIViscosity->gradx( 2 )( { 0, tSpaceDim - 1 }, { 0, 0 } ) );
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_ddivfluxdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adDivFluxdu )
        {
            // get the residual viscosity FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the viscosity property
            const std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // get the der FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the derivative of the divergence of the flux
            adDivFluxdu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // compute diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropViscosity );

            // compute diffusion coefficient space derivative
            Matrix<DDRMat > tdDiffusionCoeffdx;
            compute_ddiffusiondx(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropViscosity,
                    tdDiffusionCoeffdx );

            // get spatial dimension
            uint tSpaceDim = tFIViscosity->get_space_dim();

            // get number of space time bases
            uint tNumBases = tFIViscosity->get_number_of_space_time_bases();

            // if derivative wrt to residual dof type (here viscosity)
            if( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                // get second order shape function derivatives
                Matrix< DDRMat > td2Ndx2( 1, tNumBases, 0.0 );
                for( uint iSpace = 0; iSpace < tSpaceDim; iSpace ++ )
                {
                    td2Ndx2 += tFIViscosity->dnNdxn( 2 ).get_row( iSpace );
                }

                // add contribution to derivative
                adDivFluxdu +=
                        trans( tdDiffusionCoeffdx ) * tFIViscosity->dnNdxn( 1 ) +
                        tDiffusionCoeff * td2Ndx2;
            }

            // compute derivative of diffusion coefficient
            Matrix< DDRMat > tdDiffdu;
            compute_ddiffusiondu(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropViscosity,
                    aDofTypes,
                    tdDiffdu );

            // compute diffusion coefficient space derivative
            Matrix<DDRMat > tdDiffusionCoeffdxdu;
            compute_ddiffusiondxdu(
                    mResidualDofType( 0 ),
                    mMasterFIManager,
                    tPropViscosity,
                    aDofTypes,
                    tdDiffusionCoeffdxdu );

            // add contribution from diffusion coefficient to derivative
            adDivFluxdu +=
                    trans( tFIViscosity->gradx( 1 ) ) * tdDiffusionCoeffdxdu  +
                    sum( tFIViscosity->gradx( 2 )( { 0, tSpaceDim-1 }, { 0, 0 } ) ) * tdDiffdu ;
        }

//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_wij(
//                Matrix< DDRMat > & aWij )
//        {
//            // get the velocity FI
//            // FIXME protect dof
//            Field_Interpolator * tFIVelocity =
//                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );
//
//            // get gradient of velocity
//            Matrix< DDRMat > tGradVelocity = tFIVelocity->gradx( 1 );
//
//            // switch on space dim
//            switch ( tFIVelocity->get_number_of_fields() )
//            {
//                case 2 :
//                {
//                    // init aWij = [ w11 w12 w21 w22]
//                    aWij.set_size( 4, 1, 0.0 );
//
//                    // compute Wij
//                    aWij( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
//                    aWij( 2 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
//                    break;
//                }
//                case 3 :
//                {
//                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
//                    aWij.set_size( 9, 1, 0.0 );
//
//                    // compute Wij
//                    aWij( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
//                    aWij( 2 ) = 0.5 * ( tGradVelocity( 2, 0 ) - tGradVelocity( 0, 2 ) );
//                    aWij( 3 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
//                    aWij( 5 ) = 0.5 * ( tGradVelocity( 2, 1 ) - tGradVelocity( 1, 2 ) );
//                    aWij( 6 ) = 0.5 * ( tGradVelocity( 0, 2 ) - tGradVelocity( 2, 0 ) );
//                    aWij( 7 ) = 0.5 * ( tGradVelocity( 1, 2 ) - tGradVelocity( 2, 1 ) );
//                    break;
//                }
//                default:
//                    MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_wij - space dim can only be 2 or 3" );
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dwijdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adwijdu )
//        {
//            // get the der FI
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // get the velocity FI
//            // FIXME protect dof type
//            Field_Interpolator * tFIVelocity =
//                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );
//
//            // switch on space dim
//            switch ( tFIVelocity->get_number_of_fields() )
//            {
//                case 2 :
//                {
//                    // init aWij = [ w11 w12 w21 w22]
//                    adwijdu.set_size( 4, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//                    if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
//                    {
//                        // get gradient of velocity
//                        Matrix< DDRMat > tdNdxVelocity = tFIVelocity->dnNdxn( 1 );
//
//                        // get number of bases for displacement
//                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();
//
//                        // compute adwijdu
//                        adwijdu( { 1, 1 }, { 0, tNumBases - 1 } )             =   0.5 * tdNdxVelocity.get_row( 1 );
//                        adwijdu( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 0 );
//                        adwijdu( { 2, 2 }, { 0, tNumBases - 1 } )             = - 0.5 * tdNdxVelocity.get_row( 1 );
//                        adwijdu( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 0 );
//                    }
//                    break;
//                }
//                case 3 :
//                {
//                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
//                    adwijdu.set_size( 9, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//                    if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
//                    {
//                        // get gradient of velocity
//                        Matrix< DDRMat > tdNdxVelocity = tFIVelocity->dnNdxn( 1 );
//
//                        // get number of bases for displacement
//                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();
//
//                        // compute adwijdu
//                        adwijdu( { 1, 1 }, { 0, tNumBases - 1 } )                 =   0.5 * tdNdxVelocity.get_row( 1 );
//                        adwijdu( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = - 0.5 * tdNdxVelocity.get_row( 0 );
//
//                        adwijdu( { 2, 2 }, { 0, tNumBases - 1 } )                 =   0.5 * tdNdxVelocity.get_row( 2 );
//                        adwijdu( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 0 );
//
//                        adwijdu( { 3, 3 }, { 0, tNumBases - 1 } )                 = - 0.5 * tdNdxVelocity.get_row( 1 );
//                        adwijdu( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 0 );
//
//                        adwijdu( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 2 );
//                        adwijdu( { 5, 5 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 1 );
//
//                        adwijdu( { 6, 6 }, { 0, tNumBases - 1 } )                 = - 0.5 * tdNdxVelocity.get_row( 2 );
//                        adwijdu( { 6, 6 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 0 );
//
//                        adwijdu( { 7, 7 }, { tNumBases, 2 * tNumBases - 1 } )     = - 0.5 * tdNdxVelocity.get_row( 2 );
//                        adwijdu( { 7, 7 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 1 );
//                    }
//                    break;
//                }
//                default:
//                    MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dwijdu - space dim can only be 2 or 3" );
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_s()
//        {
//            // compute wij
//            Matrix< DDRMat > tWij;
//            this->compute_wij( tWij );
//
//            // compute WijWij
//            Matrix< DDRMat > tWijWij = trans( tWij ) * tWij;
//
//            // compute s
//            return std::sqrt( 2.0 * tWijWij( 0 ) );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dsdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >             & adsdu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tDerFI =
//                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // init dsdu
//            adsdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );
//
//            // compute sbar
//            real tS = this->compute_s();
//
//            // if s is greater than zero
//            if( tS > 0.0 )
//            {
//                // compute wij
//                Matrix< DDRMat > tWij;
//                this->compute_wij( tWij );
//
//                // compute dwijdu
//                Matrix< DDRMat > tdWijdu;
//                this->compute_dwijdu( aDofTypes, tdWijdu );
//
//                // compute dsdu
//                adsdu += 2.0 * trans( tWij ) * tdWijdu / tS;
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_chi()
//        {
//            // get the residual dof FI (here viscosity)
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the density and gravity properties
//            std::shared_ptr< Property > tPropViscosity =
//                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
//
//            // compute chi
//            return tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dchidu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >             & adchidu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tDerFI =
//                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // init adchidu
//            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the residual dof FI (here viscosity)
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the density and gravity properties
//            std::shared_ptr< Property > tPropViscosity =
//                    mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );
//
//            // compute chi
//            real tChi = tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );
//
//            // if residual dof type (here viscosity)
//            if( aDofTypes( 0 ) == mResidualDofType( 0 ) )
//            {
//                adchidu += tDerFI->N() / tPropViscosity->val()( 0 );
//            }
//
//            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
//            {
//                adchidu -= tChi * tPropViscosity->dPropdDOF( aDofTypes ) / tPropViscosity->val()( 0 );
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_fv1()
//        {
//            // compute chi, chi³
//            real tChi = this->compute_chi();
//            real tChi3 = std::pow( tChi, 3.0 );
//
//            // compute fv1
//            return tChi3 / ( tChi3 + std::pow( mCv1, 3.0 ) );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dfv1du(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >             & adfv1du )
//        {
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute dchidu
//            Matrix< DDRMat > tdchidu;
//            this->compute_dchidu( aDofTypes, tdchidu );
//
//            // compute adfv1du
//            adfv1du = 3.0 * std::pow( mCv1, 3.0 ) * std::pow( tChi, 2.0 ) * tdchidu /
//                    std::pow( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), 2.0 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_fv2()
//        {
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute fv1
//            real tFv1 = this->compute_fv1();
//
//            // compute fv2
//            return 1.0 - tChi / ( 1 + tChi * tFv1 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dfv2du(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >             & adfv2du )
//        {
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute dchidu
//            Matrix< DDRMat > tdchidu;
//            this->compute_dchidu( aDofTypes, tdchidu );
//
//            // compute fv1
//            real tFv1 = this->compute_fv1();
//
//            // compute dfv1du
//            Matrix< DDRMat > tdfv1du;
//            this->compute_dfv1du( aDofTypes, tdfv1du );
//
//            // compute adfv2du
//            adfv2du = ( std::pow( tChi, 2.0 ) * tdfv1du - tdchidu ) /
//                    ( std::pow( 1.0 + tChi * tFv1, 2.0 ) );
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_fn()
//        {
//            // compute chi, chi³
//            real tChi  = this->compute_chi();
//            real tChi3 = std::pow( tChi, 3 );
//
//            // compute fn
//            return ( mCn1 + tChi3 ) / ( mCn1 - tChi3 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dfndu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adfndu )
//        {
//            // compute chi, chi², chi³
//            real tChi = this->compute_chi();
//            real tChi2 = std::pow( tChi, 2 );
//            real tChi3 = std::pow( tChi, 3 );
//
//            // compute dchidu
//            Matrix< DDRMat > tdchidu;
//            this->compute_dchidu( aDofTypes, tdchidu );
//
//            // compute adfndu
//            adfndu = 6.0 * mCn1 * tChi2 * tdchidu / std::pow( mCn1 - tChi3, 2 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_sbar()
//        {
//            // get the viscosity FI
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the wall distance property
//            std::shared_ptr< Property > tPropWallDistance =
//                    mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );
//
//            // get wall distance
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_sbar - Negative or zero wall distance, exiting!");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // compute fv2
//            real tFv2 = this->compute_fv2();
//
//            // compute s
//            return tFv2 * tFIViscosity->val()( 0 ) / std::pow( mKappa * tWallDistance, 2.0 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dsbardu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >             & adsbardu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // init dsbardu
//            adsbardu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the viscosity FI
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the wall distance properties
//            std::shared_ptr< Property > tPropWallDistance =
//                    mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );
//
//            // get the wall distance value
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // FIXME check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dsbardu - Negative or zero wall distance, exiting!");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // compute fv2
//            real tFv2 = this->compute_fv2();
//
//            // compute dfv2du
//            Matrix< DDRMat > tdfv2du;
//            this->compute_dfv2du( aDofTypes, tdfv2du );
//
//            // compute dsbardu
//            adsbardu +=
//                    tFIViscosity->val() * tdfv2du /
//                    std::pow( mKappa * tWallDistance, 2.0 );
//
//            // if dof type is residual dof type
//            if( aDofTypes( 0 ) == mResidualDofType( 0 ) )
//            {
//                // add contribution
//                adsbardu += tFv2 * tFIViscosity->N() /
//                        std::pow( mKappa * tWallDistance, 2.0 );
//            }
//
//            // if wall distance depends on derivative dof type
//            if( ( tPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
//                    ( tWallDistance > mEpsilon ) )
//            {
//                //MORIS_LOG( "compute_dsbardu - Wall distance depends on dof type - not tested!" );
//
//                // add contribution to dsbardu
//                adsbardu -=
//                        2.0 * tFv2 * tFIViscosity->val()( 0 ) * tPropWallDistance->dPropdDOF( aDofTypes ) /
//                        ( std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ) );
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_smod()
//        {
//            // compute Sbar
//            real tS = this->compute_s();
//
//            // compute S
//            real tSBar = this->compute_sbar();
//
//            // compute s
//            real tSMod = tS * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar ) /
//                    ( ( mCv3 - 2.0 * mCv2 ) * tS - tSBar );
//
//            return tSMod;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dsmoddu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adsmoddu )
//        {
//            // compute S
//            real tS = this->compute_s();
//
//            // compute dsbardu
//            Matrix< DDRMat > tdsdu;
//            this->compute_dsdu( aDofTypes, tdsdu );
//
//            // compute SBar
//            real tSBar = this->compute_sbar();
//
//            // compute dsdu
//            Matrix< DDRMat > tdsbardu;
//            this->compute_dsbardu( aDofTypes, tdsbardu );
//
//            // compute smod num
//            real tSModNum = tS * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar );
//
//            // compute smod deno
//            real tSModDeno = ( mCv3 - 2.0 * mCv2 ) * tS - tSBar;
//
//            // compute dsmoddu
//            adsmoddu = ( ( tdsdu * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar ) +
//                    tS * ( std::pow( mCv2, 2 ) * tdsdu + mCv3 * tdsbardu ) ) * tSModDeno -
//                    tSModNum * ( ( mCv3 - 2.0 * mCv2 ) * tdsdu - tdsbardu ) ) /
//                            std::pow( tSModDeno, 2 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_stilde()
//        {
//            // compute S
//            real tS = this->compute_s();
//
//            // compute SBar
//            real tSBar = this->compute_sbar();
//
//            // init stilde
//            real tSTilde = tS;
//
//            // compute STilde
//            if( tSBar >= - mCv2 * tS )
//            {
//                tSTilde += tSBar;
//            }
//            else
//            {
//                // compute sMod
//                tSTilde += this->compute_smod();
//            }
//
//            // return STilde
//            return tSTilde;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dstildedu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adstildedu )
//        {
//            // compute SBar
//            real tS = this->compute_s();
//
//            // compute dSBardu
//            Matrix< DDRMat > tdSdu;
//            this->compute_dsdu( aDofTypes, tdSdu );
//
//            // compute SBar
//            real tSBar = this->compute_sbar();
//
//            // init dstildedu
//            adstildedu = tdSdu;
//
//            // compute dstildedu
//            if( tSBar >= - mCv2 * tS )
//            {
//                // compute dSdu
//                Matrix< DDRMat > tdSBardu;
//                this->compute_dsbardu( aDofTypes, tdSBardu );
//
//                // add dsdu
//                adstildedu += tdSBardu;
//            }
//            else
//            {
//                // compute dSModdu
//                Matrix< DDRMat > tdSModdu;
//                this->compute_dsmoddu( aDofTypes, tdSModdu );
//
//                // compute sMod
//                adstildedu += tdSModdu;
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_ft2()
//        {
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute ft2
//            real tFt2 = mCt3 * std::exp( - mCt4 * std::pow( tChi, 2.0 ) );
//
//            return tFt2;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dft2du(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adft2du )
//        {
//            // compute ft2
//            real tFt2 = this->compute_ft2();
//
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute dChidu
//            Matrix< DDRMat > tdChidu;
//            this->compute_dchidu( aDofTypes, tdChidu );
//
//            // compute dft2du
//            adft2du = - mCt4 * tFt2 * 2.0 * tChi * tdChidu;
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_r()
//        {
//            // compute stilde
//            real tSTilde = this->compute_stilde();
//
//            // get the residual dof FI (here viscosity)
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get the wall distance property
//            std::shared_ptr< Property > tPropWallDistance =
//                    mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );
//
//            // get the wall distance value
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // FIXME check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_r - Negative or zero wall distance, exiting!");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // compute viscosity / ( stilde * kappa² * d² )
//            real tR = tFIViscosity->val()( 0 ) / ( tSTilde * std::pow( mKappa * tWallDistance, 2.0 ) );
//
//            // check that r is finite and greater than zero or set it to mRLim
//            Matrix<DDRMat> tRMatrix( 1, 1, tR );
//            Matrix<DDRMat> tInvRMatrix( 1, 1, 1/tR );
//            if( !isfinite( tRMatrix ) || !isfinite(tInvRMatrix) )
//            {
//                tR = mRLim;
//            }
//
//            return std::min( tR, mRLim );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_drdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adrdu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tDerFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // init adrdu
//            adrdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );
//
//            // compute r
//            real tR = this->compute_r();
//
//            // if r < 10
//            if( tR < mRLim )
//            {
//                // get the residual dof FI (here viscosity)
//                Field_Interpolator * tFIViscosity =
//                        mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//                // get the wall distance property
//                std::shared_ptr< Property > tPropWallDistance =
//                        mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );
//
//                // get the wall distance value
//                real tWallDistance = tPropWallDistance->val()( 0 );
//
////                // FIXME check negative/zero wall distance
////                MORIS_ERROR( tWallDistance > 0.0,
////                        "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_drdu - Negative or zero wall distance, exiting!");
//
//                // threshold wall distance
//                tWallDistance = std::max(tWallDistance,mEpsilon);
//
//                // compute stilde
//                real tSTilde = this->compute_stilde();
//
//                // compute dStildedu
//                Matrix< DDRMat > tdSTildedu;
//                this->compute_dstildedu( aDofTypes, tdSTildedu );
//
//                // add contribution from dStildedu
//                adrdu -= tFIViscosity->val() * tdSTildedu / std::pow( tSTilde * mKappa * tWallDistance, 2.0 );
//
//                // if dof type is viscosity
//                if( aDofTypes( 0 ) == mResidualDofType( 0 ) )
//                {
//                    // add contribution from viscosity
//                    adrdu += tSTilde * tDerFI->N() / std::pow( tSTilde * mKappa * tWallDistance, 2.0 );
//                }
//
//                // if wall distance depends on derivative dof type
//                if( ( tPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
//                        ( tWallDistance > mEpsilon ) )
//                {
//                    //MORIS_LOG( "compute_drdu - Wall distance depends on dof type - not tested!" );
//
//                    // add contribution from wall distance
//                    adrdu -= 2.0 * tFIViscosity->val()( 0 ) * tPropWallDistance->dPropdDOF( aDofTypes ) /
//                            ( tSTilde * std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ) );
//                }
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_g()
//        {
//            // compute r
//            real tR = this->compute_r();
//
//            // compute g
//            return tR + mCw2 * ( std::pow( tR, 6.0 ) - tR );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dgdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adgdu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tDerFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // init adgdu
//            adgdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );
//
//            // compute r
//            real tR = this->compute_r();
//
//            // compute drdu
//            Matrix< DDRMat > tdrdu;
//            this->compute_drdu( aDofTypes, tdrdu );
//
//            // compute adgdu
//            adgdu += ( 1.0 + mCw2 * ( 6.0 * std::pow( tR, 5.0 ) - 1.0 ) ) * tdrdu;
//        }
//
//        //------------------------------------------------------------------------------
//
//        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_fw()
//        {
//            // compute g
//            real tG = this->compute_g();
//
//            // compute fw
//            real tFw = ( 1.0 + std::pow( mCw3, 6.0 ) ) / ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) );
//            tFw = tG * std::pow( tFw, 1.0 / 6.0 );
//
//            return tFw;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dfwdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adfwdu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tDerFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // init adfwdu
//            adfwdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );
//
//            // compute g
//            real tG = this->compute_g();
//
//            // compute dgdu
//            Matrix< DDRMat > tdgdu;
//            this->compute_dgdu( aDofTypes, tdgdu );
//
//            // compute fw
//            real tFw = this->compute_fw();
//
//            // init adfwdu
//            adfwdu += ( tFw * std::pow( mCw3, 6.0 ) * tdgdu ) /
//                    ( tG * ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) ) );
//        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
