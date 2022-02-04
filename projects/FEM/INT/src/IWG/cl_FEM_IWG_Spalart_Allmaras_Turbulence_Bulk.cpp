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

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
