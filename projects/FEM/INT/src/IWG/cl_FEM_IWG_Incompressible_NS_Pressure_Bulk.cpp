/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Pressure_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Incompressible_NS_Pressure_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Incompressible_NS_Pressure_Bulk::IWG_Incompressible_NS_Pressure_Bulk()
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Gravity" ]          = static_cast< uint >( IWG_Property_Type::GRAVITY );
        mPropertyMap[ "ThermalExpansion" ] = static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION );
        mPropertyMap[ "ReferenceTemp" ]    = static_cast< uint >( IWG_Property_Type::REF_TEMP );
        mPropertyMap[ "InvPermeability" ]  = static_cast< uint >( IWG_Property_Type::INV_PERMEABILITY );
        mPropertyMap[ "MassSource" ]       = static_cast< uint >( IWG_Property_Type::MASS_SOURCE );
        mPropertyMap[ "Load" ]             = static_cast< uint >( IWG_Property_Type::BODY_LOAD );
        mPropertyMap[ "PressureSpring" ]   = static_cast< uint >( IWG_Property_Type::PRESSURE_SPRING );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "IncompressibleFluid" ] = static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "IncompressibleFlow" ] = static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW );
    }

    //------------------------------------------------------------------------------

    void IWG_Incompressible_NS_Pressure_Bulk::compute_residual( real aWStar )
    {
        // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here pressure), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the pressure and velocity FIs
            Field_Interpolator * tPressureFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            Field_Interpolator * tVelocityFI = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX ); //FIXME this need to be protected

            // get the mass source property
            const std::shared_ptr< Property > & tMassSourceProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::MASS_SOURCE ) );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tIncFluidCM =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the pressure spring property
            const std::shared_ptr< Property > & tPressureSpringProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESSURE_SPRING ) );

            // get the incompressible flow stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPPSPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW ) );

            // get the density
            real tDensity = tIncFluidCM->get_property( "Density" )->val()(0);

            // get residual vector
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual weak form
            tRes += aWStar * ( tPressureFI->N_trans() * tVelocityFI->div() );

            // if source term is defined
            if ( tMassSourceProp != nullptr )
            {
                // add contribution of source term
                tRes -= aWStar * ( tPressureFI->N_trans() * tMassSourceProp->val()( 0 ) / tDensity );
            }

            // if pressure spring is defined
            if ( tPressureSpringProp != nullptr )
            {
                // add contribution of source term
                tRes += aWStar * (
                        tPressureSpringProp->val()( 0 ) * tPressureFI->N_trans() * tPressureFI->val() );
            }

            // if PSPG stabilization
            if( tSPPSPG != nullptr )
            {
                // compute residual strong form for momentum
                Matrix< DDRMat > tRM;
                this->compute_residual_strong_form_momentum( tRM );

                // compute the residual weak form
                tRes += aWStar * ( trans( tPressureFI->dnNdxn( 1 ) ) * tSPPSPG->val()( 0 ) * tRM );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Pressure_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here pressure), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the pressure FIs
            Field_Interpolator * tPressureFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the source property
            const std::shared_ptr< Property > & tMassSourceProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::MASS_SOURCE ) );

            // get the pressure spring property
            const std::shared_ptr< Property > & tPressureSpringProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESSURE_SPRING ) );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tIncFluidCM =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density property
            const std::shared_ptr< Property > & tDensityProp = tIncFluidCM->get_property( "Density" );

            // get the incompressible flow stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPPSPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW ) );

            // get the density
            real tDensity = tDensityProp->val()( 0 );

            // if PSPG stabilization
            Matrix< DDRMat > tRM;
            if( tSPPSPG != nullptr )
            {
                // compute residual strong form for momentum
                this->compute_residual_strong_form_momentum( tRM );
            }

            // compute the Jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if dof type is velocity
                // FIXME protect velocity dof type
                if( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    Field_Interpolator * tVelocityFI = mLeaderFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                    // compute the Jacobian
                    tJac += aWStar * ( tPressureFI->N_trans() * tVelocityFI->div_operator() );
                }

                // if source term has dependency on the dof type
                if ( tMassSourceProp != nullptr )
                {
                    // if source property depends on dof type
                    if ( tMassSourceProp->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac -= aWStar * (
                                tPressureFI->N_trans() * tMassSourceProp->dPropdDOF( tDofType ) / tDensity );
                    }

                    // if density depends on dof type
                    if( tDensityProp->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac +=  aWStar * (
                                tPressureFI->N_trans() * tMassSourceProp->val() * tDensityProp->dPropdDOF( tDofType ) / std::pow( tDensity, 2 ) );
                    }
                }

                // if pressure spring property has dependency on the dof type
                if ( tPressureSpringProp != nullptr )
                {
                    // if dof type is pressure
                    if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // add contribution of pressure spring
                        tJac += aWStar * (
                                tPressureSpringProp->val()( 0 ) * tPressureFI->N_trans() * tPressureFI->N() );
                    }

                    // if pressure spring property depends on dof type
                    if ( tPressureSpringProp->check_dof_dependency( tDofType ) )
                    {
                        // add contribution of pressure spring
                        tJac += aWStar * (
                                tPressureFI->N_trans() * tPressureFI->val() * tPressureSpringProp->dPropdDOF( tDofType ) );
                    }
                }

                // if PSPG stabilization
                if( tSPPSPG != nullptr )
                {
                    // compute the Jacobian strong form for momentum
                    Matrix< DDRMat > tJM;
                    compute_jacobian_strong_form_momentum( tDofType, tJM );

                    // compute the Jacobian contribution from strong form
                    tJac += aWStar * ( trans( tPressureFI->dnNdxn( 1 ) ) * tSPPSPG->val()( 0 ) * tJM );

                    // if stabilization parameter has dependency on the dof type
                    if ( tSPPSPG->check_dof_dependency( tDofType ) )
                    {
                        // compute the Jacobian
                        tJac += aWStar * (
                                trans( tPressureFI->dnNdxn( 1 ) ) * tRM * tSPPSPG->dSPdLeaderDOF( tDofType ).get_row( 0 ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Bulk::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Bulk::compute_residual_strong_form_momentum(
                Matrix< DDRMat > & aRM )
        {
            // get the velocity and pressure FIs
            Field_Interpolator * tVelocityFI = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX ); //FIXME

            // get the density and gravity properties
            const std::shared_ptr< Property > & tGravityProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::GRAVITY ) );

            const std::shared_ptr< Property > & tThermalExpProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION ) );

            const std::shared_ptr< Property > & tRefTempProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::REF_TEMP ) );

            // get the inverted permeability (porosity) property
            const std::shared_ptr< Property > & tInvPermeabProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::INV_PERMEABILITY ) );

            // get the body load property
            const std::shared_ptr< Property > & tLoadProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tIncFluidCM =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density property from CM
            const std::shared_ptr< Property > & tDensityProp = tIncFluidCM->get_property( "Density" );

            // get the density value
            real tDensity = tDensityProp->val()( 0 );

            // compute the residual strong form
            aRM =   tDensity * trans( tVelocityFI->gradt( 1 ) ) +
                    tDensity * trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val() -
                    tIncFluidCM->divflux();

            // if body load
            if ( tLoadProp != nullptr )
            {
                // add contribution of body load to momentum residual
                aRM -= tLoadProp->val() ;
            }

            // if permeability
            if ( tInvPermeabProp != nullptr )
            {
                // add Brinkman term to residual strong form
                aRM += tInvPermeabProp->val()( 0 ) * tVelocityFI->val() ;
            }

            // if gravity
            if ( tGravityProp != nullptr )
            {
                // add gravity to residual strong form
                aRM += tDensity * tGravityProp->val();

                // if thermal expansion and reference temperature
                if( tThermalExpProp != nullptr && tRefTempProp != nullptr )
                {
                    // get the temperature field interpolator
                    // FIXME protect FI
                    Field_Interpolator * tTempFI =
                            mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                    // add contribution to residual
                    aRM -= tDensity * tGravityProp->val() * tThermalExpProp->val() *
                            ( tTempFI->val() - tRefTempProp->val() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian_strong_form_momentum(
                const Vector< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & aJM )
        {
            // get the velocity dof and the derivative dof FIs
            Field_Interpolator * tVelocityFI = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX ); //FIXME
            Field_Interpolator * tDerFI      = mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize aJMC
            aJM.set_size(
                    tVelocityFI->get_number_of_fields(),
                    tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the gravity properties
            const std::shared_ptr< Property > & tGravityProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::GRAVITY ) );

            const std::shared_ptr< Property > & tThermalExpProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION ) );

            const std::shared_ptr< Property > & tRefTempProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::REF_TEMP ) );

            // get the inverted permeability (porosity) property
            const std::shared_ptr< Property > & tInvPermeabProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::INV_PERMEABILITY ) );

            // get the body load property
            const std::shared_ptr< Property > & tLoadProp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tIncFluidCM =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density property from CM
            const std::shared_ptr< Property > & tDensityProp = tIncFluidCM->get_property( "Density" );

            // get the density value
            real tDensity = tDensityProp->val()( 0 );

            // if derivative wrt to velocity dof type
            // FIXME protect the velocity dof type
            if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
            {
                // compute the term uj vij
                Matrix< DDRMat > tujvij;
                this->compute_ujvij( tujvij );

                // compute the term dnNdtn
                Matrix< DDRMat > tdnNdtn;
                this->compute_dnNdtn( tdnNdtn );

                // compute the Jacobian strong form
                aJM =   tDensity * tdnNdtn
                        + tDensity * trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->N()
                        + tDensity * tujvij ;
            }
            else
            {
                aJM.fill( 0.0 );
            }

            // if density depends on dof type
            if( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution to Jacobian strong form
                aJM +=  trans( tVelocityFI->gradt( 1 ) ) * tDensityProp->dPropdDOF( aDofTypes ) +
                        trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val() * tDensityProp->dPropdDOF( aDofTypes );
            }

            // if CM depends on dof type
            if( tIncFluidCM->check_dof_dependency( aDofTypes ) )
            {
                // add contribution to Jacobian
                aJM -= tIncFluidCM->ddivfluxdu( aDofTypes );
            }

            // if permeability
            if ( tInvPermeabProp != nullptr )
            {
                // if derivative dof type is residual dof type
                if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
                {
                    // add Brinkman term to Jacobian strong form
                    aJM += tInvPermeabProp->val()( 0 ) * tVelocityFI->N();
                }

                // if permeability depends on dof type
                if( tInvPermeabProp->check_dof_dependency( aDofTypes ) )
                {
                    // add Brinkman term to Jacobian strong form
                    aJM += tVelocityFI->val() * tInvPermeabProp->dPropdDOF( aDofTypes );
                }
            }

            // if body load
            if ( tLoadProp != nullptr )
            {
                // if indirect DoF dependency of body load
                if( tLoadProp->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to momentum Jacobian
                    aJM -= tLoadProp->dPropdDOF( aDofTypes );
                }
            }

            // if gravity
            if ( tGravityProp != nullptr )
            {
                // if gravity depends on dof type
                if( tGravityProp->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to Jacobian
                    aJM += tDensity * tGravityProp->dPropdDOF( aDofTypes );
                }

                // if density depends on dof type
                if( tDensityProp->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to Jacobian
                    aJM += tGravityProp->val() * tDensityProp->dPropdDOF( aDofTypes );
                }

                // if thermal expansion and reference temperature
                if( tThermalExpProp != nullptr && tRefTempProp != nullptr )
                {
                    // get the temperature field interpolator
                    // FIXME protect FI
                    Field_Interpolator * tTempFI =
                            mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                    // if dof type is temperature
                    if( aDofTypes( 0 ) == MSI::Dof_Type::TEMP )
                    {
                        // add contribution to Jacobian
                        aJM -= tDensity * tGravityProp->val() * tThermalExpProp->val() * tTempFI->N();
                    }

                    // if thermal expansion property depends on dof type
                    if( tThermalExpProp->check_dof_dependency( aDofTypes ) )
                    {
                        // add contribution to Jacobian
                        aJM -=  tDensity * tGravityProp->val() *
                                ( tTempFI->val() - tRefTempProp->val() ) * tThermalExpProp->dPropdDOF( aDofTypes );
                    }

                    // if reference temperature property depends on dof type
                    if( tRefTempProp->check_dof_dependency( aDofTypes ) )
                    {
                        // add contribution to Jacobian
                        aJM +=  tDensity * tGravityProp->val() *
                                tThermalExpProp->val() * tRefTempProp->dPropdDOF( aDofTypes );
                    }

                    // if gravity property has dependency on the dof type
                    if ( tGravityProp->check_dof_dependency( aDofTypes ) )
                    {
                        // compute the Jacobian
                        aJM -=  tDensity *
                                tThermalExpProp->val()( 0 ) * ( tTempFI->val()( 0 ) - tRefTempProp->val()( 0 ) ) *
                                tGravityProp->dPropdDOF( aDofTypes );
                    }

                    // if density depends on dof type
                    if( tDensityProp->check_dof_dependency( aDofTypes ) )
                    {
                        // add density contribution to residual strong form
                        aJM -=  tGravityProp->val() *
                                tThermalExpProp->val() * ( tTempFI->val() - tRefTempProp->val() ) *
                                tDensityProp->dPropdDOF( aDofTypes );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Bulk::compute_ujvij(
                Matrix< DDRMat > & aujvij )
        {
            // get the velocity dof type FI
            // FIXME protect velocity dof type
            Field_Interpolator * tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // init size for uj vij
            uint tNumRow = tVelocityFI->dnNdxn( 1 ).n_rows();
            uint tNumCol = tVelocityFI->dnNdxn( 1 ).n_cols();
            aujvij.set_size( tNumRow, tNumRow * tNumCol, 0.0 );

            // loop over the number of rows
            for( uint i = 0; i < tNumRow; i++ )
            {
                // compute uj vij
                aujvij( { i, i }, { i * tNumCol, ( i + 1 ) * tNumCol -1 } ) =
                        trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Bulk::compute_ujvijrm(
                Matrix< DDRMat > & aujvijrm,
                Matrix< DDRMat > & arm )
        {
            // get the velocity dof type FI
            // FIXME protect velocity dof type
            Field_Interpolator * tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // set size for uj vij rM
            uint tNumField = tVelocityFI->get_number_of_fields();
            uint tNumBases = tVelocityFI->get_number_of_space_time_bases();
            aujvijrm.set_size( tNumField * tNumBases, tNumField * tNumBases , 0.0 );

            // loop over the number of fields
            for( uint iField = 0; iField < tNumField; iField++ )
            {
                // loop over the number of fields
                for( uint iField2 = 0; iField2 < tNumField; iField2++ )
                {
                    // compute uj vij rm
                    aujvijrm(
                            { iField  * tNumBases, ( iField + 1 )  * tNumBases - 1 },
                            { iField2 * tNumBases, ( iField2 + 1 ) * tNumBases - 1 } ) =
                                    trans( tVelocityFI->dnNdxn( 1 ).get_row( iField2 ) ) * tVelocityFI->NBuild() * arm( iField );
                }
            }
        }

        //------------------------------------------------------------------------------

        // FIXME provided directly by the field interpolator?
        void IWG_Incompressible_NS_Pressure_Bulk::compute_dnNdtn( Matrix< DDRMat > & adnNdtn )
        {
            // get the velocity dof type FI
            Field_Interpolator * tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // init size for dnNdtn
            uint tNumRowt = tVelocityFI->get_number_of_fields();
            uint tNumColt = tVelocityFI->dnNdtn( 1 ).n_cols();
            adnNdtn.set_size( tNumRowt, tNumRowt * tNumColt , 0.0 );

            // loop over the fields
            for( uint iField = 0; iField < tNumRowt; iField++ )
            {
                // fill the matrix for each dimension
                adnNdtn( { iField, iField }, { iField * tNumColt, ( iField + 1 ) * tNumColt - 1 } ) =
                        tVelocityFI->dnNdtn( 1 ).matrix_data();
            }
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
