/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_SUPG_Advection.cpp
 *
 */

#include "cl_FEM_SP_SUPG_Advection.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
// LINALG/src
#include "fn_norm.hpp"
#include "fn_dot.hpp"

#include "fn_FEM_CM_Phase_State_Functions.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_SUPG_Advection::SP_SUPG_Advection()
        {
            // set the property pointer cell size
            mLeaderProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Conductivity" ]       = static_cast< uint >( Property_Type::CONDUCTIVITY );
            mPropertyMap[ "Density" ]            = static_cast< uint >( Property_Type::DENSITY );
            mPropertyMap[ "HeatCapacity" ]       = static_cast< uint >( Property_Type::HEAT_CAPACITY );
            mPropertyMap[ "LatentHeat" ]         = static_cast< uint >( Property_Type::LATENT_HEAT );
            mPropertyMap[ "PCTemp" ]             = static_cast< uint >( Property_Type::PC_TEMP );
            mPropertyMap[ "PhaseStateFunction" ] = static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION );
            mPropertyMap[ "PhaseChangeConst" ]   = static_cast< uint >( Property_Type::PHASE_CHANGE_CONST );
            mPropertyMap[ "Source" ]             = static_cast< uint >( Property_Type::SOURCE );
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Advection::set_parameters( moris::Vector< Matrix< DDRMat > > aParameters )
        {
            // FIXME not necessary
            // set mParameters
            mParameters = aParameters;

            // get number of parameters
            uint tParamSize = aParameters.size();

            // check for proper size of constant function parameters
            MORIS_ERROR( tParamSize <= 1,
                    "SP_SUPG_Advection::set_parameters - no more than one constant parameter can be set." );

            // if a parameter is specified
            if ( tParamSize > 0 )
            {
                // check for proper parameter type; here just a scalar
                MORIS_ERROR( aParameters( 0 ).numel() == 1,
                        "SP_SUPG_Advection::set_parameters - 1st parameter is not a scalar but a vector." );

                mBetaTime = aParameters( 0 )( 0 );

                // set beta time flag to true
                if ( std::abs( mBetaTime ) > MORIS_REAL_EPS )
                {
                    mSetBetaTime = true;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Advection::set_dof_type_list(
                moris::Vector< moris::Vector< MSI::Dof_Type > >& aDofTypes,
                moris::Vector< std::string >&                  aDofStrings,
                mtk::Leader_Follower                            aIsLeader )
        {
            // switch on leader follower
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    // set dof type list
                    mLeaderDofTypes = aDofTypes;

                    // loop on dof type
                    for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if ( tDofString == "Velocity" )
                        {
                            mLeaderDofVelocity = tDofType;
                        }
                        else if ( tDofString == "ScalarField" )
                        {
                            mLeaderDofScalarField = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false,
                                    "SP_SUPG_Advection::set_dof_type_list - Unknown aDofString : %s \n",
                                    tDofString.c_str() );
                        }
                    }
                    break;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // set dof type list
                    mFollowerDofTypes = aDofTypes;
                    break;
                }
                default:
                    MORIS_ERROR( false, "SP_SUPG_Advection::set_dof_type_list - unknown leader follower type." );
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Advection::build_global_dof_type_list()
        {
            // call parent implementation
            Stabilization_Parameter::build_global_dof_type_list();

            // get number of dof types
            uint tNumLeaderGlobalDofTypes = mLeaderGlobalDofTypes.size();

            // init child specific eval flags
            mdLengthScaledLeaderDofEval.set_size( tNumLeaderGlobalDofTypes, 1, true );

            // init child specific storage
            mdLengthScaledLeaderDof.resize( tNumLeaderGlobalDofTypes );
        }

        //------------------------------------------------------------------------------
        /**
         * reset evaluation flags
         */
        void
        SP_SUPG_Advection::reset_eval_flags()
        {
            // call parent implementation
            Stabilization_Parameter::reset_eval_flags();

            // reset child specific eval flags for chi
            mLengthScaleEval = true;
            mdLengthScaledLeaderDofEval.fill( true );
        }

        //------------------------------------------------------------------------------

        real
        SP_SUPG_Advection::compute_effective_conductivity()
        {
            // get the conductivity property
            const std::shared_ptr< Property >& tPropConductivity =
                    mLeaderProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // get the density property
            const std::shared_ptr< Property >& tPropDensity =
                    mLeaderProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the capacity property
            const std::shared_ptr< Property >& tPropCapacity =
                    mLeaderProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            // get the latent heat property
            const std::shared_ptr< Property >& tPropLatentHeat =
                    mLeaderProp( static_cast< uint >( Property_Type::LATENT_HEAT ) );

            // check that conductivity is set
            MORIS_ASSERT( tPropConductivity != nullptr,
                    "SP_SUPG_Advection::compute_effective_conductivity - conductivity not defined\n" );

            // get contribution of density to effective conductivity
            real tDensity = 1.0;
            if ( tPropDensity != nullptr )
            {
                tDensity = tPropDensity->val()( 0 );
            }

            // get contribution of capacity to effective conductivity
            real tCapacity = 1.0;
            if ( tPropCapacity != nullptr )
            {
                tCapacity = tPropCapacity->val()( 0 );
            }

            // get contribution of latent to effective conductivity
            real tLatentHeatContrib = 0.0;
            if ( tPropLatentHeat != nullptr )
            {
                // get the phase change properties
                const std::shared_ptr< Property >& tPropPCTemp =
                        mLeaderProp( static_cast< uint >( Property_Type::PC_TEMP ) );

                const std::shared_ptr< Property >& tPropPhaseChangeFunction =
                        mLeaderProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );

                const std::shared_ptr< Property >& tPropPhaseChangeConstant =
                        mLeaderProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );

                // check that all phase change properties are set
                MORIS_ASSERT( tPropPCTemp != nullptr and tPropPhaseChangeFunction != nullptr and tPropPhaseChangeConstant != nullptr,
                        "SP_SUPG_Advection::compute_effective_conductivity - some or all change properties are not defined\n" );

                // get the scalar field FI
                Field_Interpolator* tFIScalarField =
                        mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofScalarField );

                // compute derivative of Phase State Function
                real tdfdT = eval_dFdTemp(
                        tPropPCTemp->val()( 0 ),
                        tPropPhaseChangeConstant->val()( 0 ),
                        tPropPhaseChangeFunction->val()( 0 ),
                        tFIScalarField );

                // compute contribution of latent heat
                tLatentHeatContrib = tPropLatentHeat->val()( 0 ) * tdfdT;
            }

            // compute effective conductivity
            return tPropConductivity->val()( 0 ) / ( tDensity * ( tCapacity + tLatentHeatContrib ) );
        }

        //------------------------------------------------------------------------------

        bool
        SP_SUPG_Advection::compute_derivative_of_effective_conductivity(
                Matrix< DDRMat >&                   aEffectiveConductivitydu,
                const moris::Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the conductivity property
            const std::shared_ptr< Property >& tPropConductivity =
                    mLeaderProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // get the density property
            const std::shared_ptr< Property >& tPropDensity =
                    mLeaderProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the capacity property
            const std::shared_ptr< Property >& tPropCapacity =
                    mLeaderProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            // get the latent heat property
            const std::shared_ptr< Property >& tPropLatentHeat =
                    mLeaderProp( static_cast< uint >( Property_Type::LATENT_HEAT ) );

            // get conductivity
            real tConductivity = tPropConductivity->val()( 0 );

            // get contribution of density to effective conductivity
            real tDensity = 1.0;
            if ( tPropDensity != nullptr )
            {
                tDensity = tPropDensity->val()( 0 );
            }

            // get contribution of capacity to effective conductivity
            real tCapacity = 1.0;
            if ( tPropCapacity != nullptr )
            {
                tCapacity = tPropCapacity->val()( 0 );
            }

            // get contribution of latent to effective conductivity
            real tLatentHeatContrib = 0.0;
            if ( tPropLatentHeat != nullptr )
            {
                // get the phase change properties
                const std::shared_ptr< Property >& tPropPCTemp =
                        mLeaderProp( static_cast< uint >( Property_Type::PC_TEMP ) );

                const std::shared_ptr< Property >& tPropPhaseChangeFunction =
                        mLeaderProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );

                const std::shared_ptr< Property >& tPropPhaseChangeConstant =
                        mLeaderProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );

                // get the temperature FI
                Field_Interpolator* tFIScalarField =
                        mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofScalarField );

                // compute derivative of Phase State Function
                real tdfdT = eval_dFdTemp(
                        tPropPCTemp->val()( 0 ),
                        tPropPhaseChangeConstant->val()( 0 ),
                        tPropPhaseChangeFunction->val()( 0 ),
                        tFIScalarField );

                // add contribution of latent heat
                tLatentHeatContrib = tPropLatentHeat->val()( 0 ) * tdfdT;
            }

            // set flag for dependency to false
            bool tIsDependent = false;

            // consider dependency of conductivity on dof types
            if ( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                aEffectiveConductivitydu = 1.0 / ( tDensity * ( tCapacity + tLatentHeatContrib ) ) * tPropConductivity->dPropdDOF( aDofTypes );
                tIsDependent             = true;
            }
            else
            {
                aEffectiveConductivitydu.fill( 0.0 );
            }

            // consider dependency of density on dof types
            if ( tPropDensity != nullptr )
            {
                if ( tPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    const real tFactor = tConductivity * ( tCapacity + tLatentHeatContrib ) / std::pow( tDensity * ( tCapacity + tLatentHeatContrib ), 2.0 );

                    aEffectiveConductivitydu -= tFactor * tPropDensity->dPropdDOF( aDofTypes );

                    tIsDependent = true;
                }
            }

            // consider dependency of density on dof types
            if ( tPropCapacity != nullptr )
            {
                if ( tPropCapacity->check_dof_dependency( aDofTypes ) )
                {
                    const real tFactor = tConductivity * tDensity / std::pow( tDensity * ( tCapacity + tLatentHeatContrib ), 2.0 );

                    aEffectiveConductivitydu -= tFactor * tPropCapacity->dPropdDOF( aDofTypes );

                    tIsDependent = true;
                }
            }

            // consider dependency of latent heat contribution on dof types
            if ( tPropLatentHeat != nullptr )
            {
                // get the phase change properties
                const std::shared_ptr< Property >& tPropPCTemp =
                        mLeaderProp( static_cast< uint >( Property_Type::PC_TEMP ) );

                const std::shared_ptr< Property >& tPropPhaseChangeFunction =
                        mLeaderProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );

                const std::shared_ptr< Property >& tPropPhaseChangeConstant =
                        mLeaderProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );

                // consider dependency of phase state function on scalar field
                if ( aDofTypes( 0 ) == mLeaderDofScalarField )
                {
                    // get the scalar field FI
                    Field_Interpolator* tFIScalarField =
                            mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofScalarField );

                    const moris::Matrix< DDRMat > dfdDof = eval_dFdTempdDOF(
                            tPropPCTemp->val()( 0 ),
                            tPropPhaseChangeConstant->val()( 0 ),
                            tPropPhaseChangeFunction->val()( 0 ),
                            tFIScalarField );

                    const real tFactor = tConductivity * tDensity * tPropLatentHeat->val()( 0 ) / std::pow( tDensity * ( tCapacity + tLatentHeatContrib ), 2.0 );

                    aEffectiveConductivitydu -= tFactor * dfdDof;

                    tIsDependent = true;
                }

                // if density depends on dof type
                if ( tPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    MORIS_ERROR( false, "SP_SUPG_Advection::compute_derivative_of_effective_conductivity - %s\n", "Dof dependence of density for phase change not implemented.\n" );
                }

                // if latent heat depends on the dof type
                if ( tPropLatentHeat->check_dof_dependency( aDofTypes ) )
                {
                    MORIS_ERROR( false, "SP_SUPG_Advection::compute_derivative_of_effective_conductivity - %s\n", "Dof dependence of latent heat not implemented.\n" );
                }
            }

            // return flag whether dof dependence exists
            return tIsDependent;
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Advection::eval_SP()
        {
            // get the velocity FI
            Field_Interpolator* tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

            // get the mass source
            const std::shared_ptr< Property >& tSourceProp =
                    mLeaderProp( static_cast< uint >( Property_Type::SOURCE ) );

            // compute effective conductivity
            const real tEffectiveConductivity = this->compute_effective_conductivity();

            // compute and threshold the velocity norm (thresholding for consistency with derivatives)
            const real tNorm = std::max( norm( tVelocityFI->val() ), mEpsilon );

            // compute and threshold hugn
            const real tHugn = this->length_scale();

            // compute tau1
            const real tTau1 = 2.0 * tNorm / tHugn;

            // compute tau2
            const real tTau2 = 4.0 * tEffectiveConductivity / std::pow( tHugn, 2.0 );

            // compute sum of square terms
            real tSum = std::pow( tTau1, 2.0 ) + std::pow( tTau2, 2.0 );

            // if time solve
            if ( mSetBetaTime )
            {
                // compute time increment tDeltaT
                const real tDeltaT = mBetaTime * mLeaderFIManager->get_IP_geometry_interpolator()->get_time_step();

                // compute tau3
                const real tTau3 = 2.0 / tDeltaT;

                // add contribution from time term
                tSum += std::pow( tTau3, 2.0 );
            }

            // add contribution from source term
            if ( tSourceProp != nullptr )
            {
                tSum += std::pow( tSourceProp->val()( 0 ), 2.0 );
            }

            // threshold sum of square terms
            tSum = std::max( tSum, mEpsilon );

            // compute stabilization parameter value
            mPPVal = { { std::pow( tSum, -0.5 ) } };
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Advection::eval_dSPdLeaderDOF(
                const moris::Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            const uint tDofIndex = mLeaderGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator* tFIDer =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdLeaderDof, dTau1dDof, dTau3dDof
            mdPPdLeaderDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            Matrix< DDRMat > tdTau1dDof( 1, tFIDer->get_number_of_space_time_coefficients() );
            Matrix< DDRMat > tdTau2dDof( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the velocity FI
            Field_Interpolator* tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

            // get the source property
            const std::shared_ptr< Property >& tSourceProp =
                    mLeaderProp( static_cast< uint >( Property_Type::SOURCE ) );

            // compute effective conductivity
            const real tEffectiveConductivity = this->compute_effective_conductivity();

            // compute and threshold the velocity norm (thresholding for consistency with derivatives)
            const real tNorm = std::max( norm( tVelocityFI->val() ), mEpsilon );

            // compute and threshold hugn
            const real tHugn = this->length_scale();

            // compute tau1
            const real tTau1 = 2.0 * tNorm / tHugn;

            // compute tau2
            const real tTau2 = 4.0 * tEffectiveConductivity / std::pow( tHugn, 2.0 );

            // compute sum of square terms
            real tSum = std::pow( tTau1, 2.0 ) + std::pow( tTau2, 2.0 );

            // if time solve
            if ( mSetBetaTime )
            {
                // compute time increment tDeltaT
                const real tDeltaT = mBetaTime * mLeaderFIManager->get_IP_geometry_interpolator()->get_time_step();

                // compute tau3
                const real tTau3 = 2.0 / tDeltaT;

                // add contribution from time term
                tSum += std::pow( tTau3, 2.0 );
            }

            // add contribution from source term
            if ( tSourceProp != nullptr )
            {
                tSum += std::pow( tSourceProp->val()( 0 ), 2.0 );
            }

            // compute dSPdu
            if ( tSum > mEpsilon )
            {
                // if dof type is velocity
                if ( aDofTypes( 0 ) == mLeaderDofVelocity )
                {
                    // compute derivative of hugn wrt velocity dof
                    Matrix< DDRMat > tdNormdu( 1, tVelocityFI->get_number_of_space_time_coefficients() );
                    Matrix< DDRMat > tdAbsdu( 1, tVelocityFI->get_number_of_space_time_coefficients(), 0.0 );

                    // compute derivative of the velocity norm (compute only derivative if not thresholded)
                    if ( tNorm > mEpsilon )
                    {
                        tdNormdu = trans( tVelocityFI->val() ) * tVelocityFI->N() / tNorm;
                    }
                    else
                    {
                        tdNormdu.fill( 0.0 );
                    }

                    // compute dtau1du
                    tdTau1dDof = 2.0 * ( tHugn * tdNormdu - this->dlengthscaledleaderu( aDofTypes ) * tNorm ) / std::pow( tHugn, 2.0 );

                    // compute dtau2du
                    tdTau2dDof = -8.0 * tEffectiveConductivity * this->dlengthscaledleaderu( aDofTypes ) / std::pow( tHugn, 3.0 );
                }
                else
                {
                    tdTau1dDof.fill( 0.0 );
                    tdTau2dDof.fill( 0.0 );
                }

                // consider derivative of effective conductivity on dof types
                Matrix< DDRMat > tEffectiveConductivitydu( 1, tFIDer->get_number_of_space_time_coefficients() );

                bool tIsDependent = this->compute_derivative_of_effective_conductivity(
                        tEffectiveConductivitydu,
                        aDofTypes );

                // compute dtau2du
                if ( tIsDependent )
                {
                    tdTau2dDof += 4.0 * tEffectiveConductivitydu / std::pow( tHugn, 2.0 );
                }

                const real tPrefactor = -std::pow( tSum, -1.5 );

                mdPPdLeaderDof( tDofIndex ) = tPrefactor * ( tTau1 * tdTau1dDof + tTau2 * tdTau2dDof );

                if ( tSourceProp != nullptr )
                {
                    if ( tSourceProp->check_dof_dependency( aDofTypes ) )
                    {
                        // compute dtau3du
                        mdPPdLeaderDof( tDofIndex ) += tPrefactor * tSourceProp->val()( 0 ) * tSourceProp->dPropdDOF( aDofTypes );
                    }
                }
            }
            else
            {
                mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mdPPdLeaderDof( tDofIndex ) ),
                    "SP_SUPG_Advection::eval_dSPdLeaderDOF - mdPPdLeaderDof contains NAN or INF, exiting for tDofIndex = %d !\n",
                    tDofIndex );
        }

        //------------------------------------------------------------------------------

        real
        SP_SUPG_Advection::length_scale()
        {
            // if the length scale parameter was not evaluated
            if ( mLengthScaleEval )
            {
                // evaluate the length scale parameter
                this->eval_length_scale();

                // set bool for evaluation
                mLengthScaleEval = false;
            }
            // return the length scale parameter value
            return mLengthScale;
        }

        void
        SP_SUPG_Advection::eval_length_scale()
        {
            // get the velocity FI
            Field_Interpolator* tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

            // compute and threshold the velocity norm (thresholding for consistency with derivatives)
            const real tNorm = std::max( norm( tVelocityFI->val() ), mEpsilon );

            // get the abs term
            const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
            real       tAbs      = 0.0;
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                tAbs += std::abs( dot( tVelocityFI->val(), tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) );
            }

            // threshold tAbs
            tAbs = std::max( tAbs, mEpsilon );

            // compute and threshold hugn
            mLengthScale = std::max( 2.0 * tNorm / tAbs, mEpsilon );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        SP_SUPG_Advection::dlengthscaledleaderu(
                const moris::Vector< MSI::Dof_Type >& aDofType )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType, mtk::Leader_Follower::LEADER ),
                    "SP_SUPG_Advection::dlengthscaledleaderu - no dependency on this dof type." );

            // get the dof index
            uint tDofIndex = mLeaderGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdLengthScaledLeaderDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dlengthscaledleaderu( aDofType );

                // set bool for evaluation
                mdLengthScaledLeaderDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdLengthScaledLeaderDof( tDofIndex );
        }

        void
        SP_SUPG_Advection::eval_dlengthscaledleaderu(
                const moris::Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mLeaderGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator* tFIDer =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdLengthScaledLeaderDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the length scale value
            real tHugn = this->length_scale();

            // compute derivative of hugn (compute only derivative if not thresholded)
            if ( tHugn > mEpsilon )
            {
                // get the velocity FI
                Field_Interpolator* tVelocityFI =
                        mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

                // compute and threshold the velocity norm (thresholding for consistency with derivatives)
                const real tNorm = std::max( norm( tVelocityFI->val() ), mEpsilon );

                // compute derivative of the modified velocity norm (compute only derivative if not thresholded)
                Matrix< DDRMat > tdNormdu( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
                if ( tNorm > mEpsilon )
                {
                    if ( aDofTypes( 0 ) == mLeaderDofVelocity )
                    {
                        tdNormdu += trans( tVelocityFI->val() ) * tVelocityFI->N() / tNorm;
                    }
                }

                // get the abs term
                const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
                real       tAbs      = 0.0;
                for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                {
                    tAbs += std::abs( dot( tVelocityFI->val(), tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) );
                }

                // threshold tAbs
                tAbs = std::max( tAbs, mEpsilon );

                // compute derivative of the abs term (compute only derivative if not thresholded)
                Matrix< DDRMat > tdAbsdu( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
                if ( tAbs > mEpsilon )
                {
                    if ( aDofTypes( 0 ) == mLeaderDofVelocity )
                    {
                        uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
                        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                        {
                            real tAdd = dot( tVelocityFI->val(), tVelocityFI->dnNdxn( 1 ).get_column( iNode ) );

                            // handle case that tAdd( 0, 0 ) is smaller than threshold
                            if ( std::abs( tAdd ) > mEpsilon )
                            {
                                tdAbsdu +=
                                        tAdd * trans( tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) * tVelocityFI->N() / std::abs( tAdd );
                            }
                        }
                    }
                }

                // compute the derivative of the length scale
                mdLengthScaledLeaderDof( tDofIndex ) =
                        2.0 * ( tdNormdu * tAbs - tdAbsdu * tNorm ) / std::pow( tAbs, 2.0 );
            }
            else
            {
                mdLengthScaledLeaderDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

