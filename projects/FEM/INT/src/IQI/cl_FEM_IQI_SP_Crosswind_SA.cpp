/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Strong_Residual_SA.cpp
 *
 */

#include "cl_FEM_IQI_SP_Crosswind_SA.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_CM_Spalart_Allmaras_Turbulence.hpp"

#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_SP_Crosswind_SA::IQI_SP_Crosswind_SA()
    {
        // set fem IQI type
        mFEMIQIType = fem::IQI_Type::SP_CROSSWIND_SA;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "SpalartAllmarasTurbulence" ] = static_cast< uint >( IQI_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IQI_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "Crosswind" ] = static_cast< uint >( IQI_Stabilization_Type::CROSSWIND );
    }

    //------------------------------------------------------------------------------

    void
    IQI_SP_Crosswind_SA::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get the residual viscosity FI
        Field_Interpolator* tFIViscosity =    //
                mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VISCOSITY );

        // get the SA turbulence CM
        const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

        // cast constitutive model base class pointer to SA constitutive model
        CM_Spalart_Allmaras_Turbulence* tCMSATurbulencePtr =
                dynamic_cast< CM_Spalart_Allmaras_Turbulence* >( tCMSATurbulence.get() );

        // get the crosswind stabilization parameter
        const std::shared_ptr< Stabilization_Parameter >& tSPCrosswind =
                mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::CROSSWIND ) );

        // get zero tolerance for crosswind
        const real tEpsilon = tSPCrosswind->val()( 1 );

        // compute residual of the strong form
        Matrix< DDRMat > tR;
        this->compute_residual_strong_form( tR );

        // compute the norm of the viscosity gradient
        real tNorm = std::max( norm( tFIViscosity->gradx( 1 ) ), tEpsilon );

        // compute the abs of the strong form of the residual
        real tRAbs = std::max( std::abs( tR( 0, 0 ) ), tEpsilon );

        // compute full crosswind stabilization parameter value
        aQI = std::max( tSPCrosswind->val()( 0 ) * tRAbs / tNorm - tCMSATurbulencePtr->diffusion_coefficient()( 0 ), 0.0 );
    }

    //------------------------------------------------------------------------------

    void
    IQI_SP_Crosswind_SA::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // evaluate strong form
        Matrix< DDRMat > tQI( 1, 1 );
        this->compute_QI( tQI );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * tQI;
    }

    //------------------------------------------------------------------------------

    void
    IQI_SP_Crosswind_SA::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false,
                "IQI_SP_Crosswind_SA::compute_dQIdu - not implemented\n." );
    }

    //------------------------------------------------------------------------------

    void
    IQI_SP_Crosswind_SA::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false,
                "IQI_SP_Crosswind_SA::compute_dQIdu - not implemented\n." );
    }

    //------------------------------------------------------------------------------

    void
    IQI_SP_Crosswind_SA::compute_residual_strong_form( Matrix< DDRMat >& aR )
    {
        // get the residual viscosity FI
        // FIXME protect dof type
        Field_Interpolator* tFIViscosity =    //
                mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VISCOSITY );

        // get the SA turbulence CM
        const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

        // compute strong form of residual
        aR = tFIViscosity->gradt( 1 )
           + trans( tCMSATurbulence->modified_velocity() ) * tFIViscosity->gradx( 1 )
           - tCMSATurbulence->production_term()
           + tCMSATurbulence->wall_destruction_term()
           - tCMSATurbulence->divflux();

        MORIS_ASSERT( isfinite( aR ),
                "IQI_SP_Crosswind_SA::compute_residual_strong_form - Residual contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
