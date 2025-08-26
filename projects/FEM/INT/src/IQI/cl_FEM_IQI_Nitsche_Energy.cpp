/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Traction.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Nitsche_Energy.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Nitsche_Energy::IQI_Nitsche_Energy()
    {
        mFEMIQIType = fem::IQI_Type::NITSCHE_ENERGY;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // Set size for property model pointer cell
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // Set size for stabilization param pointer cell
        mStabilizationParam.resize( static_cast< uint >( IQI_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map, property map and stabilization map
        mConstitutiveMap[ "TractionCM" ] = (uint)IQI_Constitutive_Type::TRACTION_CM;
        mPropertyMap["Dirichlet"] = (uint)IQI_Property_Type::DIRICHLET;
        mStabilizationMap["DirichletNitsche"] = (uint)IQI_Stabilization_Type::DIRICHLET_NITSCHE;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Nitsche_Energy::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get the constitutive model for computing the fluxes/tractions
        const std::shared_ptr< Constitutive_Model >& tCMLeader =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );
        
        // get the specified displacement        
        const std::shared_ptr< Property >& tPropDirichlet =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::DIRICHLET ) );

        // get the stabilization parameter
        const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::DIRICHLET_NITSCHE ) );

        // Obtain leader and follower field interpolators to get the displacement
        Field_Interpolator* tFILeader   = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
        //Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

        MORIS_ASSERT( mNormal.numel() > 0,
                "IQI_Nitsche_Energy::compute_QI() - "
                "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                "and averaged such that there is a well-defined normal." );

        // evaluate traction
        const Matrix< DDRMat >
                tTraction = tCMLeader->traction( mNormal );

        // Get Displacement-
        const Matrix< DDRMat > tDisplacement = tFILeader->val();

        // Compute the integrand
        aQI = trans( tTraction ) * ( tDisplacement - tPropDirichlet->val()) + tSPNitsche->val() * (trans( tDisplacement - tPropDirichlet->val()) * ( tDisplacement - tPropDirichlet->val()));

    }

    //------------------------------------------------------------------------------

    void
    IQI_Nitsche_Energy::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the constitutive model for computing the fluxes/tractions
        const std::shared_ptr< Constitutive_Model >& tCMLeader =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );
        
        // get the specified displacement        
        const std::shared_ptr< Property >& tPropDirichlet =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::DIRICHLET ) );

        // get the stabilization parameter
        const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::DIRICHLET_NITSCHE ) );
        
        // Obtain leader and follower field interpolators to get the displacement
        Field_Interpolator* tFILeader   = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
        //Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

        MORIS_ASSERT( mNormal.numel() > 0,
                "IQI_Nitsche_Energy::compute_QI() - "
                "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                "and averaged such that there is a well-defined normal." );

        // evaluate traction
        const Matrix< DDRMat >
                tTraction = tCMLeader->traction( mNormal );

        // Get Displacement-
        const Matrix< DDRMat > tDisplacement = tFILeader->val();

        // Compute the integrand
        Matrix < DDRMat > tQI = tTraction * ( tDisplacement - tPropDirichlet->val()) + tSPNitsche->val() * (trans( tDisplacement - tPropDirichlet->val()) * ( tDisplacement - tPropDirichlet->val()));

        // add the contribution
        mSet->get_QI()( tQIIndex ) += aWStar * tQI;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Nitsche_Energy::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Nitsche_Energy::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
