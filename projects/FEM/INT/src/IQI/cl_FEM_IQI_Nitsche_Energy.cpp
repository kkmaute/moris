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
#include "fn_eye.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Nitsche_Energy::IQI_Nitsche_Energy( 
            enum CM_Function_Type aStressType,
            enum CM_Function_Type aStrainType,
            sint                  aNitscheType 
        )
    {
        mFEMIQIType = fem::IQI_Type::NITSCHE_ENERGY;

        mStressType = aStressType;
        mStrainType = aStrainType;
        mNitscheType = aNitscheType;

        // If Nitsche type is -1 then throw an error
        MORIS_ERROR( mNitscheType != -1, "IQI_Nitsche_Energy::IQI_Nitsche_Energy - Asymmetric Nitsche does not have an associated energy" );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // Set size for property model pointer cell
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // Set size for stabilization param pointer cell
        mStabilizationParam.resize( static_cast< uint >( IQI_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map, property map and stabilization map
        mConstitutiveMap[ "TractionCM" ] = (uint)IQI_Constitutive_Type::TRACTION_CM;
        mPropertyMap["Dirichlet"] = (uint)IQI_Property_Type::DIRICHLET;
        mPropertyMap["Select"] = (uint)IQI_Property_Type::SELECT;
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
        
        // Get the projection matrix property
        const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) );

        // get the stabilization parameter
        const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::DIRICHLET_NITSCHE ) );

        // Obtain leader and follower field interpolators to get the displacement
        Field_Interpolator* tFILeader   = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
        //Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

        Matrix< DDRMat > tM;
        if ( tPropSelect == nullptr )
        {
            // get number of fields which should equal spatial dimension
            const uint tSpaceDim = tFILeader->get_number_of_fields();

            // set selection matrix as identity
            eye( tSpaceDim, tSpaceDim, tM );
        }
        else
        {
            tM = tPropSelect->val();

            // skip computing residual if projection matrix is zero
            if ( norm( tM ) < MORIS_REAL_EPS )
            {
                aQI.set_size( 1, 1, 0.0 );
                return;
            }
        }

        const Matrix< DDRMat > tJump = tFILeader->val() - tPropDirichlet->val();

        aQI = 0.5 * tSPNitsche->val() * ( trans( tJump ) * tM * tJump );

        // Pure penalty does not require traction; zero times an invalid traction is NaN.
        if ( mNitscheType != 0 )
        {
            MORIS_ASSERT( mNormal.numel() > 0,
                    "IQI_Nitsche_Energy::compute_QI() - "
                    "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                    "and averaged such that there is a well-defined normal." );

            const Matrix< DDRMat > tTraction = tCMLeader->traction( mNormal, mStressType );
            aQI += mNitscheType * trans( tTraction ) * tM * tJump;
        }

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
        
        // Get the projection matrix property
        const std::shared_ptr< Property >& tPropSelect =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) );

        // get the stabilization parameter
        const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::DIRICHLET_NITSCHE ) );
        
        // Obtain leader and follower field interpolators to get the displacement
        Field_Interpolator* tFILeader   = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
        //Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
        
        Matrix< DDRMat > tM;
        if ( tPropSelect == nullptr )
        {
            // get number of fields which should equal spatial dimension
            const uint tSpaceDim = tFILeader->get_number_of_fields();

            // set selection matrix as identity
            eye( tSpaceDim, tSpaceDim, tM );
        }
        else
        {
            tM = tPropSelect->val();

            // skip computing residual if projection matrix is zero
            if ( norm( tM ) < MORIS_REAL_EPS ) return;
        }

        const Matrix< DDRMat > tJump = tFILeader->val() - tPropDirichlet->val();

        Matrix< DDRMat > tQI = 0.5 * tSPNitsche->val() * ( trans( tJump ) * tM * tJump );

        // Evaluate the constitutive traction only when its contribution is active.
        if ( mNitscheType != 0 )
        {
            MORIS_ASSERT( mNormal.numel() > 0,
                    "IQI_Nitsche_Energy::compute_QI() - "
                    "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                    "and averaged such that there is a well-defined normal." );

            const Matrix< DDRMat > tTraction = tCMLeader->traction( mNormal, mStressType );
            tQI += mNitscheType * trans( tTraction ) * tM * tJump;
        }

        // add the contribution
        mSet->get_QI()( tQIIndex ) += aWStar * tQI;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Nitsche_Energy::compute_dQIdu( real aWStar )
    {
        
        // sint tQIIndex = mSet->get_QI_assembly_index( mName );
        // uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();
        // // compute dQIdu for indirect dof dependencies
        // for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
        // {            
            
        //     // get the treated dof type
        //     Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDof );

        //     // get leader index for residual dof type, indices for assembly
        //     uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
        //     uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        //     uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        //     const Matrix< DDRMat >
        //         tTraction( 8, 1, 1.0 );

        //     // compute dQIdu
        //     mSet->get_residual()( tQIIndex )(
        //                 { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
        //                 aWStar * 0.5 * ( tTraction );
            
        // }

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
