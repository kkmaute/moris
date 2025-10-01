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
#include "cl_FEM_IQI_Traction_Potential.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Traction_Potential::IQI_Traction_Potential()
    {

        mFEMIQIType = fem::IQI_Type::TRACTION_POTENTIAL;

        // set size for the constitutive model pointer cell
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );


        // populate the constitutive map
        mPropertyMap[ "Traction" ] = (uint)IQI_Property_Type::TRACTION;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Traction_Potential::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get the property model for the tractions
        const std::shared_ptr< Property >& tPropLeader =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::TRACTION ) );

        MORIS_ASSERT( mNormal.numel() > 0,
                "IQI_Traction::compute_QI() - "
                "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                "and averaged such that there is a well-defined normal." );

        // Obtain leader and follower field interpolators to get the displacement
        Field_Interpolator* tFILeader   = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

        // Compute traction potential energy
        aQI = trans( tFILeader->val() ) * tPropLeader->val();

    }

    //------------------------------------------------------------------------------

    void
    IQI_Traction_Potential::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the property model for the tractions
        const std::shared_ptr< Property >& tPropLeader =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::TRACTION ) );

        MORIS_ASSERT( mNormal.numel() > 0,
                "IQI_Traction::compute_QI() - "
                "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                "and averaged such that there is a well-defined normal." );

        // Obtain leader and follower field interpolators to get the displacement
        Field_Interpolator* tFILeader   = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

        // Compute traction potential energy
        Matrix< DDRMat > tMat = trans( tFILeader->val() ) * tPropLeader->val();

        // add the contribution
        mSet->get_QI()( tQIIndex ) += aWStar * tMat;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Traction_Potential::compute_dQIdu( real aWStar )
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
        //         tTraction( 8, 1, 0.0 );

        //     // compute dQIdu
        //     mSet->get_residual()( tQIIndex )(
        //                 { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
        //                 aWStar * 0.5 * ( tTraction );
            
        // }
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Traction_Potential::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
