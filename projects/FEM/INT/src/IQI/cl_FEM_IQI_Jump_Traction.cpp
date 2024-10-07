/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Jump_Traction.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Jump_Traction.hpp"
#include "fn_norm.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Jump_Traction::IQI_Jump_Traction()
    {
        mFEMIQIType = fem::IQI_Type::JUMP_TRACTION;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "TractionCM" ] = (uint)IQI_Constitutive_Type::TRACTION_CM;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Jump_Traction::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get the constitutive model for computing the fluxes/tractions
        const std::shared_ptr< Constitutive_Model >& tCMLeader =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );
        const std::shared_ptr< Constitutive_Model >& tCMFollower =
                mFollowerCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );

        MORIS_ASSERT( mNormal.numel() > 0,
                "IQI_Jump_Traction::compute_QI() - "
                "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                "and averaged such that there is a well-defined normal." );

        // evaluate traction difference
        const Matrix< DDRMat >
                tTractionJump =
                        tCMLeader->traction( mNormal ) - tCMFollower->traction( -1.0 * mNormal );

        // based on the IQI index select the norm or the individual component
        if ( mIQITypeIndex == -1 )
        {
            // compute and return magnitude
            real tNorm = norm( tTractionJump );
            aQI        = { { tNorm } };
        }
        else
        {
            // pick the component otherwise (0,1,2)
            aQI = { tTractionJump( mIQITypeIndex ) };
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Jump_Traction::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the constitutive model for computing the fluxes/tractions
        const std::shared_ptr< Constitutive_Model >& tCMLeader =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );
        const std::shared_ptr< Constitutive_Model >& tCMFollower =
                mFollowerCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );

        MORIS_ASSERT( mNormal.numel() > 0,
                "IQI_Jump_Traction::compute_QI() - "
                "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                "and averaged such that there is a well-defined normal." );

        // evaluate traction difference
        const Matrix< DDRMat >
                tTractionJump =
                        tCMLeader->traction( mNormal ) - tCMFollower->traction( -1.0 * mNormal );

        // initialize the output matrix
        Matrix< DDRMat > tMat;

        // based on the IQI index select the norm or the individual component
        if ( mIQITypeIndex == -1 )
        {
            // compute and return magnitude
            real tNorm = norm( tTractionJump );
            tMat       = { { tNorm } };
        }
        else
        {
            // pick the component otherwise (0,1,2)
            tMat = { { tTractionJump( mIQITypeIndex ) } };
        }

        // add the contribution
        mSet->get_QI()( tQIIndex ) += aWStar * tMat;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Jump_Traction::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Jump_Traction::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------

}    // namespace moris::fem
