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
#include "cl_FEM_IQI_Traction.hpp"
#include "fn_norm.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Traction::IQI_Traction(           
            enum CM_Function_Type aStressType,
            enum CM_Function_Type aStrainType
        )
    {
        // Assign stress and strain type
        mStressType = aStressType;
        mStrainType = aStrainType;
    
        mFEMIQIType = fem::IQI_Type::TRACTION;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "TractionCM" ] = (uint)IQI_Constitutive_Type::TRACTION_CM;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Traction::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get the constitutive model for computing the fluxes/tractions
        const std::shared_ptr< Constitutive_Model >& tCMLeader =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );

        MORIS_ASSERT( mNormal.numel() > 0,
                "IQI_Traction::compute_QI() - "
                "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                "and averaged such that there is a well-defined normal." );

        // evaluate traction
        const Matrix< DDRMat >
                tTraction = tCMLeader->traction( mNormal , mStressType );

        // based on the IQI index select the norm or the individual component
        if ( mIQITypeIndex == -1 )
        {
            // compute and return magnitude
            real tNorm = norm( tTraction );
            aQI        = { tNorm };
        }
        else
        {
            // pick the component otherwise (0,1,2)
            aQI = { tTraction( mIQITypeIndex ) };
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Traction::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the constitutive model for computing the fluxes/tractions
        const std::shared_ptr< Constitutive_Model >& tCMLeader =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );

        MORIS_ASSERT( mNormal.numel() > 0,
                "IQI_Traction::compute_QI() - "
                "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                "and averaged such that there is a well-defined normal." );

        // evaluate traction
        const Matrix< DDRMat >
                tTraction = tCMLeader->traction( mNormal , mStressType );

        // initialize the output matrix
        Matrix< DDRMat > tMat;

        // based on the IQI index select the norm or the individual component
        if ( mIQITypeIndex == -1 )
        {
            // compute and return magnitude
            real tNorm = norm( tTraction );
            tMat       = { { tNorm } };
        }
        else
        {
            // pick the component otherwise (0,1,2)
            tMat = { tTraction( mIQITypeIndex ) };
        }

        // add the contribution
        mSet->get_QI()( tQIIndex ) += aWStar * tMat;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Traction::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Traction::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
