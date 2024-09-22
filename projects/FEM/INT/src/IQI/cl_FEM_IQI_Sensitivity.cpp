/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Sensitivity.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Sensitivity.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Sensitivity::IQI_Sensitivity() {}

    //------------------------------------------------------------------------------

    void
    IQI_Sensitivity::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        IQI::set_parameters( aParameters );

        // size of parameter list
        uint tParamSize = mParameters.size();

        // check that eigen vector index is defined
        MORIS_ERROR( tParamSize == 1,
                "IQI_Sensitivity::set_parameters - index of sensitivity vector needs to be specified.\n" );

        // extract eigenvalue index
        MORIS_ERROR( mParameters( 0 ).numel() == 1,
                "IQI_Sensitivity::set_parameters - size of sensitivity list needs to be one.\n" );

        mDerivativeVectorIndex = mParameters( 0 )( 0 );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Sensitivity::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // check if dof index was set (for the case of vector field)
        if ( mQuantityDofType.size() > 1 )
        {
            MORIS_ERROR( mIQITypeIndex != -1, "IQI_Sensitivity::compute_QI - mIQITypeIndex not set." );
        }
        else
        {
            mIQITypeIndex = 0;
        }

        Matrix< DDRMat > tMat;

        this->evaluate_QI( tMat );

        mSet->get_QI()( tQIIndex ) += aWStar * tMat;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Sensitivity::compute_QI( Matrix< DDRMat >& aQI )
    {
        // evaluate QI
        this->evaluate_QI( aQI );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Sensitivity::evaluate_QI( Matrix< DDRMat >& aMat )
    {
        // get field interpolator for a given dof type
        Field_Interpolator* tFI =
                mLeaderAdjointFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ), mDerivativeVectorIndex );

        // check that field interpolator exists
        MORIS_ASSERT( tFI != nullptr,
                "IQI_Sensitivity::evaluate_QI - field interpolator does not exist." );

        if ( mQuantityDofType.size() > 1 && mIQITypeIndex != -1 )
        {
            // evaluate DOF value
            aMat = { tFI->val()( mIQITypeIndex ) };
        }
        else
        {
            aMat = tFI->val();
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Sensitivity::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false, "IQI_Sensitivity::compute_dQIdu - function not implemented" );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Sensitivity::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "IQI_Sensitivity::compute_dQIdu - function not implemented" );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
