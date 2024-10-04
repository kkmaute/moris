/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Dof_Sensitivity.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Dof_Sensitivity.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Dof_Sensitivity::IQI_Dof_Sensitivity() {}

    //------------------------------------------------------------------------------

    void
    IQI_Dof_Sensitivity::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        IQI::set_parameters( aParameters );

        // size of parameter list
        uint tParamSize = mParameters.size();

        // check that eigen vector index is defined
        MORIS_ERROR( tParamSize == 1,
                "IQI_Dof_Sensitivity::set_parameters - index of sensitivity vector needs to be specified.\n" );

        // extract optimization variable index index
        MORIS_ERROR( mParameters( 0 ).numel() == 1,
                "IQI_Dof_Sensitivity::set_parameters - size of sensitivity list needs to be one.\n" );

        mOptimizationVariableId = mParameters( 0 )( 0 );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Dof_Sensitivity::compute_QI( real aWStar )
    {
        MORIS_ERROR( false, "IQI_Dof_Sensitivity::compute_QI( real aWStar ) - function not implemented" );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Dof_Sensitivity::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get field interpolator for a given dof type
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        Field_Interpolator* tFIsens =
                mLeaderAdjointFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ), mOptimizationVariableId );

        // check that field interpolator exists
        MORIS_ASSERT( tFI != nullptr,
                "IQI_Dof_Sensitivity::evaluate_QI - field interpolator does not exist." );

        if ( mQuantityDofType.size() > 1 && mIQITypeIndex != -1 )
        {
            // evaluate DOF value
            aQI = { tFIsens->val()( mIQITypeIndex ) };
        }
        else
        {
            aQI = tFIsens->val();
        }

        // check whether optimization variable influences this node
        sint tLocalAdvIndex = mSet->get_current_adv_geo_index( mOptimizationVariableId );

        if ( tLocalAdvIndex >= 0 )
        {
            // get partial derivative of node location wrt adv
            const Matrix< DDRMat >& tDxDp = mSet->get_current_adv_geo_weight();

            if ( tDxDp.numel() > 0 )
            {
                const Matrix< DDRMat >& tSpatialGradient = tFI->gradx( 1 );

                if ( mQuantityDofType.size() > 1 && mIQITypeIndex != -1 )
                {

                    aQI( 0 ) += dot( tSpatialGradient.get_column( mIQITypeIndex ), tDxDp.get_column( tLocalAdvIndex ) );
                }
                else
                {
                    aQI( 0 ) += dot( tSpatialGradient, tDxDp );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Dof_Sensitivity::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false, "IQI_Dof_Sensitivity::compute_dQIdu - function not implemented" );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Dof_Sensitivity::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "IQI_Dof_Sensitivity::compute_dQIdu - function not implemented" );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
