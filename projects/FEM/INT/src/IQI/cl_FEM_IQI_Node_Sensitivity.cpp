/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Node_Sensitivity.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Node_Sensitivity.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Node_Sensitivity::IQI_Node_Sensitivity() {}

    //------------------------------------------------------------------------------

    void
    IQI_Node_Sensitivity::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
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
    IQI_Node_Sensitivity::compute_QI( real aWStar )
    {
        MORIS_ERROR( false, "IQI_Node_Sensitivity::compute_QI( real aWStar ) - function not implemented" );
    }

    //------------------------------------------------------------------------------

    void IQI_Node_Sensitivity::compute_QI( Matrix< DDRMat >& aQI )
    {
        // check if direction index was set (needed as dxdp is a vector field)
        MORIS_ERROR( mIQITypeIndex != -1, "IQI_Node_Sensitivity::compute_QI - mIQITypeIndex not set." );

        // check whether optimization variable influences this node
        sint tLocalAdvIndex = mSet->get_current_adv_geo_index( mOptimizationVariableId );

        if ( tLocalAdvIndex >= 0 )
        {
            // get dxdp field
            const Matrix< DDRMat >& tDxDp = mSet->get_current_adv_geo_weight();

            if ( tDxDp.numel() > 0 )
            {
                // evaluate DOF value
                aQI = { tDxDp( mIQITypeIndex, tLocalAdvIndex ) };

                return;
            }
        }

        // if optimization variable has no influence, set criterion to zero
        aQI = { 0 };
    }

    //------------------------------------------------------------------------------

    void IQI_Node_Sensitivity::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false, "IQI_Node_Sensitivity::compute_dQIdu - function not implemented" );
    }

    //------------------------------------------------------------------------------

    void IQI_Node_Sensitivity::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "IQI_Node_Sensitivity::compute_dQIdu - function not implemented" );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
