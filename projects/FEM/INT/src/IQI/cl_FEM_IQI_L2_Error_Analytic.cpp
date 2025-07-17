/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_L2_Error_Analytic.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_L2_Error_Analytic.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IQI_L2_Error_Analytic::IQI_L2_Error_Analytic()
    {
        // set FEM IQI type
        mFEMIQIType = fem::IQI_Type::L2_ERROR_ANALYTIC;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "L2Check" ] = static_cast< uint >( IQI_Property_Type::L2_CHECK );
    }

    //------------------------------------------------------------------------------

    void
    IQI_L2_Error_Analytic::compute_QI( Matrix< DDRMat >& aQI )
    {
        // check that mQuantityDofType is set
        MORIS_ERROR( mQuantityDofType.size() > 0,
                "IQI_L2_Error_Analytic::compute_QI - mQuantityDofType not set." );

        // get field interpolator
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        real tJumpNorm;

        Field_Interpolator* tFIField = nullptr;
        if ( mLeaderFieldTypes.size() != 0 )
        {
            tFIField = mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            tJumpNorm = norm( tFI->val() - tFIField->val() );
        }
        else
        {
            // get analytical solution property
            std::shared_ptr< Property >& tPropL2Check =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::L2_CHECK ) );
            tJumpNorm = norm( tFI->val() - tPropL2Check->val() );
        }

        // evaluate the QI
        aQI = tJumpNorm * tJumpNorm;
    }

    //------------------------------------------------------------------------------

    void
    IQI_L2_Error_Analytic::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // check that mQuantityDofType is set
        MORIS_ERROR( mQuantityDofType.size() > 0,
                "IQI_L2_Error_Analytic::compute_QI - mQuantityDofType not set." );

        // get field interpolator
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // get analytical solution property
        std::shared_ptr< Property >& tPropL2Check =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::L2_CHECK ) );

        // compute jump
        real tJumpNorm = norm( tFI->val() - tPropL2Check->val() );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * ( tJumpNorm * tJumpNorm );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
