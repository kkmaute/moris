/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_H1_Error_Analytic.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_H1_Error_Analytic.hpp"

#include "fn_vectorize.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IQI_H1_Error_Analytic::IQI_H1_Error_Analytic()
    {
        // set FEM IQI type
        mFEMIQIType = fem::IQI_Type::H1_ERROR_ANALYTIC;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "H1Check" ] = static_cast< uint >( IQI_Property_Type::H1_CHECK );
    }

    //------------------------------------------------------------------------------

    void
    IQI_H1_Error_Analytic::compute_QI( Matrix< DDRMat >& aQI )
    {
        // check that mQuantityDofType is set
        MORIS_ERROR( mQuantityDofType.size() > 0,
                "IQI_L2_Error_Analytic::compute_QI - mQuantityDofType not set." );

        // get field interpolator
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // get analytical solution property
        std::shared_ptr< Property >& tPropH1Check =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::H1_CHECK ) );

        // get jump between value and analytic
        real tJumpNorm = norm( vectorize( tFI->gradx( 1 ) - tPropH1Check->val() ) );

        // evaluate the QI
        aQI = tJumpNorm * tJumpNorm;
    }

    //------------------------------------------------------------------------------

    void
    IQI_H1_Error_Analytic::compute_QI( real aWStar )
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
        std::shared_ptr< Property >& tPropH1Check =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::H1_CHECK ) );

        // get jump between value and analytic
        real tJumpNorm = norm( vectorize( tFI->gradx( 1 ) - tPropH1Check->val() ) );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * ( tJumpNorm * tJumpNorm );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
