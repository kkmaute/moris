/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_PDV_Sensitivity.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_PDV_Sensitivity.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_PDV_Sensitivity::IQI_PDV_Sensitivity()
    {
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        mPropertyMap[ "L2Check" ] = static_cast< uint >( IQI_Property_Type::L2_CHECK );
        mPropertyMap[ "Gradx" ]   = static_cast< uint >( IQI_Property_Type::SPATIAL_GRADIENT );

        mFEMIQIType = fem::IQI_Type::PDV_SENSITIVITY;
    }

    //------------------------------------------------------------------------------

    void
    IQI_PDV_Sensitivity::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        IQI::set_parameters( aParameters );

        // size of parameter list
        uint tParamSize = mParameters.size();

        // check that eigen vector index is defined
        MORIS_ERROR( tParamSize == 1,
                "IQI_PDV_Sensitivity::set_parameters - index of sensitivity vector needs to be specified.\n" );

        // extract optimization variable index index
        MORIS_ERROR( mParameters( 0 ).numel() == 2,
                "IQI_PDV_Sensitivity::set_parameters - size of sensitivity list needs to be 2.\n" );

        mPDVTypeIndex           = mParameters( 0 )( 0 );
        mOptimizationVariableId = mParameters( 0 )( 1 );
    }

    //------------------------------------------------------------------------------

    void
    IQI_PDV_Sensitivity::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get field interpolator for a given field type
        Field_Interpolator* tFI = mLeaderFIManager->get_field_interpolators_for_type( mLeaderDvTypes( mPDVTypeIndex ) );

        // Field_Interpolator* tFIsens = mLeaderAdjointFIManager->get_field_interpolators_for_type( mLeaderDvTypes( mPDVTypeIndex ), mOptimizationVariableId );

        // BRENDAN assumes the PDV_Type field does not depend on ADVs
        // // check that field interpolator exists
        // MORIS_ASSERT( tFI != nullptr,
        //         "IQI_PDV_Sensitivity::evaluate_QI - field interpolator does not exist." );

        // if ( mLeaderDvTypes.size() > 1 && mIQITypeIndex != -1 )
        // {
        //     // evaluate field value
        //     aQI = { tFIsens->val()( mIQITypeIndex ) };
        // }
        // else
        // {
        //     aQI = tFIsens->val();
        // }

        // check whether optimization variable influences this node
        sint tLocalAdvIndex = mSet->get_current_adv_geo_index( mOptimizationVariableId );

        real tValue = 0.0;

        if ( tLocalAdvIndex >= 0 )
        {
            // get partial derivative of node location wrt adv
            const Matrix< DDRMat >& tDxDp = mLeaderAdjointFIManager->get_field_interpolators_for_type( gen::PDV_Type::DXDP )->val();

            if ( tDxDp.numel() > 0 )
            {
                const Matrix< DDRMat >& tSpatialGradient = tFI->gradx( 1 );
                // const Matrix< DDRMat >& tSpatialGradientProp = tPropGradx->val();

                tValue = dot( tSpatialGradient, tDxDp );
            }
        }

        // See if we want to do a check on the L2 value
        std::shared_ptr< Property >& tPropL2Check = mLeaderProp( static_cast< uint >( IQI_Property_Type::L2_CHECK ) );
        if ( tPropL2Check != nullptr )
        {
            Matrix< DDRMat > tL2CheckVal = tPropL2Check->val();
            real             tJumpNorm   = std::abs( tValue - norm( tL2CheckVal ) );

            mSet->get_QI()( tQIIndex ) += aWStar * tJumpNorm * tJumpNorm;
        }
        else
        {
            mSet->get_QI()( tQIIndex ) += aWStar * tValue;
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_PDV_Sensitivity::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get field interpolator for a given field type
        Field_Interpolator* tFI = mLeaderFIManager->get_field_interpolators_for_type( mLeaderDvTypes( mPDVTypeIndex ) );

        std::shared_ptr< Property >& tPropGradx = mLeaderProp( static_cast< uint >( IQI_Property_Type::SPATIAL_GRADIENT ) );

        // Field_Interpolator* tFIsens = mLeaderAdjointFIManager->get_field_interpolators_for_type( mLeaderDvTypes( mPDVTypeIndex ), mOptimizationVariableId );

        // BRENDAN assumes the PDV_Type field does not depend on ADVs
        // // check that field interpolator exists
        // MORIS_ASSERT( tFI != nullptr,
        //         "IQI_PDV_Sensitivity::evaluate_QI - field interpolator does not exist." );

        // if ( mLeaderDvTypes.size() > 1 && mIQITypeIndex != -1 )
        // {
        //     // evaluate field value
        //     aQI = { tFIsens->val()( mIQITypeIndex ) };
        // }
        // else
        // {
        //     aQI = tFIsens->val();
        // }

        // check whether optimization variable influences this node
        sint tLocalAdvIndex = mSet->get_current_adv_geo_index( mOptimizationVariableId );

        if ( tLocalAdvIndex >= 0 )
        {
            // get partial derivative of node location wrt adv
            const Matrix< DDRMat >& tDxDp = mLeaderAdjointFIManager->get_field_interpolators_for_type( gen::PDV_Type::DXDP )->val();

            if ( tDxDp.numel() > 0 )
            {
                if ( tPropGradx == nullptr )
                {
                    const Matrix< DDRMat >& tSpatialGradient = tFI->gradx( 1 );
                    aQI( 0 ) += dot( tSpatialGradient, tDxDp );
                }
                else
                {
                    const Matrix< DDRMat >& tSpatialGradient = tPropGradx->val();
                    aQI( 0 ) += dot( tSpatialGradient, tDxDp );
                }


                if ( norm( aQI ) > 1e-12 and this->get_name() == "IQIBulkL2ErrorSens" )
                {
                    PRINT( tFI->gradx( 1 ) );
                    PRINT( tDxDp );
                    PRINT( aQI );
                    std::cout << std::endl;
                }
            }
        }

        // See if we want to do a check on the L2 value
        std::shared_ptr< Property >& tPropL2Check = mLeaderProp( static_cast< uint >( IQI_Property_Type::L2_CHECK ) );
        if ( tPropL2Check != nullptr )
        {
            Matrix< DDRMat > tL2CheckVal = tPropL2Check->val();
            real             tJumpNorm   = norm( aQI - tL2CheckVal );

            aQI = { tJumpNorm * tJumpNorm };
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_PDV_Sensitivity::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false, "IQI_PDV_Sensitivity::compute_dQIdu - function not implemented" );
        return;
    }

    //------------------------------------------------------------------------------

    void
    IQI_PDV_Sensitivity::compute_dQIdu(
            Vector< MSI::Dof_Type >& aFieldType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "IQI_PDV_Sensitivity::compute_dQIdu - function not implemented" );
        return;
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
