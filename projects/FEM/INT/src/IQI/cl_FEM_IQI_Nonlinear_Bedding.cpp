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
#include "cl_FEM_IQI_Nonlinear_Bedding.hpp"
#include "fn_norm.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Nonlinear_Bedding::IQI_Nonlinear_Bedding()
    {

        mFEMIQIType = fem::IQI_Type::NONLINEAR_BEDDING;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Bedding" ] = static_cast< uint >( IWG_Property_Type::BEDDING );
        mPropertyMap[ "Bedding_Threshold" ] = static_cast< uint >( IWG_Property_Type::BEDDING_THRESHOLD );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Nonlinear_Bedding::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get the parameter value for bedding
        const std::shared_ptr< Property >& tPropLeader =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

        // Get the parameter value for bedding threshold
        const std::shared_ptr< Property >& tPropBeddingThreshold =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING_THRESHOLD ) );

        // Obtain leader and follower field interpolators to get the displacement
        Field_Interpolator* tFILeader = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

        // Obtain displacement from the field interpolator
        const Matrix< DDRMat >& tDisplacement = tFILeader->val();

        // Transpose of displacement times displacement
        Matrix< DDRMat > tDispDotDisp = ( trans( tDisplacement ) * tDisplacement );

        // compute the bedding IQI value
        aQI = 0.5 * tPropLeader->val()( 0 ) * std::tanh( tDispDotDisp( 0 ) / ( tPropBeddingThreshold->val()( 0 ) ) ) * ( tDispDotDisp( 0 ) );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Nonlinear_Bedding::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the parameter value for bedding
        const std::shared_ptr< Property >& tPropLeader =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );\

        // Get the parameter value for bedding threshold
        const std::shared_ptr< Property >& tPropBeddingThreshold =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING_THRESHOLD ) );

        // Obtain leader and follower field interpolators to get the displacement
        Field_Interpolator* tFILeader = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

        // Obtain displacement from the field interpolator
        const Matrix< DDRMat >& tDisplacement = tFILeader->val();

        // Transpose of displacement times displacement
        Matrix< DDRMat > tDispDotDisp = ( trans( tDisplacement ) * tDisplacement );

        // compute the bedding IQI value
        Matrix< DDRMat > tQI = 0.5 * tPropLeader->val() * std::tanh( tDispDotDisp( 0 ) / ( tPropBeddingThreshold->val()( 0 ) ) ) * ( tDispDotDisp( 0 ) );

        // add the contribution
        mSet->get_QI()( tQIIndex ) += aWStar * tQI;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Nonlinear_Bedding::compute_dQIdu( real aWStar )
    {
        // check the point is inside the bounded box
        if ( !this->is_within_box_bounds() )
        {
            return;
        }

        // get the column index to assemble in residual
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the parameter value for bedding
        const std::shared_ptr< Property >& tPropLeader =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

        // Get the parameter value for bedding threshold
        const std::shared_ptr< Property >& tPropBeddingThreshold =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING_THRESHOLD ) );

        // get the number of leader dof type dependencies
        uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // compute dQIdu for indirect dof dependencies
        for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
        {
            // get the treated dof type
            Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDof );

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            if ( tDofType( 0 ) == MSI::Dof_Type::UX )
            {
                // get field interpolator for displacements
                Field_Interpolator* tDisplacementFI =
                        mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                // Get displacement
                const Matrix< DDRMat >& tDisplacement = tDisplacementFI->val();

                // Get basis functions
                const Matrix< DDRMat >& tN = tDisplacementFI->N();

                // Displacment transpose times displacement
                Matrix< DDRMat > tDispDotDisp = trans( tDisplacement ) * tDisplacement;

                // compute bedding contribution - displacements

                mSet->get_residual()( tQIIndex )(
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * 0.5 * tPropLeader->val() * ( ( tDispDotDisp( 0 ) ) * ( trans( tDisplacement ) * tN ) * ( 2.0 / ( 1.0 + std::pow( tDispDotDisp( 0 ) / ( tPropBeddingThreshold->val()( 0 ) ), 2 ) ) )
                        + 2.0 * std::tanh( tDispDotDisp( 0 ) / ( tPropBeddingThreshold->val()( 0 ) ) ) * ( trans( tDisplacement ) * tN ) );

            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Nonlinear_Bedding::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {

        // check the point is inside the bounded box
        if ( !this->is_within_box_bounds() )
        {
            return;
        }

        // get the bedding stabilization parameter
        const std::shared_ptr< Property >& tPropLeader =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );
        
                // Get the parameter value for bedding threshold
        const std::shared_ptr< Property >& tPropBeddingThreshold =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING_THRESHOLD ) );

        // initialize derivative
        adQIdu.fill( 0.0 );

        // if DoF type is displacements
        if ( aDofType( 0 ) == MSI::Dof_Type::UX )
        {
            // get field interpolator for displacements
            Field_Interpolator* tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofType( 0 ) );

            // Get displacement
            const Matrix< DDRMat >& tDisplacement = tDisplacementFI->val();

            // Get basis functions
            const Matrix< DDRMat >& tN = tDisplacementFI->N();

            // Displacment transpose times displacement
            Matrix< DDRMat > tDispDotDisp = trans( tDisplacement ) * tDisplacement;

            // compute dQIdu
            adQIdu = 0.5 * tPropLeader->val()( 0 ) * ( ( tDispDotDisp( 0 ) ) * ( trans( tDisplacement ) * tN ) * ( 2.0 / ( 1.0 + std::pow( tDispDotDisp( 0 ) / ( tPropBeddingThreshold->val()( 0 ) ), 2 ) ) )
                     + 2.0 * std::tanh( tDispDotDisp( 0 ) / ( tPropBeddingThreshold->val()( 0 ) ) ) * ( trans( tDisplacement ) * tN ) );

        }

    }

    //-------------------------------------------------------------------------------------------------------------------------------------

    bool IQI_Nonlinear_Bedding::is_within_box_bounds()
    {
        // check if the box bounds are empty then skip
        if ( mParameters.empty() )
        {
            return true;
        }

        // if the box bounds are not empty then check if it is inside the box
        else
        {
            // get the coordinate
            const Matrix< DDRMat >& tGaussPoint = mLeaderFIManager->get_IG_geometry_interpolator()->valx();

            // check if the calculation point coordinates are more then lower corner of the box
            bool tLowerBound = std::equal( mParameters( 0 ).begin(), mParameters( 0 ).end(), tGaussPoint.begin(), tGaussPoint.end(), []( real aA, real aB ) -> bool { return aA < aB; } );

            // check if the calculation point coordinates are less then upper corner of the box
            bool tUpperBound = std::equal( tGaussPoint.begin(), tGaussPoint.end(), mParameters( 1 ).begin(), mParameters( 1 ).end(), []( real aA, real aB ) -> bool { return aA < aB; } );

            // combine the two bounds that satisfy both
            return tUpperBound and tLowerBound;
        }
    }


    //------------------------------------------------------------------------------

}    // namespace moris::fem
