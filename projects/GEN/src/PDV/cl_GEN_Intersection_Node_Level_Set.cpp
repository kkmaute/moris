/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Level_Set.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Level_Set.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"

#include "fn_eye.hpp"
#include "fn_norm.hpp"

namespace moris::ge
{
    Intersection_Node_Level_Set::Intersection_Node_Level_Set(
            uint                                  aNodeIndex,
            const Cell< Node* >&                  aBaseNodes,
            const Parent_Node&                    aFirstParentNode,
            const Parent_Node&                    aSecondParentNode,
            real                                  aLocalCoordinate,
            mtk::Geometry_Type                    aBaseGeometryType,
            std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry )
            : Intersection_Node(
                    aNodeIndex,
                    aBaseNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aLocalCoordinate,
                    aBaseGeometryType,
                    aInterfaceGeometry )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Intersection_Node_Level_Set::determine_is_intersected()
    {
        // Lock the interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // Get the difference between the level set value of the parents and the isocontour threshold
        real tFirstDiffFromThreshold = tLockedInterfaceGeometry->get_field_value( this->get_first_parent_node().get_index(), this->get_first_parent_node().get_global_coordinates() );
        real tSecondDiffFromThreshold = tLockedInterfaceGeometry->get_field_value( this->get_second_parent_node().get_index(), this->get_second_parent_node().get_global_coordinates() );

        // get the isocontour thresholds from the geometry
        real tLocalCoordinate = this->get_local_coordinate();

        // Determine if edge is intersected
        bool tIsIntersected;
        if ( this->is_first_parent_on_interface() or this->is_second_parent_on_interface() )
        {
            return true;
        }
        // FIXME: This check should be unnecessary as the local edge coordinate should be sufficient
        // to determine whether edge is intersected; it is only "useful" if parent node's level set value
        // is determined by method that is different from intersection nodes; for example level set value child node
        // of child node is computed via analytic field and intersection node via bi-linear interpolation
        else if ( tFirstDiffFromThreshold * tSecondDiffFromThreshold > 0 )
        {
            tIsIntersected = false;

            // check for consistency of parent values and local coordinate
            MORIS_ASSERT( std::abs( tLocalCoordinate ) > 1,
                    "Intersection_Node::Intersection_Node - inconsistent parent level set values versus local coordinate - p1 %e p2 %e loc %e.",
                    tFirstDiffFromThreshold,
                    tSecondDiffFromThreshold,
                    tLocalCoordinate );
        }
        else
        {
            tIsIntersected = ( std::abs( tLocalCoordinate ) <= 1.0 );

            // check for consistency with parent values
            // this check is currently useless but should be performed is inconsistency issue (see comment above) is resolved
            MORIS_ASSERT( tIsIntersected ? tFirstDiffFromThreshold * tSecondDiffFromThreshold < 0 : tFirstDiffFromThreshold * tSecondDiffFromThreshold > 0,
                    "Intersection_Node::Intersection_Node - inconsistent parent level set values - p1 %e p2 %e loc %e.",
                    tFirstDiffFromThreshold,
                    tSecondDiffFromThreshold,
                    tLocalCoordinate );
        }

        return tIsIntersected;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Intersection_Node_Level_Set::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor )
    {
        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // Get sensitivity values from other ancestors
        Matrix< DDRMat > tSensitivitiesToAdd;
        const Cell< Basis_Node >& tLocatorNodes = this->get_locator_nodes();
        for ( uint iLocatorNode = 0; iLocatorNode < tLocatorNodes.size(); iLocatorNode++ )
        {
            // Get geometry field sensitivity with respect to ADVs
            const Matrix< DDRMat >& tFieldSensitivities = tLockedInterfaceGeometry->get_dfield_dadvs(
                    tLocatorNodes( iLocatorNode ).get_index(),
                    tLocatorNodes( iLocatorNode ).get_global_coordinates() );

            // Ancestor sensitivities
            tSensitivitiesToAdd =
                    0.5 * aSensitivityFactor * this->get_dxi_dfield_from_ancestor( iLocatorNode ) * mParentVector * tFieldSensitivities;

            // Resize sensitivities
            uint tJoinedSensitivityLength = aCoordinateSensitivities.n_cols();
            aCoordinateSensitivities.resize( tSensitivitiesToAdd.n_rows(),
                    tJoinedSensitivityLength + tSensitivitiesToAdd.n_cols() );

            // Join sensitivities
            for ( uint tCoordinateIndex = 0; tCoordinateIndex < tSensitivitiesToAdd.n_rows(); tCoordinateIndex++ )
            {
                for ( uint tAddedSensitivity = 0; tAddedSensitivity < tSensitivitiesToAdd.n_cols(); tAddedSensitivity++ )
                {
                    aCoordinateSensitivities( tCoordinateIndex, tJoinedSensitivityLength + tAddedSensitivity ) =
                            tSensitivitiesToAdd( tCoordinateIndex, tAddedSensitivity );
                }
            }
        }

        // Add first parent sensitivities
        Matrix< DDRMat > tLocCoord = ( 1.0 - this->get_local_coordinate() ) *
                                     eye( mParentVector.n_rows(), mParentVector.n_rows() );

        Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + mParentVector * this->get_dxi_dcoordinate_first_parent() );
        this->get_first_parent_node().append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );

        // Add second parent sensitivities
        tLocCoord = ( 1.0 + this->get_local_coordinate() ) *
                                     eye( mParentVector.n_rows(), mParentVector.n_rows() );

        tSensitivityFactor = 0.5 * ( tLocCoord + mParentVector * this->get_dxi_dcoordinate_second_parent() );
        this->get_second_parent_node().append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Intersection_Node_Level_Set::get_coordinate_determining_adv_ids()
    {
        // Set size
        mCoordinateDeterminingADVIDs.set_size( 0, 0 );

        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // Get sensitivity values from other ancestors
        const Cell< Basis_Node >& tLocatorNodes = this->get_locator_nodes();
        for ( uint iLocatorNode = 0; iLocatorNode < tLocatorNodes.size(); iLocatorNode++ )
        {
            // Get geometry field sensitivity with respect to ADVs
            const Matrix< DDSMat >& tAncestorADVIDs = tLockedInterfaceGeometry->get_determining_adv_ids(
                    tLocatorNodes( iLocatorNode ).get_index(),
                    tLocatorNodes( iLocatorNode ).get_global_coordinates() );

            // Join IDs
            this->join_adv_ids( tAncestorADVIDs );
        }

        // Add parent IDs
        this->join_adv_ids( this->get_first_parent_node().get_coordinate_determining_adv_ids() );
        this->join_adv_ids( this->get_second_parent_node().get_coordinate_determining_adv_ids() );

        // Return joined ADV IDs
        return mCoordinateDeterminingADVIDs;
    }

    //--------------------------------------------------------------------------------------------------------------

}
