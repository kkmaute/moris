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
#include "cl_GEN_Basis_Node.hpp"

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

    void
    Intersection_Node_Level_Set::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor )
    {
        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // Get sensitivity values from other ancestors
        Matrix< DDRMat > tSensitivitiesToAdd;
        const Cell< Basis_Node >& tLocatorNodes = this->get_field_basis_nodes();
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

        // Add first parent coordinate sensitivities
        Matrix< DDRMat > tLocCoord = ( 1.0 - this->get_local_coordinate() ) * eye( mParentVector.n_rows(), mParentVector.n_rows() );
        Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + mParentVector * this->get_dxi_dcoordinate_first_parent() );
        this->get_first_parent_node().append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );

        // Add second parent coordinate sensitivities
        tLocCoord = ( 1.0 + this->get_local_coordinate() ) * eye( mParentVector.n_rows(), mParentVector.n_rows() );
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
