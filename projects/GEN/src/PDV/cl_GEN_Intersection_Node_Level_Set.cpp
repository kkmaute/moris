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
            mtk::Geometry_Type                    aBackgroundGeometryType,
            mtk::Interpolation_Order              aBackgroundInterpolationOrder,
            std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry )
            : Intersection_Node(
                    aNodeIndex,
                    aBaseNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aLocalCoordinate,
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    aInterfaceGeometry )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Intersection_Node_Level_Set::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // Get parent nodes
        const Basis_Node& tFirstParentNode = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();

        // Compute parent vector
        Matrix< DDRMat > tParentVector = trans( tSecondParentNode.get_global_coordinates() - tFirstParentNode.get_global_coordinates() );

        // Get sensitivity values from other ancestors
        Matrix< DDRMat > tSensitivitiesToAdd;
        const Cell< Basis_Node >& tFieldBasisNodes = this->get_field_basis_nodes();
        real tConstant = 1.0 / tFieldBasisNodes.size();
        for ( uint iLocatorNode = 0; iLocatorNode < tFieldBasisNodes.size(); iLocatorNode++ )
        {
            // Get geometry field sensitivity with respect to ADVs
            const Matrix< DDRMat >& tFieldSensitivities = tLockedInterfaceGeometry->get_dfield_dadvs(
                    tFieldBasisNodes( iLocatorNode ).get_index(),
                    tFieldBasisNodes( iLocatorNode ).get_global_coordinates() );

            // Ancestor sensitivities
            tSensitivitiesToAdd =
                    tConstant * aSensitivityFactor * this->get_dxi_dfield_from_ancestor( iLocatorNode ) * tParentVector * tFieldSensitivities;

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
        if ( tFirstParentNode.depends_on_advs() )
        {
            Matrix< DDRMat > tLocCoord = ( 1.0 - this->get_local_coordinate() ) * eye( tParentVector.n_rows(), tParentVector.n_rows() );
            Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_first_parent() );
            tFirstParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
        }

        // Add second parent coordinate sensitivities
        if ( tSecondParentNode.depends_on_advs() )
        {
            Matrix< DDRMat > tLocCoord = ( 1.0 + this->get_local_coordinate() ) * eye( tParentVector.n_rows(), tParentVector.n_rows() );
            Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_second_parent() );
            tSecondParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Intersection_Node_Level_Set::get_coordinate_determining_adv_ids() const
    {
        // Initialize ADV IDs
        Matrix< DDSMat > tCoordinateDeterminingADVIDs;

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
            Intersection_Node::join_adv_ids( tCoordinateDeterminingADVIDs, tAncestorADVIDs );
        }

        // Get parent nodes
        const Basis_Node& tFirstParentNode = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();

        // Add parent IDs
        if ( tFirstParentNode.depends_on_advs() )
        {
            Intersection_Node::join_adv_ids( tCoordinateDeterminingADVIDs, tFirstParentNode.get_coordinate_determining_adv_ids() );
        }
        if ( tSecondParentNode.depends_on_advs() )
        {
            Intersection_Node::join_adv_ids( tCoordinateDeterminingADVIDs, tSecondParentNode.get_coordinate_determining_adv_ids() );
        }

        // Return joined ADV IDs
        return tCoordinateDeterminingADVIDs;
    }

    //--------------------------------------------------------------------------------------------------------------

}