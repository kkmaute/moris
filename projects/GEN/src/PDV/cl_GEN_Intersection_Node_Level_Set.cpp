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

namespace moris::gen
{
    Intersection_Node_Level_Set::Intersection_Node_Level_Set(
            uint                     aNodeIndex,
            const Cell< Background_Node* >& aBackgroundNodes,
            const Parent_Node&       aFirstParentNode,
            const Parent_Node&       aSecondParentNode,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder,
            Level_Set_Geometry&      aInterfaceGeometry )
            : Intersection_Node(
                    aNodeIndex,
                    aBackgroundNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aInterfaceGeometry.compute_intersection_local_coordinate( aBackgroundNodes, aFirstParentNode, aSecondParentNode ),
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Intersection_Node_Level_Set::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        // Get parent nodes
        const Basis_Node& tFirstParentNode = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();

        // Compute parent vector
        Matrix< DDRMat > tParentVector = trans( tSecondParentNode.get_global_coordinates() - tFirstParentNode.get_global_coordinates() );

        // Get sensitivity values from other ancestors
        Matrix< DDRMat > tSensitivitiesToAdd;
        const Cell< Basis_Node >& tFieldBasisNodes = this->get_field_basis_nodes();
        for ( uint iFieldBasisNode = 0; iFieldBasisNode < tFieldBasisNodes.size(); iFieldBasisNode++ )
        {
            // Get geometry field sensitivity with respect to ADVs
            const Matrix< DDRMat >& tFieldSensitivities = mInterfaceGeometry.get_dfield_dadvs(
                    tFieldBasisNodes( iFieldBasisNode ).get_index(),
                    tFieldBasisNodes( iFieldBasisNode ).get_global_coordinates() );

            // Ancestor sensitivities
            tSensitivitiesToAdd =
                    0.5 * aSensitivityFactor * this->get_dxi_dfield_from_ancestor( iFieldBasisNode ) * tParentVector * tFieldSensitivities;

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

        // Get sensitivity values from other ancestors
        const Cell< Basis_Node >& tFieldBasisNodes = this->get_field_basis_nodes();
        for ( uint iFieldBasisNode = 0; iFieldBasisNode < tFieldBasisNodes.size(); iFieldBasisNode++ )
        {
            // Get geometry field sensitivity with respect to ADVs
            const Matrix< DDSMat >& tAncestorADVIDs = mInterfaceGeometry.get_determining_adv_ids(
                    tFieldBasisNodes( iFieldBasisNode ).get_index(),
                    tFieldBasisNodes( iFieldBasisNode ).get_global_coordinates() );

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

    Geometry& Intersection_Node_Level_Set::get_interface_geometry()
    {
        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Geometry& Intersection_Node_Level_Set::get_interface_geometry() const
    {
        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

}
