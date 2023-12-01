/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Surface_Mesh.cpp
 *
 */

#include <limits>
#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"

namespace moris::ge
{
    Intersection_Node_Surface_Mesh::Intersection_Node_Surface_Mesh(
            std::shared_ptr< Intersection_Node >     aFirstParentNode,
            std::shared_ptr< Intersection_Node >     aSecondParentNode,
            uint                                     aFirstParentNodeIndex,
            uint                                     aSecondParentNodeIndex,
            const Matrix< DDRMat >&                  aFirstParentNodeCoordinates,
            const Matrix< DDRMat >&                  aSecondParentNodeCoordinates,
            const Matrix< DDUMat >&                  aAncestorNodeIndices,
            const Cell< Matrix< DDRMat > >&          aAncestorNodeCoordinates,
            std::shared_ptr< Surface_Mesh_Geometry > aInterfaceGeometry )
            : Intersection_Node(
                    compute_local_coordinate(
                            aFirstParentNodeCoordinates,
                            aSecondParentNodeCoordinates,
                            aAncestorNodeCoordinates ),
                    aFirstParentNode,
                    aSecondParentNode,
                    aFirstParentNodeIndex,
                    aSecondParentNodeIndex,
                    { { -1 } },
                    { { 1 } },
                    aAncestorNodeIndices,
                    aAncestorNodeCoordinates,
                    Element_Intersection_Type::Linear_1D )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
        // Ensure that only two ancestor nodes are provided
        MORIS_ERROR( aAncestorNodeCoordinates.size() != 2, "GEN: Intersection_Node_Surface_Mesh - Exactly two surface mesh nodes must be supplied." );
    }

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::compute_global_coordinates()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - compute_global_coordinates() not implemented yet." );
        return { { std::numeric_limits< double >::quiet_NaN() } };
    }

    bool
    Intersection_Node_Surface_Mesh::determine_is_intersected(
            const Element_Intersection_Type aAncestorBasisFunction,
            const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
            const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates )
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - determine_is_intersected() not implemented yet." );
        return std::numeric_limits< double >::quiet_NaN();
    }

    Matrix< DDRMat > compute_raycast_rotation(
            const Matrix< DDRMat >& aFirstParentNodeGlobalCoordinates,
            const Matrix< DDRMat >& aSecondParentNodeGlobalCoordinates )
    {
        
    }

    real
    Intersection_Node_Surface_Mesh::compute_local_coordinate(
            const Matrix< DDRMat >&         aFirstParentNodeCoordinates,
            const Matrix< DDRMat >&         aSecondParentNodeCoordinates,
            const Cell< Matrix< DDRMat > >& aAncestorNodeCoordinates )
    {

    }

    void
    Intersection_Node_Surface_Mesh::append_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor )
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dcoordinate_dadv() not implemented yet." );

        return;
    };

    Matrix< DDSMat >
    Intersection_Node_Surface_Mesh::get_coordinate_determining_adv_ids()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_coordinate_determining_adv_ids() not implemented yet." );
        return { { -1 } };
    }

    real
    Intersection_Node_Surface_Mesh::get_dxi_dfield_from_ancestor( uint aAncestorIndex )
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dxi_dfield_from_ancestor() not implemented yet." );

        return std::numeric_limits< double >::quiet_NaN();
    }

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::get_dxi_dcoordinate_first_parent()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dxi_dcoordinate_first_parent() not implemented yet." );
        return { { std::numeric_limits< double >::quiet_NaN() } };
    }

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::get_dxi_dcoordinate_second_parent()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dxi_dcoordinate_second_parent() not implemented yet." );
        return { { std::numeric_limits< double >::quiet_NaN() } };
    }
}    // namespace moris::ge
