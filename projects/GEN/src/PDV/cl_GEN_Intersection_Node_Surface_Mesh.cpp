/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Surface_Mesh.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"
#include "cl_GEN_Parent_Node.hpp"

#include "fn_norm.hpp"

namespace moris::ge
{
    Intersection_Node_Surface_Mesh::Intersection_Node_Surface_Mesh(
            uint                     aNodeIndex,
            const Cell< Node* >&     aBackgroundNodes,
            const Parent_Node&       aFirstParentNode,
            const Parent_Node&       aSecondParentNode,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder,
            Surface_Mesh_Geometry&   aInterfaceGeometry )
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

    Geometry& Intersection_Node_Surface_Mesh::get_interface_geometry()
    {
        
        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Geometry& Intersection_Node_Surface_Mesh::get_interface_geometry() const
    {
        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Intersection_Node_Surface_Mesh::append_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dcoordinate_dadv() not implemented yet." );
    }
    
    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Intersection_Node_Surface_Mesh::get_coordinate_determining_adv_ids() const
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_coordinate_determining_adv_ids() not implemented yet." );
        return { { -1 } };
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
}    // namespace moris::ge
