/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Voxel.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Voxel.hpp"
#include "cl_GEN_Voxel_Geometry.hpp"
#include "cl_GEN_Parent_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node_Voxel::Intersection_Node_Voxel(
            uint                              aNodeIndex,
            const Cell< Node* >&              aBaseNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder,
            std::shared_ptr< Voxel_Geometry > aInterfaceGeometry )
            : Intersection_Node(
                    aNodeIndex,
                    aBaseNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    Intersection_Node_Voxel::get_local_coordinate( aFirstParentNode, aSecondParentNode, aInterfaceGeometry ),
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    aInterfaceGeometry )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    real Intersection_Node_Voxel::get_local_coordinate(
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode,
            std::shared_ptr< Voxel_Geometry > aInterfaceGeometry )
    {
        // Get parent geometric regions
        Geometric_Region tFirstParentGeometricRegion = aInterfaceGeometry->get_geometric_region( aFirstParentNode.get_index(), aFirstParentNode.get_global_coordinates() );
        Geometric_Region tSecondParentGeometricRegion = aInterfaceGeometry->get_geometric_region( aSecondParentNode.get_index(), aSecondParentNode.get_global_coordinates() );

        // Return local coordinate based on geometric regions of the parent nodes
        if ( tFirstParentGeometricRegion == tSecondParentGeometricRegion and tFirstParentGeometricRegion not_eq Geometric_Region::INTERFACE )
        {
            // Geometric regions are the same; no intersection
            return MORIS_REAL_MAX;
        }
        else if ( tFirstParentGeometricRegion == Geometric_Region::INTERFACE )
        {
            // First parent on interface
            return -1.0;
        }
        else if ( tSecondParentGeometricRegion == Geometric_Region::INTERFACE )
        {
            // Second parent on interface
            return 1.0;
        }
        else
        {
            // Interface is midway between the parent nodes
            return 0.0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
