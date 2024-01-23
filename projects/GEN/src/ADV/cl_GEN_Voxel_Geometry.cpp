/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Voxel_Geometry.cpp
 *
 */

#include "cl_GEN_Voxel_Geometry.hpp"
#include "cl_GEN_Voxel_Input.hpp"
#include "cl_GEN_Intersection_Node_Voxel.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Voxel_Geometry::Voxel_Geometry(
            std::shared_ptr< Voxel_Input > aVoxelInput,
            uint                           aIndex )
            : Geometry( Design_Parameters(), 1e-12 )
            , mVoxelInput( aVoxelInput )
            , mIndex( aIndex )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Voxel_Geometry::depends_on_advs()
    {
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Voxel_Geometry::get_geometric_region(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        if ( mVoxelInput->get_voxel_ID( aNodeCoordinates ) == mIndex )
        {
            return Geometric_Region::POSITIVE;
        }
        else
        {
            return Geometric_Region::NEGATIVE;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node* Voxel_Geometry::create_intersection_node(
            uint                     aNodeIndex,
            const Cell< Node* >&     aBaseNodes,
            const Parent_Node&       aFirstParentNode,
            const Parent_Node&       aSecondParentNode,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder )
    {
        return new Intersection_Node_Voxel(
                aNodeIndex,
                aBaseNodes,
                aFirstParentNode,
                aSecondParentNode,
                aBackgroundGeometryType,
                aBackgroundInterpolationOrder,
                shared_from_this() );
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< mtk::Field > > get_mtk_fields()
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

}
