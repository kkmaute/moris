/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Single_Grain.cpp
 *
 */

#include "cl_GEN_Voxel_Geometry.hpp"
#include "cl_GEN_Voxel_Input.hpp"

#include "cl_Ascii.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Voxel_Geometry::Voxel_Geometry(
            std::shared_ptr< Voxel_Input > aVoxelField,
            uint                           aIndex )
            : mVoxelInput( aVoxelField )
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
        // TODO
        MORIS_ERROR( false, "A voxel geometry currently cannot create intersection nodes." );
        return nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< mtk::Field > > get_mtk_fields()
    {
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

}
