/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Surface_Mesh_Geometry.cpp
 *
 */

#include "cl_GEN_Surface_Mesh_Geometry.hpp"
#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Parameters::Surface_Mesh_Parameters( const ParameterList& aParameterList )
            : Field_Parameters( aParameterList )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Geometry::Surface_Mesh_Geometry(
            const std::string&      aFilePath,
            Surface_Mesh_Parameters aParameters )
            : Object( aFilePath, aParameters.mOffsets )
            , mParameters( aParameters )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Surface_Mesh_Geometry::get_geometric_region(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        sdf::Object_Region tRegion = raycast_point( *this, aNodeCoordinates );

        switch ( tRegion )
        {
            case sdf::Object_Region::INSIDE:
            {
                return Geometric_Region::NEGATIVE;
                break;
            }
            case sdf::Object_Region::OUTSIDE:
            {
                return Geometric_Region::POSITIVE;
                break;
            }
            default:
            {
                return Geometric_Region::INTERFACE;
                break;
            }
        }
        
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< Intersection_Node > Surface_Mesh_Geometry::create_intersection_node(
            uint                                 aEdgeFirstNodeIndex,
            uint                                 aEdgeSecondNodeIndex,
            std::shared_ptr< Intersection_Node > aEdgeFirstIntersectionNode,
            std::shared_ptr< Intersection_Node > aEdgeSecondIntersectionNode,
            const Matrix< DDRMat >&              aEdgeFirstNodeLocalCoordinates,
            const Matrix< DDRMat >&              aEdgeSecondNodeLocalCoordinates,
            const Matrix< DDRMat >&              aEdgeFirstNodeGlobalCoordinates,
            const Matrix< DDRMat >&              aEdgeSecondNodeGlobalCoordinates,
            const Matrix< DDUMat >&              aBackgroundElementNodeIndices,
            const Cell< Matrix< DDRMat > >&      aBackgroundElementNodeCoordinates )
    {
        // Create linear intersection node
        return std::make_shared< Intersection_Node_Surface_Mesh >(
                aEdgeFirstIntersectionNode,
                aEdgeSecondIntersectionNode,
                aEdgeFirstNodeIndex,
                aEdgeSecondNodeIndex,
                aEdgeFirstNodeGlobalCoordinates,
                aEdgeSecondNodeGlobalCoordinates,
                shared_from_this() );
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< mtk::Field > > Surface_Mesh_Geometry::get_mtk_fields()
    {
        // TODO: maybe?
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::ge
