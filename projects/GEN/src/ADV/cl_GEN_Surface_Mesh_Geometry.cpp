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
            : mFilePath( aParameterList.get( "file_path" ) )
    {
        string_to_cell_mat( aParameterList.get< std::string >( "offsets" ), mOffsets );
    }

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Geometry::Surface_Mesh_Geometry( Surface_Mesh_Parameters aParameters )
            : Object( aParameters.mFilePath, aParameters.mOffsets )
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

    Intersection_Node* Surface_Mesh_Geometry::create_intersection_node(
            uint                 aNodeIndex,
            const Cell< Node* >& aBaseNodes,
            const Parent_Node&   aFirstParentNode,
            const Parent_Node&   aSecondParentNode,
            mtk::Geometry_Type   aBaseGeometryType )
    {
        // Create linear intersection node
        return std::make_shared< Intersection_Node_Surface_Mesh >(
                
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
