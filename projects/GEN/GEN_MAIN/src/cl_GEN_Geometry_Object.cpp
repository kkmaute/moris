/*
 * cl_GEN_Geometry_Object.cpp
 *
 *  Created on: Dec 19, 2019
 *      Author: sonne
 */

#include "cl_GEN_Geometry_Object.hpp"

namespace moris
{
    namespace ge
    {

    //------------------------------------------------------------------------------
    void GEN_Geometry_Object::set_parent_entity_index( moris::moris_index aEntityIndex )
    {
        mParentEntityIndex = aEntityIndex;
    }

    moris::moris_index const & GEN_Geometry_Object::get_parent_entity_index()
    {
        return mParentEntityIndex;
    }

    //------------------------------------------------------------------------------

    void GEN_Geometry_Object::set_interface_loc_coord( moris::real const & aLclCoord )
    {
        mInterfaceLclCoords = aLclCoord;
    }

    moris::real const & GEN_Geometry_Object::get_interface_lcl_coord()
    {
        return mInterfaceLclCoords;
    }
    //------------------------------------------------------------------------------

    void GEN_Geometry_Object::set_interface_glb_coord( moris::Matrix< moris::DDRMat > const & aGlbCoord )
    {
        mInterfaceGlbCoords = aGlbCoord.copy();
    }

    moris::Matrix< moris::DDRMat > const & GEN_Geometry_Object::get_interface_glb_coord()
    {
        return mInterfaceGlbCoords;
    }

    void GEN_Geometry_Object::mark_all_nodes_as_on_interface( )
    {
        mAllParentNodesOnInterface = true;
        mHasParentNodesOnInterface = true;
    }
    //------------------------------------------------------------------------------

    void GEN_Geometry_Object::mark_node_as_on_interface( moris::moris_index aNodeOrdinal )
    {
        mHasParentNodesOnInterface = true;
    }

    //------------------------------------------------------------------------------

    bool GEN_Geometry_Object::all_parent_nodes_on_interface( )
    {
        return mAllParentNodesOnInterface;
    }
    //------------------------------------------------------------------------------

    bool GEN_Geometry_Object::has_parent_nodes_on_interface( )
    {
        return mHasParentNodesOnInterface;
    }
    //------------------------------------------------------------------------------

    }   // end ge namespace
}       // end moris namespace


