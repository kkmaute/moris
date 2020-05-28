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
    void GEN_Geometry_Object::set_phase_val_row( moris::moris_index aPhaseValRowIndex )
    {
        mPhaseValIndex = aPhaseValRowIndex;
    }

    moris::moris_index GEN_Geometry_Object::get_phase_val_row() const
    {
        return mPhaseValIndex;
    }

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
    //------------------------------------------------------------------------------

    void GEN_Geometry_Object::set_sensitivity_dx_dp( moris::Matrix< moris::DDRMat > const & aSensitivitydxdp )
    {
        mSensitivityDxDp = aSensitivitydxdp.copy();
    }

    moris::Matrix< moris::DDRMat > const & GEN_Geometry_Object::get_sensitivity_dx_dp() const
    {
        return mSensitivityDxDp;
    }
    //------------------------------------------------------------------------------

    void GEN_Geometry_Object::set_node_adv_indices( moris::Matrix< moris::IndexMat > const & aNodeADVIndices )
    {
        mNodeADVIndices = aNodeADVIndices.copy();
    }

    moris::Matrix< moris::IndexMat > const & GEN_Geometry_Object::get_node_adv_indices() const
    {
        return mNodeADVIndices;
    }
    //------------------------------------------------------------------------------

    void GEN_Geometry_Object::mark_all_nodes_as_on_interface( )
    {
        mAllParentNodesOnInterface = true;
        mHasParentNodesOnInterface = true;
    }
    //------------------------------------------------------------------------------

    void GEN_Geometry_Object::mark_node_as_on_interface( moris::moris_index aNodeOrdinal )
    {
        mNodesOnInterface.push_back(aNodeOrdinal);
        mHasParentNodesOnInterface = true;
    }
    //------------------------------------------------------------------------------

    void GEN_Geometry_Object::mark_nodes_as_not_on_interface( )
    {
        mHasParentNodesOnInterface = false;
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

    void GEN_Geometry_Object::set_parent_entity_topology( std::shared_ptr<xtk::Topology> tParentTopo )
    {
        MORIS_ASSERT(mParentTopology == nullptr,
                     "Geometry object parent entity topology has already been set");
        mParentTopology = tParentTopo;
    }
    //------------------------------------------------------------------------------

    xtk::Topology const & GEN_Geometry_Object::get_parent_entity_topology( )
    {
        MORIS_ASSERT( mParentTopology != nullptr,
                     "Geometry object parent entity topology has not been set, either this is not an interface geometry object or set_parent_entity_topology was not called" );

        return (*mParentTopology);
    }

    //------------------------------------------------------------------------------

    }   // end ge namespace
}       // end moris namespace


