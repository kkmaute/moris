/*
 * cl_Geometry_Object.hpp
 *
 *  Created on: Jun 21, 2017
 *      Author: ktdoble
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMENG_CL_GEN_GEOMETRY_OBJECT_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMENG_CL_GEN_GEOMETRY_OBJECT_HPP_

// Standard library includes
#include <memory> // for shared_ptr

#include "cl_Matrix.hpp"

// GE includes
#include "cl_GEN_Property.hpp"
#include "cl_GEN_Pdv_Enums.hpp"

// XTK includes
#include "cl_XTK_Topology.hpp"

namespace moris
{
namespace ge
{
class GEN_Geometry_Object
{
public:
    GEN_Geometry_Object():
        mAllParentNodesOnInterface( false ),
        mHasParentNodesOnInterface( false )
    {    }


    GEN_Geometry_Object( moris::moris_index aParentEntityIndex )
    {
        mParentEntityIndex = aParentEntityIndex;
    }

    ~GEN_Geometry_Object()
    {    }

    //------------------------------------------------------------------------------
    void set_phase_val_row( moris::moris_index aPhaseValRowIndex );

    moris::moris_index get_phase_val_row() const;
    //------------------------------------------------------------------------------

    /**
     * This tells the geometry object which entity index it is associated with. At this point, the dimension of this parent
     * entity is up to the user to keep track of.
     * @param[in] aEntityIndex - Parent entity index
     */
    void set_parent_entity_index( moris::moris_index aEntityIndex );

    moris::moris_index const & get_parent_entity_index();
    //------------------------------------------------------------------------------

    /**
     * Currently set_interface_lcl_coord is only needed for an edge and requires only 1 value,
     * In future  want to extend to different dimension entities
     */
    void set_interface_loc_coord( moris::real const & aLclCoord );

    moris::real const & get_interface_lcl_coord( );
    //------------------------------------------------------------------------------

    /**
     * Global coordinate of interface point
     */
    void set_interface_glb_coord( moris::Matrix< moris::DDRMat > const & aGlbCoord );

    moris::Matrix< moris::DDRMat > const & get_interface_glb_coord( );
    //------------------------------------------------------------------------------

    /**
     * Sensitivity with respect to relevant design variables hosted in the geometry engine
     */
    void set_sensitivity_dx_dp( moris::Matrix< moris::DDRMat > const & aSensitivitydxdp );

    moris::Matrix< moris::DDRMat > const & get_sensitivity_dx_dp() const;
    //------------------------------------------------------------------------------

    void set_node_adv_indices( moris::Matrix< moris::IndexMat > const & aNodeADVIndices );

    moris::Matrix< moris::IndexMat > const & get_node_adv_indices() const;
    //------------------------------------------------------------------------------

    void mark_all_nodes_as_on_interface( );
    //------------------------------------------------------------------------------

    void mark_node_as_on_interface( moris::moris_index aNodeOrdinal );
    //------------------------------------------------------------------------------

    void mark_nodes_as_not_on_interface( );
    //------------------------------------------------------------------------------

    bool all_parent_nodes_on_interface( );
    //------------------------------------------------------------------------------

    bool has_parent_nodes_on_interface( );
    //------------------------------------------------------------------------------

    void set_parent_entity_topology( std::shared_ptr<xtk::Topology> tParentTopo );
    //------------------------------------------------------------------------------

    xtk::Topology const & get_parent_entity_topology( );
    //------------------------------------------------------------------------------
    void register_pdv_type( enum PDV_Type aPdvType );

    void set_pdv_list( moris::Cell< enum PDV_Type > aPdvList );

    moris::Cell< enum PDV_Type > get_pdv_list( );

    //------------------------------------------------------------------------------

    real get_pdv_val( const enum PDV_Type aPdvType );


    //------------------------------------------------------------------------------
private:
    moris::moris_index mPhaseValIndex;
    //------------------------------------------------------------------------------
    moris::real                      mInterfaceLclCoords;
    moris::moris_index               mParentEntityIndex;
    moris::Matrix< moris::DDRMat >   mSensitivityDxDp;
    moris::Matrix< moris::IndexMat > mNodeADVIndices;
    moris::Matrix< moris::DDRMat >   mInterfaceGlbCoords;

    //------------------------------------------------------------------------------
    // Parent topology
    std::shared_ptr< xtk::Topology > mParentTopology = nullptr;
    //------------------------------------------------------------------------------
    // Information about coincidence (along an edge)
    bool                            mAllParentNodesOnInterface;
    bool                            mHasParentNodesOnInterface;
    moris::Cell<moris::moris_index> mNodesOnInterface;
};
}
}

#endif /* PROJECTS_GEN_SRC_NEW_GEOMENG_CL_GEN_GEOMETRY_OBJECT_HPP_ */
