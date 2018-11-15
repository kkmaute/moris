/*
 * cl_Geometry_Object.hpp
 *
 *  Created on: Jun 21, 2017
 *      Author: ktdoble
 */

#ifndef SRC_GEOMENG_CL_MGE_GEOMETRY_OBJECT_HPP_
#define SRC_GEOMENG_CL_MGE_GEOMETRY_OBJECT_HPP_

// Standard library includes
#include <memory> // for shared_ptr

// Matrix Include
#include "linalg/cl_XTK_Matrix.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Geometry_Object
{
public:
    Geometry_Object():
        mAllParentNodesOnInterface(false),
        mHasParentNodesOnInterface(false)
    {
    }


    Geometry_Object(Integer aParentEntityIndex)
    {
        mParentEntityIndex = aParentEntityIndex;
    }

    ~Geometry_Object()
    {

    }


    void
    set_phase_val_row(Integer aPhaseValRowIndex)
    {
        mPhaseValIndex = aPhaseValRowIndex;
    }

    Integer
    get_phase_val_row() const
    {
        return mPhaseValIndex;
    }


    /**
     * This tells the geometry object which entity index it is associated with. At this point, the dimension of this parent
     * entity is up to the user to keep track of.
     * @param[in] aEntityIndex - Parent entity index
     */
    void set_parent_entity_index(Integer aEntityIndex)
    {
        mParentEntityIndex = aEntityIndex;
    }

    Integer const & get_parent_entity_index()
    {
        return mParentEntityIndex;
    }

    /** Currently set_interface_lcl_coord is only needed for an edge and requires only 1 value,
     * In future  want to extend to different dimension entities
     */
    void set_interface_loc_coord(Real const & aLclCoord)
    {
        mInterfaceLclCoords = aLclCoord;
    }

    /**
     * Global coordinate of interface point
     */
    void set_interface_glb_coord(moris::Matrix< Real_Matrix > const & aGlbCoord)
    {
        mInterfaceGlbCoords = aGlbCoord.copy();
    }

    /**
     * Sensitivity with respect to design relevant design variables hosted in the geometry engine
     */
    void set_sensitivity_dx_dp(moris::Matrix< Real_Matrix > const & aSensitivitydxdp)
    {
        mSensitivityDxDp = aSensitivitydxdp.copy();
    }

    moris::Matrix< Real_Matrix > const &
    get_sensitivity_dx_dp() const
    {
        return mSensitivityDxDp;
    }

    void set_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeADVIndices)
    {
        mNodeADVIndices = aNodeADVIndices.copy();
    }

    moris::Matrix< moris::IndexMat > const &
    get_node_adv_indices() const
    {
        return mNodeADVIndices;
    }

    Real const & get_interface_lcl_coord()
    {
        return mInterfaceLclCoords;
    }

    moris::Matrix< Real_Matrix > const &
    get_interface_glb_coord()
    {
        return mInterfaceGlbCoords;
    }

    void
    mark_all_nodes_as_on_interface()
    {
        mAllParentNodesOnInterface = true;
        mHasParentNodesOnInterface = true;
    }

    //
    void
    mark_node_as_on_interface(moris::moris_index aNodeOrdinal)
    {
        mNodesOnInterface.push_back(aNodeOrdinal);
        mHasParentNodesOnInterface = true;
    }

    void
    mark_nodes_as_not_on_interface()
    {
        mHasParentNodesOnInterface = false;
    }

    bool
    all_parent_nodes_on_interface()
    {
        return mAllParentNodesOnInterface;
    }

    bool
    has_parent_nodes_on_interface()
    {
        return mHasParentNodesOnInterface;
    }



private:
    Integer mPhaseValIndex;

    Real                             mInterfaceLclCoords;
    Integer                          mParentEntityIndex;
    moris::Matrix< Real_Matrix >     mSensitivityDxDp;
    moris::Matrix< moris::IndexMat > mNodeADVIndices;
    moris::Matrix< Real_Matrix >     mInterfaceGlbCoords;

    // Information about coincidence (along an edge)
    bool                     mAllParentNodesOnInterface;
    bool                     mHasParentNodesOnInterface;
    Cell<moris::moris_index> mNodesOnInterface;
};
}


#endif /* SRC_GEOMENG_CL_MGE_GEOMETRY_OBJECT_HPP_ */
