/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Pending_Node.hpp
 *
 */

#ifndef SRC_XTK_CL_XTK_PENDING_NODE_HPP_
#define SRC_XTK_CL_XTK_PENDING_NODE_HPP_

#include<memory>

#include "cl_Matrix.hpp"
// XTKL: Linear Algebra Includes
#include "cl_XTK_Topology.hpp"

#include "cl_MTK_Enums.hpp"
#include "cl_Vector.hpp"

namespace xtk
{

class Pending_Node
{
public:
    Pending_Node()
    {
        mHasFields = false;
        mHasDxDp = false;
        mSparseDxDp = false;
    }

    ~Pending_Node()
    {
    }

    void set_pending_node_info(moris::moris_index* aNodeIndLocation,
                               moris::moris_index* aNodeIdLocation,
                               moris::Matrix< moris::DDRMat > const & aCoordinates,
                               Topology const & aParentTopology,
                               moris::Matrix< moris::DDRMat > const & aLocalCoordinates)
    {
        mHasDxDp = false;
        mSparseDxDp = false;
        mNodeInd = aNodeIndLocation;
        mNodeId = aNodeIdLocation;
        mCoordinates = aCoordinates.copy();
        mLocalCoordinates = aLocalCoordinates.copy();
        mParentTopology = aParentTopology.copy();
    }

    // Node itself Functions
    moris::Matrix< moris::DDRMat > const &
    get_coordinates() const
    {
        return mCoordinates;
    }

    moris::moris_index get_node_index() const
    {
        return *mNodeInd;
    }

    moris::moris_index get_node_id() const
    {
        return *mNodeId;
    }

    // Parent Related Functions
    Topology const &
    get_parent_topology()
    {
        return (*mParentTopology);
    }

    std::shared_ptr< Topology >
    get_parent_topology_ptr()
    {
        return (mParentTopology);
    }

    moris::Matrix< moris::DDRMat > const &
    get_local_coordinate_relative_to_parent()
    {
        return mLocalCoordinates;
    }

    // Field related functions
    bool
    has_fields() const
    {
        return mHasFields;
    }

    void add_field_values(moris::Matrix< moris::DDRMat > const & aFieldValues)
    {
        //MORIS_ASSERT(!mHasFields,"Field values have already been defined on this pending node (please add all fields on a pending node at one time)");
        mFieldValues = aFieldValues.copy();
        mHasFields = true;
    }

    moris::Matrix< moris::DDRMat > const & get_field_data() const
    {
        //MORIS_ASSERT(mHasFields,"This pending node does not have a field");
        return mFieldValues;
    }

    /**
     * Sensitivity with respect to design variables
     */
    void set_sensitivity_dx_dp(moris::Matrix< moris::DDRMat > const & aSensitivitydxdp)
    {
        MORIS_ASSERT(!mHasDxDp,"This pending node already has dxdp information");
        mSensitivityDxDp = aSensitivitydxdp.copy();
        mHasDxDp = true;
    }

    moris::Matrix< moris::DDRMat > const &
    get_sensitivity_dx_dp()
    {
        MORIS_ASSERT(mHasDxDp,"This pending node does not have dxdp information");
        return mSensitivityDxDp;
    }

    void set_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeADVIndices)
    {
        MORIS_ASSERT(!mSparseDxDp,"This pending node already has sparse dxdp information");

        mNodeADVIndices = aNodeADVIndices.copy();
        mSparseDxDp = true;
    }

    moris::Matrix< moris::IndexMat > const &
    get_node_adv_indices()
    {
        MORIS_ASSERT(mSparseDxDp,"This pending node does not have sparse dxdp information");
        return mNodeADVIndices;
    }

private:
    moris::moris_index* mNodeId;
    moris::moris_index* mNodeInd;
    moris::Matrix< moris::DDRMat > mCoordinates;

    //Parent Entity information
    std::shared_ptr< Topology > mParentTopology;
    moris::Matrix< moris::DDRMat > mLocalCoordinates;

    // Field information
    bool mHasFields;
    moris::Matrix< moris::DDRMat > mFieldValues;

    // has dxdp
    bool mHasDxDp;
    bool mSparseDxDp;
    moris::Matrix< moris::DDRMat >   mSensitivityDxDp;
    moris::Matrix< moris::IndexMat > mNodeADVIndices;
};

}

#endif /* SRC_XTK_CL_XTK_PENDING_NODE_HPP_ */

