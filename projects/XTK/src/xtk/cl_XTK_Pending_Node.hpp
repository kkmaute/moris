/*
 * cl_XTK_Pending_Node.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_PENDING_NODE_HPP_
#define SRC_XTK_CL_XTK_PENDING_NODE_HPP_

#include<memory>

#include "linalg/cl_XTK_Matrix_Base.hpp"
// XTKL: Linear Algebra Includes
#include "topology/cl_XTK_Topology.hpp"


#include "mesh/cl_Mesh_Enums.hpp"
#include "containers/cl_XTK_Cell.hpp"
#include "assert/fn_xtk_assert.hpp"
namespace xtk
{

template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
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


    void set_pending_node_info(Integer* aNodeIndLocation,
                               Integer* aNodeIdLocation,
                               Mat<Real, Real_Matrix> const & aCoordinates,
                               Topology<Real,Integer,Real_Matrix,Integer_Matrix> const & aParentTopology,
                               Mat<Real, Real_Matrix> const & aLocalCoordinates)
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
    Mat<Real, Real_Matrix> const &
    get_coordinates() const
    {
        return mCoordinates;
    }

    Integer get_node_index() const
    {
        return *mNodeInd;
    }

    Integer get_node_id() const
    {
        return *mNodeId;
    }

    // Parent Related Functions
    Topology<Real,Integer,Real_Matrix,Integer_Matrix> const &
    get_parent_topology()
    {
        return (*mParentTopology);
    }

    Mat<Real, Real_Matrix> const &
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

    void add_field_values(Mat<Real, Real_Matrix> const & aFieldValues)
    {
        //XTK_ASSERT(!mHasFields,"Field values have already been defined on this pending node (please add all fields on a pending node at one time)");
        mFieldValues = aFieldValues.copy();
        mHasFields = true;
    }

    Mat<Real, Real_Matrix> const & get_field_data() const
    {
        //XTK_ASSERT(mHasFields,"This pending node does not have a field");
        return mFieldValues;
    }


    /**
     * Sensitivity with respect to design variables
     */
    void set_sensitivity_dx_dp(Mat<Real, Real_Matrix> const & aSensitivitydxdp)
    {
        XTK_ASSERT(!mHasDxDp,"This pending node already has dxdp information");
        mSensitivityDxDp = aSensitivitydxdp.copy();
        mHasDxDp = true;
    }

    Mat<Real, Real_Matrix> const &
    get_sensitivity_dx_dp()
    {
        XTK_ASSERT(mHasDxDp,"This pending node does not have dxdp information");
        return mSensitivityDxDp;
    }

    void set_node_adv_indices(Mat<Integer, Integer_Matrix> const & aNodeADVIndices)
    {
        XTK_ASSERT(!mSparseDxDp,"This pending node already has sparse dxdp information");

        mNodeADVIndices = aNodeADVIndices.copy();
        mSparseDxDp = true;
    }

    Mat<Integer, Integer_Matrix> const &
    get_node_adv_indices()
    {
        XTK_ASSERT(mSparseDxDp,"This pending node does not have sparse dxdp information");
        return mNodeADVIndices;
    }


private:
    Integer* mNodeId;
    Integer* mNodeInd;
    Mat<Real, Real_Matrix> mCoordinates;

    //Parent Entity information
    std::shared_ptr< xtk::Topology<Real,Integer,Real_Matrix,Integer_Matrix>> mParentTopology;
    Mat<Real, Real_Matrix> mLocalCoordinates;

    // Field information
    bool mHasFields;
    Mat<Real, Real_Matrix> mFieldValues;

    // has dxdp
    bool mHasDxDp;
    bool mSparseDxDp;
    Mat<Real, Real_Matrix>      mSensitivityDxDp;
    Mat<Integer,Integer_Matrix> mNodeADVIndices;
};

}

#endif /* SRC_XTK_CL_XTK_PENDING_NODE_HPP_ */
