/*
 * cl_XTK_Node.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_NODE_HPP_
#define SRC_XTK_CL_XTK_NODE_HPP_

// Standard includes
#include <limits>

// XTKL: Logging and Assertion Includes
#include "assert/fn_xtk_assert.hpp"
#include "mesh/cl_Mesh_Enums.hpp"

#include "cl_Matrix.hpp"
#include "geomeng/cl_MGE_Geometry_Object.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Node
{
public:
    Node() :
            mInterfaceFlag(false),
            mNodeInd(std::numeric_limits<Integer>::max()),
            mNodeId(std::numeric_limits<Integer>::max()),
            mDxDp(NULL)
    {
    }



    ~Node()
    {

    }

    /**
     * Set global node identifier for this node
     * The global node identifier lives in the parent mesh (mesh)
     */
    void set_node_id(Integer const & aNodeId)
    {
//        XTK_ASSERT(mNodeId == std::numeric_limits<Integer>::max(), "Node Id has already been set");
        mNodeId = aNodeId;
    }

    /**
     * NOTE: At this level, the global node ID is used for sorting and selecting tabulated connectivities for node hierarchical subdivision
     * @return Global node identifier
     */
    Integer const & get_node_id() const
    {
        if(mNodeId == std::numeric_limits<Integer>::max())
        {
            XTK_ERROR<< "Node Id has not been set";
        }
        return mNodeId;
    }

    /**
     * Set the node processor local identifier (contiguous)
     * @param aNodeInd - Node index
     */
    void set_node_index(Integer const & aNodeInd)
    {
        XTK_ASSERT(mNodeInd == std::numeric_limits<Integer>::max(), "Node Index has already been set");
        mNodeInd = aNodeInd;
    }

    /**
     * NOTE: In general, this should be used over node ids because of its contiguous behavior.
     * @return Nodes processor local identifier
     *TODO: MOVE TOWARDS USING IDs
     */
    Integer const & get_node_index() const
    {
        XTK_ASSERT(mNodeInd != std::numeric_limits<Integer>::max(), "Node Index has not been set");
        return mNodeInd;
    }

    void
    set_geometry_object(Geometry_Object<Real,Integer,Real_Matrix,Integer_Matrix> & aGeometryObject)
    {
        mpGeometryObject = std::make_shared<Geometry_Object<Real,Integer,Real_Matrix,Integer_Matrix>>(aGeometryObject);
    }



    /**
     * set_interface_node tells the node class that this node is located on an interface
     *
     */
    void set_interface_node()
    {
        XTK_ASSERT(!mInterfaceFlag, "Interface bool has already been set");
        mInterfaceFlag = true;
    }

    void add_dx_dp(moris::Matrix<Real,Real_Matrix> * aDxDp)
    {
        mDxDp = aDxDp;
    }

    moris::Matrix<Real,Real_Matrix> *
    get_dx_dp() const
    {
        XTK_ASSERT(mDxDp,"No Sensitivity Dx Dp has been set in this node, This is either a mistake or not an interface node");
        return mDxDp;
    }

    void add_adv_indices(moris::Matrix<Integer, Integer_Matrix> * aADVIndices)
    {
        mADVIndices = aADVIndices;
    }



    moris::Matrix<Integer, Integer_Matrix> *
    get_adv_indices() const
    {
        XTK_ASSERT(mADVIndices,"No ADV indices have been set in this node, This is either a mistake or not an interface node");
        return mADVIndices;
    }

    /**
     * Tells the asker if this node lives on the interface or not
     * @return
     */
    bool is_interface_node() const
    {
        return mInterfaceFlag;
    }

private:
    bool    mInterfaceFlag;
    Integer mNodeInd;
    Integer mNodeId;
    std::shared_ptr<Geometry_Object<Real, Integer, Real_Matrix,Integer_Matrix>> mpGeometryObject;

    // Remove the following because we access them but do not store them here
    Real    mFieldValue;
    Integer mParentEntityIndex;
    moris::Matrix<Real,Real_Matrix> * mDxDp;
    moris::Matrix<Integer,Integer_Matrix> * mADVIndices;
    enum EntityRank mParentEntityRank;
};
}


#endif /* SRC_XTK_CL_XTK_NODE_HPP_ */
