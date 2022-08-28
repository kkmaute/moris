/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Quad_4_Topology.hpp
 *
 */

#ifndef SRC_TOPOLOGY_CL_XTK_QUAD_4_TOPOLOGY_HPP_
#define SRC_TOPOLOGY_CL_XTK_QUAD_4_TOPOLOGY_HPP_

#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Basis_Function.hpp"

// Basis Functions
#include "cl_XTK_Quad_4_Basis_Function.hpp"

namespace xtk
{
class Quad_4_Topology : public Topology
{

public:
    Quad_4_Topology()
    {

    }

    Quad_4_Topology(moris::Matrix< moris::IndexMat > const & aNodeIndices)

    {
        this->set_node_indices(aNodeIndices);

    }

    // Required Interface Functions
    enum Topology_Type get_topology_type() const
    {
        return Topology_Type::QUAD_4;
    }

    moris::Matrix< moris::IndexMat > const & get_node_indices() const
    {
        return mNodeIndices;
    }

    Basis_Function const & get_basis_function() const
    {
        return mBasisFunction;
    }

    void set_node_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
    {
        mNodeIndices.resize(1,4);

        mNodeIndices(0) = aNodeIndices(0);
        mNodeIndices(1) = aNodeIndices(1);
        mNodeIndices(2) = aNodeIndices(2);
        mNodeIndices(3) = aNodeIndices(3);
    }

    std::shared_ptr<Topology> copy() const
    {
        std::shared_ptr<Topology> tTopologyCopy;
        tTopologyCopy = std::make_shared<Quad_4_Topology>(mNodeIndices);
        return tTopologyCopy;
    }

private:
    moris::Matrix< moris::IndexMat > mNodeIndices;
    Quad_4_Basis_Function mBasisFunction;
};
}

#endif /* SRC_TOPOLOGY_CL_XTK_QUAD_4_TOPOLOGY_HPP_ */

