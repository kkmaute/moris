/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Edge_Topology.hpp
 *
 */

#ifndef PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_EDGE_TOPOLOGY_HPP_
#define PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_EDGE_TOPOLOGY_HPP_

#include "../../new/additional/cl_GEN_Basis_Function.hpp"
#include "../../new/additional/cl_GEN_Linear_Basis_Functions.hpp"
#include "../../new/additional/cl_GEN_Topology.hpp"

namespace moris
{
namespace ge
{
class Edge_Topology : public Topology
{

public:
    Edge_Topology()
    {

    }

    Edge_Topology(moris::Matrix< moris::IndexMat > const & aNodeIndices)

    {
        this->set_node_indices(aNodeIndices);

    }

    // Required Interface Functions
    enum Topology_Type get_topology_type() const
    {
        return Topology_Type::EDGE;
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
        MORIS_ASSERT(aNodeIndices.n_cols()==2,"Should be 2 associated with a edge topology");

        mNodeIndices = aNodeIndices.copy();
    }

    std::shared_ptr<Topology> copy() const
    {
        std::shared_ptr<Topology> tTopologyCopy;
        tTopologyCopy = std::make_shared<Edge_Topology>(mNodeIndices);
        return tTopologyCopy;
    }
private:
    moris::Matrix< moris::IndexMat > mNodeIndices;
    Linear_Basis_Function mBasisFunction;

};
}
}

#endif /* PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_EDGE_TOPOLOGY_HPP_ */

