/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Hexahedron_8_Topology.hpp
 *
 */

#ifndef SRC_TOPOLOGY_CL_XTK_HEXAHEDRON_8_TOPOLOGY_HPP_
#define SRC_TOPOLOGY_CL_XTK_HEXAHEDRON_8_TOPOLOGY_HPP_

#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Basis_Function.hpp"
#include "fn_isvector.hpp"

// Basis Functions
#include "cl_XTK_Hexahedron_8_Basis_Function.hpp"

namespace xtk
{
class Hexahedron_8_Topology : public Topology
{

public:
    Hexahedron_8_Topology()
    {

    }

    Hexahedron_8_Topology(moris::Matrix< moris::IndexMat > const & aNodeIndices)
    {
        this->set_node_indices(aNodeIndices);
    }

    // Required Interface Functions
    enum Topology_Type get_topology_type() const
    {
        return Topology_Type::HEXA_8;
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
        MORIS_ASSERT(aNodeIndices.numel()>=8 && moris::isvector(aNodeIndices),"Should be 8 associated with a HEX8 topology");

        mNodeIndices.resize(1,8);

        mNodeIndices(0) = aNodeIndices(0);
        mNodeIndices(1) = aNodeIndices(1);
        mNodeIndices(2) = aNodeIndices(2);
        mNodeIndices(3) = aNodeIndices(3);
        mNodeIndices(4) = aNodeIndices(4);
        mNodeIndices(5) = aNodeIndices(5);
        mNodeIndices(6) = aNodeIndices(6);
        mNodeIndices(7) = aNodeIndices(7);
    }

    std::shared_ptr<Topology> copy() const
    {
        std::shared_ptr<Topology> tTopologyCopy;
        tTopologyCopy = std::make_shared<Hexahedron_8_Topology>(mNodeIndices);
        return tTopologyCopy;
    }
private:
    moris::Matrix< moris::IndexMat > mNodeIndices;
    Hexahedron_8_Basis_Function mBasisFunction;

};
}

#endif /* SRC_TOPOLOGY_CL_XTK_HEXAHEDRON_8_TOPOLOGY_HPP_ */

