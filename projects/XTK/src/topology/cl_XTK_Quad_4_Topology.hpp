/*
 * cl_XTK_Quad_4_Topology.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef SRC_TOPOLOGY_CL_XTK_QUAD_4_TOPOLOGY_HPP_
#define SRC_TOPOLOGY_CL_XTK_QUAD_4_TOPOLOGY_HPP_


#include "topology/cl_XTK_Topology.hpp"
#include "topology/cl_XTK_Basis_Function.hpp"

// Basis Functions
#include "topology/cl_XTK_Quad_4_Basis_Function.hpp"

#include"assert/fn_xtk_assert.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Quad_4_Topology : public Topology<Real, Integer, Real_Matrix, Integer_Matrix>
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

    Basis_Function<Real,Real_Matrix> const & get_basis_function() const
    {
        return mBasisFunction;
    }

    void set_node_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
    {
        XTK_ASSERT(aNodeIndices.numel()==4,"Should be 4 associated with a quad 4 topology");
        mNodeIndices = aNodeIndices.copy();
    }

    std::shared_ptr<Topology<Real, Integer, Real_Matrix, Integer_Matrix>> copy() const
    {
        std::shared_ptr<Topology<Real, Integer, Real_Matrix, Integer_Matrix>> tTopologyCopy;
        tTopologyCopy = std::make_shared<Quad_4_Topology<Real, Integer, Real_Matrix, Integer_Matrix>>(mNodeIndices);
        return tTopologyCopy;
    }

private:
    moris::Matrix< moris::IndexMat > mNodeIndices;
    Quad_4_Basis_Function<Real,Real_Matrix> mBasisFunction;
};
}


#endif /* SRC_TOPOLOGY_CL_XTK_QUAD_4_TOPOLOGY_HPP_ */
