/*
 * cl_XTK_Hexahedron_8.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef SRC_TOPOLOGY_CL_XTK_HEXAHEDRON_8_TOPOLOGY_HPP_
#define SRC_TOPOLOGY_CL_XTK_HEXAHEDRON_8_TOPOLOGY_HPP_

#include "topology/cl_XTK_Topology.hpp"
#include "topology/cl_XTK_Basis_Function.hpp"

// Basis Functions
#include "topology/cl_XTK_Hexahedron_8_Basis_Function.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Hexahedron_8_Topology : public Topology<Real, Integer, Real_Matrix, Integer_Matrix>
{

public:
    Hexahedron_8_Topology()
    {

    }

    Hexahedron_8_Topology(moris::Matrix<Integer, Integer_Matrix> const & aNodeIndices)
    {
        this->set_node_indices(aNodeIndices);
    }

    // Required Interface Functions
    enum Topology_Type get_topology_type() const
    {
        return Topology_Type::HEXA_8;
    }
    moris::Matrix<Integer, Integer_Matrix> const & get_node_indices() const
    {
        return mNodeIndices;
    }

    Basis_Function<Real,Real_Matrix> const & get_basis_function() const
    {
        return mBasisFunction;
    }

    void set_node_indices(moris::Matrix<Integer, Integer_Matrix> const & aNodeIndices)
    {
        XTK_ASSERT(aNodeIndices.n_cols()==8,"Should be 8 associated with a HEX8 topology");
        mNodeIndices = aNodeIndices.copy();
    }


    std::shared_ptr<Topology<Real, Integer, Real_Matrix, Integer_Matrix>> copy() const
    {
        std::shared_ptr<Topology<Real, Integer, Real_Matrix, Integer_Matrix>> tTopologyCopy;
        tTopologyCopy = std::make_shared<Hexahedron_8_Topology<Real, Integer, Real_Matrix, Integer_Matrix>>(mNodeIndices);
        return tTopologyCopy;
    }
private:
    moris::Matrix<Integer, Integer_Matrix> mNodeIndices;
    Hexahedron_8_Basis_Function<Real,Real_Matrix> mBasisFunction;


};
}



#endif /* SRC_TOPOLOGY_CL_XTK_HEXAHEDRON_8_TOPOLOGY_HPP_ */
