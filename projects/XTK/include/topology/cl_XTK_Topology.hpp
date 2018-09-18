/*
 * cl_XTK_Topology.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_TOPOLOGY_HPP_
#define SRC_XTK_CL_XTK_TOPOLOGY_HPP_


#include"xtk/cl_XTK_Enums.hpp" // For topology type enums

#include "topology/cl_XTK_Basis_Function.hpp"

namespace xtk
{

    
    enum class Topology_Type
{
    EDGE, // Edge with 2 Node
    TRI_3,
    QUAD_4,
    TET_4,
    TET_10,
    HEXA_8, // hexahedron with 8 nodes topology

};
    
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Topology
{
public:
    ~Topology()
    {

    }


    virtual enum Topology_Type get_topology_type() const = 0;

    virtual moris::Matrix< Integer_Matrix > const & get_node_indices() const = 0;

    virtual Basis_Function<Real,Real_Matrix> const & get_basis_function() const = 0;

    virtual std::shared_ptr<Topology<Real, Integer, Real_Matrix, Integer_Matrix>> copy() const = 0;

};
}



#endif /* SRC_XTK_CL_XTK_TOPOLOGY_HPP_ */
