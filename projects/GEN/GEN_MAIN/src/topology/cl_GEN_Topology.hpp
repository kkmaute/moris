/*
 * cl_XTK_Topology.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_TOPOLOGY_HPP_
#define PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_TOPOLOGY_HPP_


#include "linalg_typedefs.hpp"
#include "../additional/cl_GEN_Basis_Function.hpp"
#include "../additional/cl_GEN_Enums.hpp" // For topology type enums

namespace moris
{
namespace ge
{

    
    enum class GEN_Topology_Type
{
    EDGE, // Edge with 2 Node
    TRI_3,
    TRI_6,
    QUAD_4,
    TET_4,
    TET_10,
    HEXA_8, // hexahedron with 8 nodes topology

};
    
class GEN_Topology
{
public:
    virtual
    ~GEN_Topology()
    {

    }


    virtual enum GEN_Topology_Type get_topology_type() const = 0;

    virtual moris::Matrix< moris::IndexMat > const & get_node_indices() const = 0;

    virtual Basis_Function const & get_basis_function() const = 0;

    virtual std::shared_ptr<GEN_Topology> copy() const = 0;

};
}
}


#endif /* PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_TOPOLOGY_HPP_ */
