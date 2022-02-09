/*
 * cl_XTK_Topology.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_TOPOLOGY_HPP_
#define SRC_XTK_CL_XTK_TOPOLOGY_HPP_

#include "cl_XTK_Enums.hpp"    // For topology type enums
#include "cl_XTK_Basis_Function.hpp"

#include "linalg_typedefs.hpp"

namespace xtk
{
    enum class Topology_Type
    {
        EDGE,    // Edge with 2 Node
        TRI_3,
        TRI_6,
        QUAD_4,
        TET_4,
        TET_10,
        HEXA_8,    // hexahedron with 8 nodes topology

    };

    class Topology
    {
      public:
        virtual ~Topology()
        {
        }

        virtual enum Topology_Type get_topology_type() const = 0;

        virtual moris::Matrix< moris::IndexMat > const & get_node_indices() const = 0;

        virtual Basis_Function const & get_basis_function() const = 0;

        virtual std::shared_ptr< Topology > copy() const = 0;
    };
}    // namespace xtk
#endif /* SRC_XTK_CL_XTK_TOPOLOGY_HPP_ */
