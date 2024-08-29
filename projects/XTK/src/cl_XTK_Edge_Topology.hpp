/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Edge_Topology.hpp
 *
 */

#ifndef SRC_TOPOLOGY_CL_XTK_EDGE_TOPOLOGY_HPP_
#define SRC_TOPOLOGY_CL_XTK_EDGE_TOPOLOGY_HPP_

#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Basis_Function.hpp"

// Basis Functions
#include "cl_XTK_Linear_Basis_Functions.hpp"

namespace moris::xtk
{
    class Edge_Topology : public Topology
    {

      public:
        Edge_Topology()
        {
        }

        Edge_Topology( Matrix< IndexMat > const &aNodeIndices )

        {
            this->set_node_indices( aNodeIndices );
        }

        // Required Interface Functions
        enum Topology_Type get_topology_type() const override
        {
            return Topology_Type::EDGE;
        }
        Matrix< IndexMat > const &get_node_indices() const override
        {
            return mNodeIndices;
        }

        Basis_Function const &get_basis_function() const override
        {
            return mBasisFunction;
        }

        void set_node_indices( Matrix< IndexMat > const &aNodeIndices )
        {
            MORIS_ASSERT( aNodeIndices.n_cols() == 2, "Should be 2 associated with a edge topology" );

            mNodeIndices = aNodeIndices.copy();
        }

        std::shared_ptr< Topology > copy() const override
        {
            std::shared_ptr< Topology > tTopologyCopy;
            tTopologyCopy = std::make_shared< Edge_Topology >( mNodeIndices );
            return tTopologyCopy;
        }

      private:
        Matrix< IndexMat >    mNodeIndices;
        Linear_Basis_Function mBasisFunction;
    };
}    // namespace moris::xtk

#endif /* SRC_TOPOLOGY_CL_XTK_EDGE_TOPOLOGY_HPP_ */
