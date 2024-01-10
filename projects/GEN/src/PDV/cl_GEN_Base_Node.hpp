/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Base_Node.hpp
 *
 */

#pragma once

#include "cl_GEN_Node.hpp"

namespace moris::ge
{
    class Base_Node : public Node
    {
      private:
        Matrix< DDRMat > mCoordinates;
        static Cell< Basis_Node > mDummyLocatorNodes;

      public:
        /**
         * Base node constructor. Requires a node index as well as coordinates.
         *
         * @param aIndex Node index
         * @param aCoordinates Node coordinates
         */
        Base_Node(
                uint                    aIndex,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * Gets the coordinates of this node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_global_coordinates() override;

        /**
         * Gets the locator nodes of this node.
         * Locator nodes are the most derived basis nodes that can determine the location of this node.
         * Background nodes do not have locator nodes.
         *
         * @return Locator nodes (empty)
         */
        const Cell< Basis_Node >& get_locator_nodes() override;
    };
}
