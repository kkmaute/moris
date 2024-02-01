/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Background_Node.hpp
 *
 */

#pragma once

#include "cl_GEN_Node.hpp"

namespace moris::ge
{
    class Background_Node : public Node
    {
      private:
        Matrix< DDRMat > mCoordinates;
        static Cell< Basis_Node > mDummyLocatorNodes;

      public:
        /**
         * Background node constructor. Requires a node index as well as coordinates.
         *
         * @param aIndex Node index
         * @param aCoordinates Node coordinates
         */
        Background_Node(
                uint                    aIndex,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * Gets the coordinates of this node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_global_coordinates() const override;

        /**
         * Get the value of a coordinate of this node
         *
         * @param aCoordinateIndex index of the coordinate, obtained from casting the related PDV coordinate type
         * @return Coordinate value
         */
        real get_coordinate_value( uint aCoordinateIndex ) const override;

        /**
         * Gets the locator nodes of this node.
         * Locator nodes are the most derived basis nodes that can determine the location of this node.
         * Background nodes do not have locator nodes.
         *
         * @return Locator nodes (empty)
         */
        const Cell< Basis_Node >& get_locator_nodes() const override;
    };
}
