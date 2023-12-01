/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Parent_Base_Node.hpp
 *
 */

#pragma once

#include "cl_GEN_Parent_Node.hpp"
#include "cl_Matrix.hpp"

namespace moris::ge
{
    // Forward declare base node FIXME make base node once a base node is specifically requested from the node manager
    class Base_Node;

    class Parent_Base_Node : public Parent_Node
    {
      private:
        Base_Node*       mBaseNode;
        Matrix< DDRMat > mParametricCoordinates;

      public:
        /**
         * Constructor for a parent node given a base node and parametric coordinates
         *
         * @param aBaseNode Base node
         * @param aParametricCoordinates Parametric coordinates
         */
        Parent_Base_Node(
                Base_Node*              aBaseNode,
                const Matrix< DDRMat >& aParametricCoordinates );

        /**
         * Gets the index of the underlying node
         *
         * @return Node index
         */
        uint get_index() const override;

        /**
         * Gets the coordinates of the underlying node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_global_coordinates() const override;

        /**
         * Gets the parametric coordinates of this derived node relative to its locators
         *
         * @return Parametric coordinates
         */
        const Matrix< DDRMat >& get_parametric_coordinates() const override;

      private:
        /**
         * Gets the underlying node, specifically used to create a basis node from a parent node.
         *
         * @return Underlying node object
         */
        Node* get_node() const override;
    };
}
