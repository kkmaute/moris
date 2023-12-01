/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Parent_Derived_Node.hpp
 *
 */

#pragma once

#include "cl_GEN_Parent_Node.hpp"
#include "cl_Matrix.hpp"

namespace moris::ge
{
    // Forward declare derived node
    class Derived_Node;

    class Parent_Derived_Node : public Parent_Node
    {
      private:
        Derived_Node* mDerivedNode;

      public:
        /**
         * Constructor for a parent node given a derived node
         *
         * @param aDerivedNode Derived node
         */
        Parent_Derived_Node( Derived_Node* aDerivedNode );

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
