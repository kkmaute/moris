/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Parent_Node.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"
#include "cl_GEN_Basis_Node.hpp"

namespace moris::ge
{
    class Parent_Node
    {
        // Basis node needs to be able to be constructed using the node object from this class
        friend Basis_Node::Basis_Node( const Parent_Node& aParentNode, real aBasis );

      public:
        /**
         * Default constructor
         */
        Parent_Node() = default;

        /**
         * Default destructor
         */
        virtual ~Parent_Node() = default;

        /**
         * Gets the index of the underlying node
         *
         * @return Node index
         */
        virtual uint get_index() const = 0;

        /**
         * Gets the coordinates of the underlying node
         *
         * @return Node coordinates
         */
        virtual const Matrix< DDRMat >& get_global_coordinates() const = 0;

        /**
         * Gets the parametric coordinates of this derived node relative to its locators
         *
         * @return Parametric coordinates
         */
        virtual const Matrix< DDRMat >& get_parametric_coordinates() const = 0;

      private:
        /**
         * Gets the underlying node, specifically used to create a basis node from a parent node.
         *
         * @return Underlying node object
         */
        virtual Node* get_node() const = 0;
    };
}
