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

#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"
#include "cl_Matrix.hpp"

namespace moris::ge
{
    class Parent_Node
    {
        // Basis node needs to be able to be constructed using the node object from this class
        friend Basis_Node::Basis_Node( const Parent_Node& aParentNode, real aBasis );

      private:
        Node*            mNode;
        Matrix< DDRMat > mParametricCoordinates;

      public:
        /**
         * Constructor for a parent node given a base node and parametric coordinates
         *
         * @param aNode Base node
         * @param aParametricCoordinates Parametric coordinates
         */
        Parent_Node(
                Node*                   aNode,
                const Matrix< DDRMat >& aParametricCoordinates );

        /**
         * Gets the index of the underlying node
         *
         * @return Node index
         */
        uint get_index() const;

        /**
         * Gets the coordinates of the underlying node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_global_coordinates() const;

        /**
         * Gets the parametric coordinates of this derived node relative to its locators
         *
         * @return Parametric coordinates
         */
        const Matrix< DDRMat >& get_parametric_coordinates() const;
    };
}
