/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Node.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"

namespace moris::ge
{
    class Node
    {
      private:
        uint mIndex;

      public:
        /**
         * Constructor with a given node index
         *
         * @param aIndex Index of this node
         */
        explicit Node( uint aIndex );

        /**
         * Default destructor
         */
        virtual ~Node() = default;

        /**
         * Gets the index of this node
         *
         * @return Stored node index
         */
        uint get_index();

        /**
         * Gets the coordinates of this node
         *
         * @return Node coordinates
         */
        virtual const Matrix< DDRMat >& get_coordinates() = 0;
    };

    class Locator
    {
      private:
        Node* mNode;
        real mBasis;

      public:
        /**
         * Constructor with a generic node and basis value.
         *
         * @param aNode GEN node
         * @param aBasis Basis value for locating a derived node
         */
        Locator(
                Node* aNode,
                real aBasis );

        /**
         * Gets the index of the underlying node.
         *
         * @return Node index
         */
        uint get_index();

        /**
         * Gets the coordinates of the underlying node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_coordinates();

        /**
         * Gets the basis of this locator.
         *
         * @return Basis value
         */
        real get_basis();
    };
}
