/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Derived_Node.hpp
 *
 */

#pragma once

#include "cl_GEN_Node.hpp"

namespace moris::ge
{
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

    class Derived_Node : public Node
    {
      private:
        Cell< Locator > mLocators;
        Matrix< DDRMat > mCoordinates;

      public:
        /**
         * Derived node constructor, using other nodal information with locators.
         *
         * @param aIndex Node index
         * @param aBaseNodes Base nodes
         */
        Derived_Node(
                uint            aIndex,
                Cell< Locator > aLocators );

        /**
         * Gets the coordinates of this node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_coordinates() override;

        /**
         * Gets the locators of this derived node
         *
         * @return Locators
         */
        const Cell< Locator >& get_locators();
    };
}
