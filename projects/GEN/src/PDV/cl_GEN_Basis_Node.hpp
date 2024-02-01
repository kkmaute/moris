/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Basis_Node.hpp
 *
 */

#pragma once

#include "cl_GEN_Node.hpp"

namespace moris::gen
{
    // Forward declare parent node
    class Parent_Node;

    class Basis_Node
    {
      private:
        const Node& mNode;
        real mBasis;

      public:
        /**
         * Constructor with a generic node and basis value.
         *
         * @param aNode GEN node
         * @param aBasis Basis value for locating a derived node
         */
        Basis_Node(
                const Node& aNode,
                real        aBasis );

        /**
         * Constructor with a parent node and basis value.
         *
         * @param aParentNode Parent node
         * @param aBasis Basis value for locating a derived node
         */
        Basis_Node(
                const Parent_Node& aParentNode,
                real               aBasis );

        /**
         * Gets the index of the underlying node.
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
         * Gets the basis of this locator.
         *
         * @return Basis value
         */
        real get_basis() const;

        /**
         * Gets the locator nodes of the underlying node.
         *
         * @return Locator nodes
         */
        const Cell< Basis_Node >& get_locator_nodes() const;

        /**
         * Gets if the underlying node's position depends on ADVs.
         *
         * @return ADV dependence
         */
        bool depends_on_advs() const;

        /**
         * Appends the sensitivities of the underlying node's global coordinates with respect to ADVs.
         *
         * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
         * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
         */
        void append_dcoordinate_dadv(
                Matrix< DDRMat >&       aCoordinateSensitivities,
                const Matrix< DDRMat >& aSensitivityFactor ) const;

        /**
         * Gets the coordinate determining ADV IDs for the underlying node.
         *
         * @return ADV ID vector
         */
        Matrix< DDSMat > get_coordinate_determining_adv_ids() const;
    };
}
