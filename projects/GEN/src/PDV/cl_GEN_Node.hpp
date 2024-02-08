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

namespace moris::gen
{
    // Forward declare basis node
    class Basis_Node;

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
        uint get_index() const;

        /**
         * Gets the coordinates of this node
         *
         * @return Node coordinates
         */
        virtual const Matrix< DDRMat >& get_global_coordinates() const = 0;

        /**
         * Get the value of a coordinate of this node
         *
         * @param aCoordinateIndex index of the coordinate, obtained from casting the related PDV coordinate type
         * @return Coordinate value
         */
        virtual real get_coordinate_value( uint aCoordinateIndex ) const = 0;

        /**
         * Gets the locator nodes of this node.
         * Locator nodes are the most derived basis nodes that can determine the location of this node.
         *
         * @return Locator nodes
         */
        virtual const Vector< Basis_Node >& get_locator_nodes() const = 0;

        /**
         * Gets if this node's position depends on ADVs.
         *
         * @return ADV dependence (default false)
         */
        virtual bool depends_on_advs() const;

        /**
         * Appends the sensitivities of this node's global coordinates with respect to ADVs. By default, does nothing.
         *
         * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
         * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
         */
        virtual void append_dcoordinate_dadv(
                Matrix< DDRMat >&       aCoordinateSensitivities,
                const Matrix< DDRMat >& aSensitivityFactor ) const;

        /**
         * Gets the ADV IDs that determine the coordinates of this node. By default, returns an empty vector.
         *
         * @return ADV ID vector
         */
        virtual Matrix< DDSMat > get_coordinate_determining_adv_ids() const;
    };
}
