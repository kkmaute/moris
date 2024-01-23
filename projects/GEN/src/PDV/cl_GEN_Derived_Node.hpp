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

namespace moris::mtk
{
    enum class Geometry_Type;
    enum class Interpolation_Order;
}

namespace moris::ge
{
    // Forward declare basis node
    class Basis_Node;
    class Geometry;

    class Derived_Node : public Node
    {
      private:
        Cell< Basis_Node > mBackgroundNodes;
        Matrix< DDRMat > mGlobalCoordinates;
        Matrix< DDRMat > mParametricCoordinates;
        static inline bool gOverrideLinearInterpolation = false;

      public:
        /**
         * Derived node constructor, using other nodal information with background nodes.
         *
         * @param aIndex Node index
         * @param aBackgroundNodes Background nodes
         * @param aParametricCoordinates Parametric coordinates inside the background element
         * @param aGeometryType Geometry type of the background element
         * @param aInterpolationOrder Interpolation order of the background element
         */
        Derived_Node(
                uint                     aIndex,
                const Cell< Node* >&     aBackgroundNodes,
                const Matrix< DDRMat >&  aParametricCoordinates,
                mtk::Geometry_Type       aGeometryType,
                mtk::Interpolation_Order aInterpolationOrder );

        /**
         * Gets the global coordinates of this node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_global_coordinates() const override;

        /**
         * Gets the parametric coordinates of this derived node relative to its locators
         *
         * @return Parametric coordinates
         */
        const Matrix< DDRMat >& get_parametric_coordinates() const;

        /**
         * Gets the basis nodes of this derived node
         *
         * @return Basis nodes
         */
        const Cell< Basis_Node >& get_background_nodes() const;

        /**
         * Gets the locator nodes of this derived node.
         * Locator nodes are the most derived basis nodes that can determine the location of this node.
         * For derived nodes, these are the background nodes, and for intersection nodes these are its parents.
         *
         * @return Locator nodes
         */
        const Cell< Basis_Node >& get_locator_nodes() const override;

        /**
         * Gets if this derived node can be determined that it is on a specific interface without any field evaluation.
         *
         * @param aGeometry Potential interface geometry
         * @return If this node is on the requested interface
         */
        virtual bool is_on_interface( Geometry* aGeometry ) const;

        /**
         * Sets the flag for overriding linear interpolation, for when multilinear intersections are being used.
         */
        static void set_override_linear_interpolation();
    };
}
