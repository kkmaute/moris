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
}

namespace moris::ge
{
    // Forward declare basis node
    class Basis_Node;

    class Derived_Node : public Node
    {
      private:
        Cell< Basis_Node > mBasisNodes;
        Matrix< DDRMat > mGlobalCoordinates;
        Matrix< DDRMat > mParametricCoordinates;

      public:
        /**
         * Derived node constructor, using other nodal information with locators.
         *
         * @param aIndex Node index
         * @param aBaseNodes Base nodes
         */
        Derived_Node(
                uint                    aIndex,
                const Cell< Node* >&    aBaseNodes,
                const Matrix< DDRMat >& aParametricCoordinates,
                mtk::Geometry_Type      aGeometryType );

        /**
         * Gets the global coordinates of this node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_global_coordinates() override;

        /**
         * Gets the parametric coordinates of this derived node relative to its locators
         *
         * @return Parametric coordinates
         */
        const Matrix< DDRMat >& get_parametric_coordinates();

        /**
         * Gets the locators of this derived node
         *
         * @return Locators
         */
        const Cell< Basis_Node >& get_basis_nodes();
    };
}
