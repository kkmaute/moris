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

#include "cl_GEN_Node.hpp"

namespace moris::ge
{
    class Base_Node : public Node
    {
      private:
        Matrix< DDRMat > mCoordinates;

      public:
        /**
         * Base node constructor. Requires a node index as well as coordinates.
         *
         * @param aIndex Node index
         * @param aCoordinates Node coordinates
         */
        Base_Node(
                uint                    aIndex,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * Gets the coordinates of this node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_coordinates() override;
    };
}
