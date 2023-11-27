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
}
