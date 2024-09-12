/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Edge.hpp
 *
 */

#ifndef PROJECTS_MTK_CL_MTK_EDGE_HPP_
#define PROJECTS_MTK_CL_MTK_EDGE_HPP_

#include "cl_MTK_Cell.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------
    class Edge : public Cell
    {
        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        Edge(){};

        //------------------------------------------------------------------------------

        ~Edge() override{};

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* PROJECTS_MTK_CL_MTK_EDGE_HPP_ */
