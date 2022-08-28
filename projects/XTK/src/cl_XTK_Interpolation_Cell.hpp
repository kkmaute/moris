/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Interpolation_Cell.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_CELL_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_CELL_HPP_

#include "cl_MTK_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
using namespace moris;

namespace moris
{
    namespace mtk
    {
    class Vertex;
    }
}

namespace xtk
{

class Interpolation_Cell: public mtk::Cell
{
public:
    Interpolation_Cell(moris_id                  aCellId,
                       moris_index               aCellIndex,
                       moris_id                  aCellOwner,
                       std::shared_ptr<moris::mtk::Cell_Info> aConnectivity)
                       :Cell(aCellId,aCellIndex,aCellOwner,aConnectivity)
    {
    }

    Interpolation_Cell(){};

    virtual moris::Cell< mtk::Vertex* >
    get_vertex_pointers() const = 0;

    virtual Matrix< DDRMat >
    get_vertex_coords() const = 0;

};

}

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_CELL_HPP_ */

