/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Info_Factory.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_Cell_Info_Factory_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_Cell_Info_Factory_HPP_

#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
namespace mtk
{
class Cell_Info;
class Cell_Info_Factory
{
public:
    moris::mtk::Cell_Info*
    create_cell_info(enum CellTopology aCellTopo);

    std::shared_ptr<moris::mtk::Cell_Info>
    create_cell_info_sp(enum CellTopology aCellTopo);

    moris::mtk::Cell_Info*
    create_cell_info(enum Geometry_Type       aCellGeom,
                     enum Interpolation_Order aInterpOrder);

    std::shared_ptr<moris::mtk::Cell_Info>
    create_cell_info_sp(enum Geometry_Type       aCellGeom,
                        enum Interpolation_Order aInterpOrder);

};
}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_Cell_Info_Factory_HPP_ */

