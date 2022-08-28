/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MTK_Cell_Info_Quad16.cpp
 *
 */

#include "catch.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Quad16.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

namespace moris
{

TEST_CASE("MTK Quad16 Cell Info","[Quad_16]")
{

    moris::real tEpsilon = 1e-8;

    // construct a cell info
    mtk::Cell_Info_Quad16 tCellInfo;

    CHECK(tCellInfo.get_num_verts() == 16);

    // matrix for shape functions
    Matrix<DDRMat> tN;

    for(moris::uint i = 0; i < tCellInfo.get_num_verts(); i++)
    {

        Matrix<DDRMat> tCoords = tCellInfo.get_vertex_loc_coord((moris_index)i);

        tCellInfo.eval_N(tCoords,tN);

        CHECK(std::abs(tN(i)) - 1.0 < tEpsilon);
        CHECK(moris::norm(tN) - 1.0 < tEpsilon);

    }

}

}

