/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MTK_Cell_Info_Hex64.cpp
 *
 */

#include "catch.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Hex64.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

namespace moris
{

TEST_CASE("MTK Hex64 Cell Info","[HEX_64]")
{
    moris::real tEpsilon = 1e-8;

    // construct a cell info
    mtk::Cell_Info_Hex64 tHex64ConnInfo;

    CHECK(tHex64ConnInfo.get_num_verts() == 64);

    // matrix for shape functions
    Matrix<DDRMat> tN;

    for(moris::uint i = 0; i < tHex64ConnInfo.get_num_verts(); i++)
    {

        Matrix<DDRMat> tCoords = tHex64ConnInfo.get_vertex_loc_coord((moris_index)i);

        tHex64ConnInfo.eval_N(tCoords,tN);

        CHECK(std::abs(tN(i)) - 1.0 < tEpsilon);
        CHECK(moris::norm(tN) - 1.0 < tEpsilon);

    }

}

}

