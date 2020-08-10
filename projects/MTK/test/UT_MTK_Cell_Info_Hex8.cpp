/*
 * ut_MTK_Cell_Info_Hex27.cpp
 *
 *  Created on: Aug 7, 2020
 *      Author: kedo3694
 */
#include "catch.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

namespace moris
{

TEST_CASE("MTK Hex8 Cell Info","[HEX_8]")
{
    moris::real tEps = 1e-8;

    // construct a Hex27 info
    mtk::Cell_Info_Hex8 tHex8ConnInfo;

    CHECK(tHex8ConnInfo.get_num_verts() == 8);

    // matrix for shape functions
    Matrix<DDRMat> tN(1,8);

    // check the local coords of verts using eval_N
    for(moris::uint i = 0; i < tHex8ConnInfo.get_num_verts(); i++)
    {

        Matrix<DDRMat> tCoords = tHex8ConnInfo.get_vertex_loc_coord((moris_index)i);

        tHex8ConnInfo.eval_N(tCoords,tN);

        CHECK(std::abs(tN(i)) - 1.0 < tEps);
        CHECK(moris::norm(tN) - 1.0 < tEps);

    }

}


}
