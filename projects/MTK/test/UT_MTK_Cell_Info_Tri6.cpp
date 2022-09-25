/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MTK_Cell_Info_Tri6.cpp
 *
 */

#include "catch.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Tri6.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
namespace moris
{

    TEST_CASE( "MTK Tri6 Cell Info", "[TRI_6]" )
    {
        moris::real tEps = 1e-8;

        // construct a Hex27 info
        mtk::Cell_Info_Tri6 tTri6ConnInfo;

        CHECK( tTri6ConnInfo.get_num_verts() == 6 );

        // matrix for shape functions
        Matrix< DDRMat > tN;

        for ( moris::uint i = 0; i < tTri6ConnInfo.get_num_verts(); i++ )
        {
            Matrix< DDRMat > tCoords = tTri6ConnInfo.get_vertex_loc_coord( (moris_index)i );

            tTri6ConnInfo.eval_N( tCoords, tN );

            CHECK( std::abs( tN( i ) ) - 1.0 < tEps );
            CHECK( moris::norm( tN ) - 1.0 < tEps );
        }
    }
}    // namespace moris
