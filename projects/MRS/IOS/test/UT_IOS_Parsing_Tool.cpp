/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_IOS_Parsing_Tool.cpp
 *
 */

#include <catch.hpp>

// MORIS header files.
#include "fn_Parsing_Tools.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

// ----------------------------------------------------------------------------

namespace moris
{
    TEST_CASE( "String_To_Cell_Mat", "[moris],[ios],[cl_Logger],[String_to_Cell_Mat]")
    {
        Cell< Matrix< DDUMat > > tCellMat;

        std::string tString = "1, 2; 3, 11, 5; 7; 8, 10";

        string_to_cell_mat( tString, tCellMat );

        CHECK( tCellMat(0)(0,0) == 1  );        CHECK( tCellMat(0)(1,0) == 2  );
        CHECK( tCellMat(1)(0,0) == 3  );        CHECK( tCellMat(1)(1,0) == 11 );
        CHECK( tCellMat(1)(2,0) == 5  );        CHECK( tCellMat(2)(0,0) == 7  );
        CHECK( tCellMat(3)(0,0) == 8  );        CHECK( tCellMat(3)(1,0) == 10 );
    }

    TEST_CASE( "String_To_Cell_Mat_2", "[moris],[ios],[cl_Logger],[String_to_Cell_Mat_2]")
    {
        Cell< Matrix< DDRMat > > tCellMat;

        std::string tString = "1.0, 2.0; 3.0, 11.0 / 1.0 , 4.0 , 3.0 ; 2.0 , 5.0 ,8.0 ";

        string_to_cell_mat_2( tString, tCellMat );

        CHECK( tCellMat(0)(0,0) == 1.0  );        CHECK( tCellMat(1)(0,0) == 1.0  );
        CHECK( tCellMat(0)(0,1) == 2.0  );        CHECK( tCellMat(1)(0,2) == 3.0 );
        CHECK( tCellMat(0)(1,0) == 3.0  );        CHECK( tCellMat(1)(1,0) == 2.0  );
        CHECK( tCellMat(0)(1,1) == 11.0  );       CHECK( tCellMat(1)(1,1) == 5.0 );
    }
}

