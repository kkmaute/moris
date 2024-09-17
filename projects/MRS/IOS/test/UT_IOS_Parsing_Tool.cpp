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
#include "cl_Vector.hpp"

// ----------------------------------------------------------------------------

namespace moris
{
    TEST_CASE( "String to vector of matrices", "[moris],[ios],[cl_Logger],[String_to_vector_of_matrices]")
    {
        Vector< Matrix< DDRMat > > tCellMat;

        std::string tString = "1.0, 2.0; 3.0, 11.0 / 1.0 , 4.0 , 3.0 ; 2.0 , 5.0 ,8.0 ";

        string_to_vector_of_matrices( tString, tCellMat );

        CHECK( tCellMat(0)(0,0) == 1.0  );        CHECK( tCellMat(1)(0,0) == 1.0  );
        CHECK( tCellMat(0)(0,1) == 2.0  );        CHECK( tCellMat(1)(0,2) == 3.0 );
        CHECK( tCellMat(0)(1,0) == 3.0  );        CHECK( tCellMat(1)(1,0) == 2.0  );
        CHECK( tCellMat(0)(1,1) == 11.0  );       CHECK( tCellMat(1)(1,1) == 5.0 );
    }
}

