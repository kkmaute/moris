/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_IOS_File_To_Array.cpp
 *
 */

#include <catch.hpp>

// MORIS header files.
#include "cl_Matrix.hpp"
#include "cl_Ascii.hpp"
#include "paths.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_assert.hpp"

// ----------------------------------------------------------------------------

namespace moris
{

    TEST_CASE( "file to array", "[moris],[ios],[ascii],[file_to_array]")
    {
        // Input file
        std::string tPrefix = moris::get_base_moris_dir();
        std::string tFileName = tPrefix + "/projects/MRS/IOS/test/input_file.txt";

        // create ascii class
        Ascii tAscii= Ascii( tFileName, FileMode::OPEN_RDONLY);

        // matrix used to compare against
        Matrix< DDUMat > tRefMatrix = {{11,12,13,14},
                                       {21,22,23,24},
                                       {31,32,33,34}};

        uint tNumRows = tAscii.length();
        Matrix< DDUMat > tRow;
        string_to_mat( tAscii.line( 0 ), tRow );
        uint tNumCols = tRow.n_cols();

        Matrix< DDUMat > tMatrix( tNumRows,tNumCols );

        for( uint i = 0; i < tNumRows; i++ )
        {

            // storing ascii row to matrix
            string_to_mat( tAscii.line( i ), tRow );

            // are the number of columns the same?
            MORIS_ASSERT(tRow.n_cols() == tNumCols, "Inconsistent number of columns in input file");

            tMatrix({i,i},{0,tNumCols-1}) = tRow({0,0},{0,tNumCols-1});

            // checking stored values against reference matrix
            for ( uint j = 0; j < tNumCols; j++)
            {
                CHECK( tMatrix(i,j) == tRefMatrix(i,j) );
            }
        }
    }
}

