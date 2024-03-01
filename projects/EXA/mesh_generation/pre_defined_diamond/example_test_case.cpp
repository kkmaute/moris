/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 */

#include <catch.hpp>
// #include <filesystem> // FIXME: use filesystem with C++17 enabled

#include "cl_Logger.hpp"                  // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src
#include "cl_Communication_Tools.hpp"     // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

#include "HDF5_Tools.hpp"
#include "paths.hpp"

using namespace moris;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

TEST_CASE( "pre_defined_diamond",
        "[moris],[example],[mesh_generation],[pre_defined_diamond]" )
{
    // get the root path for MORIS
    std::string tMorisRoot = moris::get_base_moris_dir();

    // FIXME: use std::filesystem paths/copies instead, but requires C++17 with gcc version >=8
    // get the relevant file paths
    std::string tExamplesDir = tMorisRoot + "/share/doc/mesh_generation/examples/";
    std::string tXmlSource   = tExamplesDir + "Rotated_Square_Example.xml";
    std::string tTargetDir   = std::getenv( "PWD" );

    // copy relevant files from example directory into working directory
    std::string tMoveStringXml = "cp -f " + tXmlSource + " " + tTargetDir;
    
   // fixme: should be changed to C++ commands rather than system calls
    MORIS_ERROR( std::system( tMoveStringXml.c_str() ) > 0, "failure" );

    // define command line call
    int argc = 3;

    char tString1[] = "";
    char tString2[] = "--meshgen";
    char tString3[] = "Rotated_Square_Example.xml";

    char* argv[ 3 ] = { tString1, tString2, tString3 };

    if ( par_size() == 1 )
    {
        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // TODO: perform sanity checks on output
        // check_results( ... );
    }
}
