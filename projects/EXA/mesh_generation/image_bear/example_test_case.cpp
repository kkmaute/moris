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

TEST_CASE( "image_bear",
        "[moris],[example],[mesh_generation],[image_bear]" )
{
    // get the root path for MORIS
    std::string tMorisRoot = moris::get_base_moris_dir();

    // get the relevant file paths
    std::string tExamplesDir = tMorisRoot + "/share/doc/mesh_generation/examples/";
    std::string tXmlSource   = tExamplesDir + "Bear_Example.xml";
    std::string tHdf5Source  = tExamplesDir + "Bear.hdf5";
    std::string tTargetDir   = std::getenv( "PWD" );

    // copy relevant files from example directory into working directory
    std::string tMoveStringXml = "cp -f " + tXmlSource + " " + tTargetDir;
    std::string tMoveStringHdf5 = "cp -f " + tHdf5Source + " " + tTargetDir;
    
    MORIS_ERROR( std::system( tMoveStringXml.c_str() ) > 0, "failure" );
    MORIS_ERROR( std::system( tMoveStringHdf5.c_str() )  > 0, "failure" );

    // FIXME: remove the previous and use the below to use std::filesystem paths/copies, requires C++17 with gcc version >=8
    // // get all the relevant file paths
    // std::filesystem::path tExamplesDir = tMorisRoot + "/share/doc/mesh_generation/examples/";
    // std::filesystem::path tXmlSource   = tExamplesDir + "Bear_Example.xml";
    // std::filesystem::path tHdf5Source  = tExamplesDir + "Bear.hdf5";
    // std::filesystem::path tTargetDir   = std::getenv( "PWD" );
    // std::filesystem::path tXmlTarget   = tTargetDir / tXmlSource.filename();
    // std::filesystem::path tHdf5Target  = tTargetDir / tHdf5Source.filename();

    // // copy relevant files from example directory into working directory
    // std::filesystem::copy_file( tXmlSource, tXmlTarget, fs::copy_options::overwrite_existing );
    // std::filesystem::copy_file( tHdf5Source, tHdf5Target, fs::copy_options::overwrite_existing );

    // define command line call
    int argc = 3;

    char tString1[] = "";
    char tString2[] = "--meshgen";
    char tString3[] = "Bear_Example.xml";

    char* argv[ 3 ] = { tString1, tString2, tString3 };

    if ( par_size() == 1 )
    {
        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // TODO: perform sanity checks on output
        // check_results( ... );
    }
}
