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

#if __has_include( <filesystem>)
#include <filesystem>
namespace filesystem = std::filesystem;
#else
#include <experimental/filesystem>
namespace filesystem = std::experimental::filesystem;
#endif

using namespace moris;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

TEST_CASE( "obj_sphere",
        "[moris],[example],[mesh_generation],[obj_sphere]" )
{
    // get the root path for MORIS
    std::string tMorisRoot = moris::get_base_moris_dir();

    // get all the relevant file paths
    filesystem::path tMorisRootPath( tMorisRoot );
    filesystem::path tExamplesDir = tMorisRootPath / "projects/EXA/mesh_generation/obj_sphere/";
    filesystem::path tXmlSource   = tExamplesDir / "OBJ_Sphere.xml";
    filesystem::path tHdf5Source  = tExamplesDir / "Sphere.obj";
    filesystem::path tTargetDir   = filesystem::current_path();
    filesystem::path tXmlTarget   = tTargetDir / tXmlSource.filename();
    filesystem::path tHdf5Target  = tTargetDir / tHdf5Source.filename();

    // copy relevant files from example directory into working directory
    filesystem::copy_file( tXmlSource, tXmlTarget, filesystem::copy_options::overwrite_existing );
    filesystem::copy_file( tHdf5Source, tHdf5Target, filesystem::copy_options::overwrite_existing );

    // define command line call
    int argc = 3;

    char tString1[] = "";
    char tString2[] = "--meshgen";
    char tString3[] = "OBJ_Sphere.xml";

    char* argv[ 3 ] = { tString1, tString2, tString3 };

    if ( par_size() == 1 )
    {
        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // TODO: perform sanity checks on output
        // check_results( ... );
    }
}
