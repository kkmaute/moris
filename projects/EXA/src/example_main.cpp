/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_main.cpp
 *
 */

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include <filesystem>
#include <iostream>
#include <regex>

#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Logger.hpp"                   // MRS/IOS/src
#include "banner.hpp"                      // COR/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

using namespace moris;
namespace fs = std::filesystem;

//---------------------------------------------------------------

int main( int argc, char *argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // set severity level 0 - all outputs
    gLogger.initialize( 2 );

    // print banner
    moris::print_banner( argc, argv );

    // remove files from previous test runs
    if ( par_rank() == 0 )
    {
        std::string path = "./";
        std::regex  pattern1( ".*exo.*" );
        std::regex  pattern2( ".*hdf5.*" );

        for ( const auto &entry : fs::directory_iterator( path ) )
        {
            const std::string filename = entry.path().filename().string();

            if ( fs::is_regular_file( entry ) && ( std::regex_match( filename, pattern1 ) || std::regex_match( filename, pattern2 ) ) )
            {
                fs::remove( entry.path() );
            }
        }
    }

    // Run Tests
    int tRet = Catch::Session().run( argc, argv );

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return tRet;
}
