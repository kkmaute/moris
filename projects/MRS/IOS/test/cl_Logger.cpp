/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Logger.cpp
 *
 */

#include <catch.hpp>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>
#include "cl_Communication_Tools.hpp"

// MORIS header files.
#include "ios.hpp"

#include "cl_Logger.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "MORIS_LOG_TEST",
        "[moris],[ios],[cl_Logger],[MORIS_LOG],[log_test]" )
{
    MORIS_SECTION( "                           Test Section Number  %-5i  \n", 1 );

    MORIS_SECTION( " Test check severity level warning " );

    MORIS_LOG_WARNING( "Test check severity level warning  || %-5i   ", 44 );
}

// Parsers need iteration records even when expensive timing output is disabled.
// Each MPI size/rank uses its own file, including concurrent CTest size variants.
TEST_CASE( "Logger iteration metadata is independent of timing", "[moris][ios][logger_metadata]" )
{
    const std::string tBase = "IOS_iteration_metadata_" + std::to_string( moris::par_size() ) + ".log";
    const std::string tFile = tBase + "." + std::to_string( moris::par_rank() );
    {
        moris::Logger tLogger;

        // Logger emits iteration records on its configured output rank.
        tLogger.mOutputRank = moris::par_rank();
        tLogger.initialize( tBase, 3, 2 );
        tLogger.sign_in( "Regression", "Metadata", "Iteration" );
        tLogger.iterate();
        tLogger.sign_out();
    }

    std::ifstream tInput( tFile );
    REQUIRE( tInput.good() );
    const std::string tText( ( std::istreambuf_iterator< char >( tInput ) ), std::istreambuf_iterator< char >() );
    CHECK( tText.find( ";Iteration;1" ) != std::string::npos );
    tInput.close();
    std::remove( tFile.c_str() );
}
