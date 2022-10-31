/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * banner.hpp
 *
 */

#ifndef SRC_CORE_BANNER_HPP_
#define SRC_CORE_BANNER_HPP_

#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>

#include "cl_Communication_Tools.hpp"    // COM/src

namespace moris
{
    //------------------------------------------------------------------------------
    /*
     * prints a logo. Created using figlet.
     */
    inline void
    print_logo()
    {
        std::fprintf( stdout, "\n" );
        std::fprintf( stdout, "            .___  ___.   ______   .______       __       _______.\n" );
        std::fprintf( stdout, "            |   \\/   |  /  __  \\  |   _  \\     |  |     /       |\n" );
        std::fprintf( stdout, "            |  \\  /  | |  |  |  | |  |_)  |    |  |    |   (----`\n" );
        std::fprintf( stdout, "            |  |\\/|  | |  |  |  | |      /     |  |     \\   \\\n" );
        std::fprintf( stdout, "            |  |  |  | |  `--'  | |  |\\  \\----.|  | .----)   |\n" );
        std::fprintf( stdout, "            |__|  |__|  \\______/  | _| `._____||__| |_______/\n" );
        std::fprintf( stdout, "\n" );
        std::fprintf( stdout, "            Copyright (c) 2022 University of Colorado Boulder\n" );
        std::fprintf( stdout, "                               Aerospace Mechanics Research Center\n" );
        std::fprintf( stdout, "\n" );
        std::fprintf( stdout, "\n" );
    }

    //------------------------------------------------------------------------------
    /*
     * tries to get the cpu type, returns unknown if it fails
     */
    inline std::string
    get_cpu_info()
    {
        // test if proc/cpuinfo exists
        std::ifstream tProcCpuInfo( "/proc/cpuinfo" );

        // test if file exists
        if ( tProcCpuInfo )
        {
            // temporary string
            std::string tString;

            // size for buffer
            const int tBufferSize = 128;

            // temporary buffer
            char tBuffer[ tBufferSize ];

            // create pointer for stream object
            std::shared_ptr< FILE > tStream(
                    popen( "cat /proc/cpuinfo | grep \"model name\" | head -n 1 2>&1", "r" ),
                    pclose );

            if ( tStream )
            {
                // read result from command line
                while ( !feof( tStream.get() ) )
                {
                    if ( fgets( tBuffer, tBufferSize, tStream.get() ) != nullptr )
                    {
                        tString.append( tBuffer );
                    }
                }

                // trim string
                auto tStart = tString.find_first_of( ':' ) + 2;
                auto tEnd   = tString.find_last_not_of( '\n' );

                // return info
                return tString.substr( tStart, ( tEnd - tStart + 1 ) );
            }
            else
            {
                return "unknown";
            }
        }
        else
        {
            return "unknown";
        }
    }

    //------------------------------------------------------------------------------

    /*
     * Prints a welcome banner similar to FEMDOC.
     */
    inline void
    print_banner(
            int&  argc,
            char* argv[] )
    {
        // banner is only printed by first proc
        if ( par_rank() == 0 )
        {
            print_logo();

#if defined( MORIS_HAVE_DEBUG )
            std::fprintf( stdout, "     DEBUG flags are on.\n\n" );
#endif

            // Who?
            std::string tUser     = std::getenv( "USER" ) == nullptr ? "NA" : std::getenv( "USER" );
            std::string tHostName = std::getenv( "HOSTNAME" ) == nullptr ? "NA" : std::getenv( "HOSTNAME" );
            std::string tOsType   = std::getenv( "OSTYPE" ) == nullptr ? "NA" : std::getenv( "OSTYPE" );
            std::string tHostType = std::getenv( "HOSTTYPE" ) == nullptr ? "NA" : std::getenv( "HOSTTYPE" );

            std::fprintf( stdout, "     User/Host   : %s at %s ( %s / %s )\n",    //
                    tUser.c_str(),
                    tHostName.c_str(),
                    tOsType.c_str(),
                    tHostType.c_str() );

            std::string tCpuInfo = get_cpu_info();

            std::fprintf( stdout, "     CPU Info    : %s \n", tCpuInfo.c_str() );
            std::fprintf( stdout, "     Procs Used  : %i \n", (int)par_size() );

            // insert blank line
            std::fprintf( stdout, "\n" );

            // When built?
            std::fprintf( stdout, "     Build Date  : %s at %s\n", __DATE__, __TIME__ );

            // When run
            time_t tTimeStamp = time( NULL );
            std::fprintf( stdout, "     Date of Run : %s \n", ctime( &tTimeStamp ) );

            // What Matrix lib?
#ifdef MORIS_USE_ARMA
            std::fprintf( stdout, "     Matrix Lib  : Armadillo\n" );
#else
            std::fprintf( stdout, "     Matrix Lib  : Eigen\n" );
#endif

            // insert blank line
            std::fprintf( stdout, "\n" );

            // What?
            std::fprintf( stdout, "     Executable  : %s\n", argv[ 0 ] );

            std::fprintf( stdout, "     Arguments   :" );

            for ( int ia = 1; ia < argc; ++ia )
            {
                std::fprintf( stdout, " %s", argv[ ia ] );
            }
            std::fprintf( stdout, "\n" );

            // Where?
            std::fprintf( stdout, "     Run Dir     : %s\n", std::getenv( "PWD" ) );
            std::fprintf( stdout, "\n" );
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris

#endif /* SRC_CORE_BANNER_HPP_ */
