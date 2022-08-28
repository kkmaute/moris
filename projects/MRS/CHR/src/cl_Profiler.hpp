/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Profiler.hpp
 *
 */

#ifndef PROJECTS_MRS_CHR_SRC_CL_PROFILER_HPP_
#define PROJECTS_MRS_CHR_SRC_CL_PROFILER_HPP_

#include <iostream>
#include <string>

#ifdef WITHGPERFTOOLS
#include "gperftools/profiler.h"
#endif

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src

/*
 * Instruction :
 *
 * - checkout google  performance tools from
 *   https://github.com/gperftools/gperftools
 *
 *   install using the commands
 *   sh autogen.sh
 *   ./configure --prefix= ${APPS}/gperftools
 *   make
 *   make install
 *
 * - install kcachegrind
 *
 * - have USE_GPERFTOOLS on during compilation
 *
 * - create profile using
 *    Profiler tProf;
 *
 *    < code you want to profile >
 *
 *    tProf.stop();
 *
 * - after run, call
 *   kcachegrind /tmp/gprofmoris.callgrind"
 * /
 */
namespace moris
{
    class Profiler
    {
        std::string mLogFile;
        std::string mCallgrindFile;
//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        Profiler( const std::string & aLogFile )
        {
#ifdef WITHGPERFTOOLS
            // fixme: replace this line with moris logger
            std::cout << "Starting MORIS Profiler ..." << std::endl;

            // grab path without file extension
            std::string tBasepath = aLogFile.substr( 0, aLogFile.find_last_of(".") );

            if( par_size() == 1 )
            {
                mLogFile       = aLogFile;
                mCallgrindFile = tBasepath + ".callgrind";
            }
            else
            {
                // make paths parallel
                std::string tSuffix = "." + std::to_string( par_size() ) + "." + std::to_string( par_rank() );

                // path for parallel logfile
                mLogFile       = tBasepath
                        + tSuffix + aLogFile.substr( aLogFile.find_last_of("."), aLogFile.length() );

                // path for callgrind file
                mCallgrindFile =  tBasepath + tSuffix + ".callgrind";
            }

            // start google profiler
            ProfilerStart( mLogFile.c_str() );
#endif
        }

//------------------------------------------------------------------------------

        Profiler() : Profiler( "/tmp/gprofmoris.log" )
        {

        }

//------------------------------------------------------------------------------

        void
        stop()
        {
#ifdef WITHGPERFTOOLS

            // stop the google profiler
            ProfilerStop();

            // get path to moris executable
            const std::string & tExec = gMorisComm.get_exec_path();

            // get app path
            std::string tApps = getenv("APPS");

            // assemble command line
            std::string tCommand = tApps + "/gperftools/bin/pprof --callgrind "
                    + tExec + " " + mLogFile + " > " + mCallgrindFile;

            // fixme: replace this line with moris logger
            std::cout << "Creating callgrind file " <<  mCallgrindFile << " ..." << std::endl;

            // convert logfile to callgrind
            system( tCommand.c_str() );

            // fixme: replace this line with moris logger
            std::cout << "... Stopped MORIS Profiler" << std::endl;
#endif
        }

//------------------------------------------------------------------------------
    };

}

#endif /* PROJECTS_MRS_CHR_SRC_CL_PROFILER_HPP_ */

