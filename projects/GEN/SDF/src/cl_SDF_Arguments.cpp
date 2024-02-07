/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Arguments.cpp
 *
 */

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src

#include "cl_SDF_Arguments.hpp"

namespace moris
{
    namespace sdf
    {
        Arguments::Arguments(
                int  & argc,
                char * argv[] )
                {

            if( argc == 1 )
            {
                // print usage and exit
                mState = State::PRINT_USAGE;
            }
            else
            {
                bool tSdfSwitch      = false;
                bool tArgumentsError = false;

                // assume usage switch by default
                mState = State::PRINT_USAGE;

                // loop over all arguments
                for( int k=0; k<argc; ++k )
                {
                    if (   std::string( argv[ k ] ) == "--version"
                            || std::string( argv[ k ] ) == "-v" )
                    {
                        mState = State::PRINT_VERSION;
                        break;
                    }
                    else if ( ( std::string( argv[ k ] ) == "--help" )
                            || std::string( argv[ k ] ) == "-h" )
                    {
                        mState = State::PRINT_HELP;
                        break;
                    }
                    else if (
                            std::string( argv[ k ] ) == "--parameters"
                                    || std::string( argv[ k ] ) == "-p" )
                    {
                        if( k<argc-1 )
                        {
                            // return parameter path as output
                            mParameterPath = std::string( argv[ k+1 ] );
                        }
                        else
                        {
                            if( par_rank() == 0 )
                            {
                                std::cout << "No file path provided." << std::endl;
                                tArgumentsError = true;
                                break;
                            }
                        }
                    }
                    else if (
                            std::string( argv[ k ] ) == "--mesh"
                                    || std::string( argv[ k ] ) == "-m" )
                    {
                        if( k<argc-1 )
                        {
                            // return parameter path as output
                            mInputMeshPath = std::string( argv[ k+1 ] );
                        }
                        else
                        {
                            if( par_rank() == 0 )
                            {
                                std::cout << "No file path for input mesh provided." << std::endl;
                                tArgumentsError = true;
                                break;
                            }
                        }
                    }
                    else if (
                            std::string( argv[ k ] ) == "--exodus"
                                    || std::string( argv[ k ] ) == "-e" )
                    {
                        if( k<argc-1 )
                        {
                            // return parameter path as output
                            mOutputMeshPath = std::string( argv[ k+1 ] );
                        }
                        else
                        {
                            if( par_rank() == 0 )
                            {
                                std::cout << "No file path for output mesh provided." << std::endl;
                                tArgumentsError = true;
                                break;
                            }
                        }
                    }
                    else if (
                            std::string( argv[ k ] ) == "--raycast"
                                    || std::string( argv[ k ] ) == "-r" )
                    {
                        if ( ! tSdfSwitch )
                        {
                            mState = State::CALCULATE_RAYCAST;
                        }
                    }
                    else if (
                            std::string( argv[ k ] ) == "--sdf"
                                    || std::string( argv[ k ] ) == "-s" )
                    {
                        mState = State::CALCULATE_RAYCAST_AND_SDF;
                        tSdfSwitch = true;
                    }
                    else if (   std::string( argv[ k ] ) == "--timestep"
                            || std::string( argv[ k ] ) == "-t" )
                    {
                        if( k<argc-1 )
                        {
                            // return parameter path as output
                            mTimestep = std::atof( argv[ k+1 ] );
                        }
                        else
                        {
                            if( par_rank() == 0 )
                            {
                                std::cout << "No timestep provided." << std::endl;
                                tArgumentsError = true;
                                break;
                            }
                        }
                    }
                }
                if( ( (    mState == State::CALCULATE_RAYCAST
                        || mState == State::CALCULATE_RAYCAST_AND_SDF )
                        && mParameterPath.size() == 0 ) || tArgumentsError )
                {
                    mState = State::PRINT_USAGE;
                }
            }
                }
        //--------------------------------------------------------------------------------

        void
        Arguments::print_usage()
        {
            if( par_rank() == 0 )
            {
                std::cout << "Usage: sdf [option] <file>..." << std::endl;
                std::cout << std::endl;
                std::cout<< "run sdf --help to show options" << std::endl;
            }
        }

        //--------------------------------------------------------------------------------

        void
        Arguments::print_help()
        {
            if( par_rank() == 0 )
            {
                std::cout << "Usage: sdf [option] <file> ..." << std::endl;
                std::cout << std::endl;
                std::cout<< "Options:" << std::endl;

                std::cout<< "--exodus                  Path to output mesh                     ( short -e )" << std::endl;
                std::cout<< "--help                    Print this help screen                  ( short -h )" << std::endl;
                std::cout<< "--mesh                    Path to background mesh                 ( short -m )" << std::endl;
                std::cout<< "--parameters <xmlfile>    Process parameters from <xmlfile>       ( short -p )" << std::endl;
                std::cout<< "--raycast                 Perform raycast only                    ( short -r )" << std::endl;
                std::cout<< "--sdf                     Perform raycast and calculate sdf       ( short -s )" << std::endl;
                std::cout<< "--timestep   <double>     Sets a timestep for the exo-file        ( short -t )" << std::endl;
                std::cout<< "--version                 Print banner and exit                   ( short -v )" << std::endl;
            }
        }

        //--------------------------------------------------------------------------------
    }
}

