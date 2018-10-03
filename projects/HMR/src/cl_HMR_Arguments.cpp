#include <cstdlib>

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src

#include "HMR_Tools.hpp"
#include "cl_HMR_Arguments.hpp"

namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------

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
                // assume refine step by default
                mState = State::REFINE_MESH;

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
                                break;
                            }
                        }
                    }
                    else if (   std::string( argv[ k ] ) == "--in"
                             || std::string( argv[ k ] ) == "-i" )
                    {
                        if( k<argc-1 )
                        {
                            // return parameter path as output
                            mHdf5InputPath = parallelize_path( std::string( argv[ k+1 ] ) );
                        }
                        else
                        {
                            if( par_rank() == 0 )
                            {
                                std::cout << "No file path provided." << std::endl;
                                break;
                            }
                        }
                    }
                    else if (   std::string( argv[ k ] ) == "--out"
                            || std::string( argv[ k ] ) == "-o" )
                    {
                        if( k<argc-1 )
                        {
                            // return parameter path as output
                            mHdf5OutputPath = parallelize_path( std::string( argv[ k+1 ] ) );

                        }
                        else
                        {
                            if( par_rank() == 0 )
                            {
                                std::cout << "No file path provided." << std::endl;
                                break;
                            }
                        }
                    }
                    else if (   std::string( argv[ k ] ) == "--exodus"
                            || std::string( argv[ k ] ) == "-e" )
                    {
                        if( k<argc-1 )
                        {
                            // return parameter path as output
                            mExodusPath = std::string( argv[ k+1 ] );
                        }
                        else
                        {
                            if( par_rank() == 0 )
                            {
                                std::cout << "No file path provided." << std::endl;
                                break;
                            }
                        }

                    }
                    else if (   std::string( argv[ k ] ) == "--coeffs"
                            || std::string( argv[ k ] ) == "-c" )
                    {
                        if( k<argc-1 )
                        {
                            // return parameter path as output
                            mBinaryPath = std::string( argv[ k+1 ] );
                        }
                        else
                        {
                            if( par_rank() == 0 )
                            {
                                std::cout << "No file path provided." << std::endl;
                                break;
                            }
                        }
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
                                break;
                            }
                        }
                    }
                    else if
                    (   std::string( argv[ k ] ) == "--init"
                            || std::string( argv[ k ] ) == "-n" )
                    {
                        mState = State::INITIALIZE_MESH;
                    }
                }

                // detect invalid input
                if
                (    ( mState == State::REFINE_MESH || mState == State::INITIALIZE_MESH )
                  && ( mParameterPath.size() == 0 ) )
                {
                    mState = State::PRINT_USAGE;
                }
            }
        }

//---------------------------------------------------------------------------------


        void
        Arguments::print_usage()
        {
            if( par_rank() == 0 )
            {
                std::cout << "Usage: hmr [option] <file>..." << std::endl;
                std::cout << std::endl;
                std::cout<< "run hmr --help to show options" << std::endl;
            }
        }

//---------------------------------------------------------------------------------

        void
        Arguments::print_help()
        {
            if( par_rank() == 0 )
            {
                std::cout << "Usage: hmr [option] <file> ..." << std::endl;
                std::cout << std::endl;
                std::cout<< "Options:" << std::endl;

                std::cout<< "--coeffs     <binaryfile> Dump coefficients into binary file ( short -c )" << std::endl;
                std::cout<< "--exodus     <exofile>    Dump output mesh into exodus file  ( short -e )" << std::endl;
                std::cout<< "--in         <infile>     Load existing mesh from HDF5 file  ( short -i )" << std::endl;
                std::cout<< "--init                    Create a tensor field and quit     ( short -n )" << std::endl;
                std::cout<< "--help                    Print this help screen             ( short -h )" << std::endl;
                std::cout<< "--out        <infile>     Save refined  mesh into HDF5 file  ( short -o )" << std::endl;
                std::cout<< "--parameters <xmlfile>    Process parameters from <xmlfile>  ( short -p )" << std::endl;
                std::cout<< "--timestep   <double>     Sets a timestep for the exo-file   ( short -t )" << std::endl;
                std::cout<< "--version                 Print banner and exit              ( short -v )" << std::endl;
            }
        }

//---------------------------------------------------------------------------------
    }
}
