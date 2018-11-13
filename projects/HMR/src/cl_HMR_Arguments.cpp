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

                bool tArgumentsError = false;

                // assume refine step by default
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
                                std::cout << "No file path provided for parameters file." << std::endl;
                                tArgumentsError = true;
                                break;
                            }
                        }
                    }
                    else if
                    (   std::string( argv[ k ] ) == "--init"
                            || std::string( argv[ k ] ) == "-n" )
                    {
                        if( mState == State::MAP_FIELDS || mState == State::REFINE_MESH )
                        {
                            tArgumentsError = true;
                        }
                        else
                        {
                            mState = State::INITIALIZE_MESH;
                        }
                    }
                    else if
                    (   std::string( argv[ k ] ) == "--map"
                            || std::string( argv[ k ] ) == "-m" )
                    {
                        if( mState == State::INITIALIZE_MESH || mState == State::REFINE_MESH )
                        {
                            tArgumentsError = true;
                        }
                        else
                        {
                            mState = State::MAP_FIELDS;
                        }
                    }
                    else if
                    (   std::string( argv[ k ] ) == "--refine"
                            || std::string( argv[ k ] ) == "-r" )
                    {
                        if( mState == State::INITIALIZE_MESH || mState == State::MAP_FIELDS )
                        {
                            tArgumentsError = true;
                        }
                        else
                        {
                            mState = State::REFINE_MESH;
                        }
                    }

                }

                // detect invalid input
                if
                ( (   ( mState == State::REFINE_MESH || mState == State::INITIALIZE_MESH )
                  && ( mParameterPath.size() == 0 ) ) || tArgumentsError )
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
                //std::cout<< "--bincoeffs  <binaryfile> Dump coefficients into binary file      ( short -b )" << std::endl;
                //std::cout<< "--coeffs     <hdf5file>   Dump coefficients into hdf5 file        ( short -c )" << std::endl;
                //std::cout<< "--exodus     <exofile>    Dump output mesh into exodus file       ( short -e )" << std::endl;
                std::cout<< "--help                    Print this help screen                  ( short -h )" << std::endl;
                //std::cout<< "--in         <infile>     Load existing database from HDF5 file   ( short -i )" << std::endl;
                std::cout<< "--init                    Create a tensor field and quit          ( short -n )" << std::endl;
                std::cout<< "--laststep                Dump unrefined step into exodus         ( short -l )" << std::endl;
                std::cout<< "--map                     Map fields from input database to out   ( short -m )" << std::endl;
                //std::cout<< "--out        <outfile>    Save refined  datanbase into HDF5 file  ( short -o )" << std::endl;
                std::cout<< "--parameters <xmlfile>    Process parameters from <xmlfile>       ( short -p )" << std::endl;
                //std::cout<< "--timestep   <double>     Sets a timestep for the exo-file        ( short -t )" << std::endl;
                std::cout<< "--version                 Print banner and exit                   ( short -v )" << std::endl;
            }
        }

//---------------------------------------------------------------------------------
    }
}
