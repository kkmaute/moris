#include <cstdlib>

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src

#include "cl_MTK_Mapper_Arguments.hpp"


namespace moris
{
    namespace mapper
    {
// -----------------------------------------------------------------------------

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
                mState = State::MAP_FIELDS;

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
                }

                if( tArgumentsError )
                {
                    // print usage and exit
                    mState = State::PRINT_USAGE;
                }
            }
        }
// -----------------------------------------------------------------------------

        void
        Arguments::print_usage()
        {
            if( par_rank() == 0 )
            {
                std::cout << "Usage: mapper [option] <file>..." << std::endl;
                std::cout << std::endl;
                std::cout<< "run mapper --help to show options" << std::endl;
            }
        }

// -----------------------------------------------------------------------------

        void
        Arguments::print_help()
        {
            if( par_rank() == 0 )
            {
                std::cout << "Usage: mapper [option] <file> ..." << std::endl;
                std::cout << std::endl;
                std::cout<< "--help                    Print this help screen                  ( short -h )" << std::endl;
                std::cout<< "--parameters <xmlfile>    Process parameters from <xmlfile>       ( short -p )" << std::endl;
                std::cout<< "--version                 Print banner and exit                   ( short -v )" << std::endl;
            }
        }

// -----------------------------------------------------------------------------
    } /* namespace mapper */
} /* namespace moris */
