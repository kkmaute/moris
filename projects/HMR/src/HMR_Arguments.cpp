#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "HMR_Arguments.hpp"

namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------

        void
        process_arguments(
                int    argc,
                char * argv[],
                std::string & aParameterPath,
                std::string & aInputPath,
                std::string & aOutputPath,
                std::string & aExodusPath )
        {
            // return empty path by default
            aParameterPath = "";
            aInputPath = "";
            aOutputPath = "";
            aExodusPath = "";

            if( argc == 1 )
            {
                // print usage and exit
                print_usage();
            }
            else
            {
                // loop over all arguments
                for( int k=0; k<argc; ++k )
                {
                    if ( std::string( argv[ k ] ) == "--version")
                    {
                        moris::print_banner( argc, argv );
                        break;
                    }
                    else if ( std::string( argv[ k ] ) == "--help" )
                    {
                        print_help();
                        break;
                    }
                    else if (
                            std::string( argv[ k ] ) == "--parameters"
                         || std::string( argv[ k ] ) == "-p" )
                    {
                        if( k<argc-1 )
                        {
                            // return parameter path as output
                            aParameterPath = std::string( argv[ k+1 ] );
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
                            aInputPath = std::string( argv[ k+1 ] );
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
                            aOutputPath = std::string( argv[ k+1 ] );
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
                            aExodusPath = std::string( argv[ k+1 ] );
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
                }
            }
        }

//---------------------------------------------------------------------------------


        void
        print_usage()
        {
            if( par_rank() == 0 )
            {
                std::cout << "Usage: hmr [options] file..." << std::endl;
                std::cout << std::endl;
                std::cout<< "run hmr --help to show options" << std::endl;
            }
        }

//---------------------------------------------------------------------------------

        void
        print_help()
        {
            if( par_rank() == 0 )
            {
                std::cout << "Usage: hmr [options] file..." << std::endl;
                std::cout<< "Options:" << std::endl;

                std::cout<< "--exodus     <exofile>   Dump output mesh into exodus file ( short -e )" << std::endl;
                std::cout<< "--in         <infile>    Load existing mesh from HDF5 file ( short -i )" << std::endl;
                std::cout<< "--out        <infile>    Save refined  mesh into HDF5 file ( short -o )" << std::endl;
                std::cout<< "--parameters <xmlfile>   Process parameters from <xmlfile> ( short -p )" << std::endl;
                std::cout<< "--version                Print banner and exit" << std::endl;
            }
        }

//---------------------------------------------------------------------------------
    }
}
