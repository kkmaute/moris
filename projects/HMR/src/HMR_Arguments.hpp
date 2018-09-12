#include <string>

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
                std::string & aExodusPath );

//---------------------------------------------------------------------------------

        void
        print_usage();

//---------------------------------------------------------------------------------

        void
        print_help();

//---------------------------------------------------------------------------------
    }
}
