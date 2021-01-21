#include "cl_Library_IO.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    Library_IO::Library_IO( const std::string & aPath ) : mPath( std::getenv( "PWD" ) )
    {
        // get first letter of aPath
        if( aPath.at( 0 ) == '/' )
        {
            // this is an absolute path
            mPath = aPath;
        }
        else
        {
            // this is a relative path
            mPath = mPath + "/" + aPath;
        }

        // try to open library file
        mLibraryHandle = dlopen( mPath.c_str(), RTLD_NOW );

        // test if loading succeeded
        if( ! mLibraryHandle )
        {
            // get error string
            std::string tError = dlerror();

            // throw error
            MORIS_ERROR( mLibraryHandle, tError.c_str() );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO::~Library_IO()
    {
        // close handle to library
        dlclose( mLibraryHandle );
    }

    //------------------------------------------------------------------------------------------------------------------
}