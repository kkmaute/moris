#ifndef MORIS_CL_LIBRARY_IO_HPP
#define MORIS_CL_LIBRARY_IO_HPP

#include <string>
#include "dlfcn.h"
#include "assert.hpp"

namespace moris
{
    /**
     * Wrapper class for shared object file
     */
    class Library_IO
    {
        // path to library file
        std::string mPath;

        // handle to shared object
        void *  mLibraryHandle;

    public:

        /**
         * Constructor
         *
         * @param aPath Path to shared object file
         */
        Library_IO( const std::string & aPath );

        /**
         * Destructor
         */
        ~Library_IO();

        /**
         * Loads a function from the shared object file.
         *
         * @tparam Function_Type Function pointer type
         * @param aFunctionName Function name to look for in the file
         * @param aThrowError parameter to check if the list exists
         * @return Function pointer
         */
        template <typename Function_Type>
        Function_Type load_function(std::string aFunctionName, bool aThrowError = true)
        {
            Function_Type aUserFunction
                    = reinterpret_cast<Function_Type>
                    ( dlsym( mLibraryHandle, aFunctionName.c_str() ) );

            //depending on the flag throw an error
            if( aThrowError )
            {
                // create error message
                std::string tError = "Could not find symbol " + aFunctionName
                                     + "  within file " + mPath;

                // make sure that loading succeeded
                MORIS_ERROR( aUserFunction, tError.c_str() );
            }

            // return function handle
            return aUserFunction;
        }
    };
}

#endif //MORIS_CL_LIBRARY_IO_HPP
