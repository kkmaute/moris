/*
 * fn_hmr_load_user_library.hpp
 *
 *  Created on: Nov 13, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_LOAD_USER_LIBRARY_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_LOAD_USER_LIBRARY_HPP_

#include <string>

// dynamik linker function
#include "dlfcn.h"

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "HMR_Globals.hpp"
#include "cl_HMR_Element.hpp"

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        /**
         * Interface for user defined function
         */
        typedef bool ( *MORIS_HMR_USER_FUNCTION )
                (
                        const Element                  * aElement,
                        const Cell< Matrix< DDRMat > > & aElementLocalValues,
                        ParameterList                  & aParameters
                );

// -----------------------------------------------------------------------------

        /**
         * Wrapper class for shared object file
         */
        class Library
        {
            // path to library file
            const std::string mPath;

            // handle to shared object
            void *  mLibraryHandle;
// -----------------------------------------------------------------------------
        public:
// -----------------------------------------------------------------------------

            Library( const std::string & aPath ) : mPath( aPath )
            {
                // try to open library file
                mLibraryHandle = dlopen( aPath.c_str(), RTLD_NOW );

                // test if loading succeeded
                if( ! mLibraryHandle )
                {
                    // get error string
                    std::string tError = dlerror();

                    // throw error
                    MORIS_ERROR( mLibraryHandle, tError.c_str() );
                }
            }

// -----------------------------------------------------------------------------

            ~Library()
            {
                // close handle to library
                dlclose( mLibraryHandle );
            }

// -----------------------------------------------------------------------------

            MORIS_HMR_USER_FUNCTION
            load_function( const std::string & aFunctionName )
            {
                MORIS_HMR_USER_FUNCTION aUserFunction
                    = reinterpret_cast<MORIS_HMR_USER_FUNCTION>
                    ( dlsym( mLibraryHandle, aFunctionName.c_str() ) );

                // create error message
                std::string tError =  "Could not find symbol " + aFunctionName
                        + "  within file " + mPath;

                // make sure that loading succeeded
                MORIS_ERROR( aUserFunction, tError.c_str() );

                // return function handle
                return aUserFunction;
            }

// -----------------------------------------------------------------------------
        };
// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */



#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_LOAD_USER_LIBRARY_HPP_ */
