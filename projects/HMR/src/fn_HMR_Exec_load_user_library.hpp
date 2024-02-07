/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_HMR_Exec_load_user_library.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_LOAD_USER_LIBRARY_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_LOAD_USER_LIBRARY_HPP_

#include <string>

// dynamic linker function
#include "dlfcn.h"

#include "cl_HMR_Element.hpp"
#include "HMR_Globals.hpp"
#include "assert.hpp"
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris::hmr
{
    // -----------------------------------------------------------------------------

    /**
     * Interface for user defined function
     */
    typedef int ( *MORIS_HMR_USER_FUNCTION )(
            Element                        *aElement,
            const Vector< Matrix< DDRMat > > &aElementLocalValues,
            ParameterList                  &aParameters );

    // -----------------------------------------------------------------------------

    /**
     * Wrapper class for shared object file
     */
    class Library
    {
        // path to library file
        std::string mPath;

        // handle to shared object
        void *mLibraryHandle;
        // -----------------------------------------------------------------------------

      public:
        // -----------------------------------------------------------------------------

        Library( const std::string &aPath )
                : mPath( std::getenv( "PWD" ) )
        {
            // get first letter of aPath
            if ( aPath.at( 0 ) == '/' )
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
            if ( !mLibraryHandle )
            {
                // get error string
                std::string tError = dlerror();

                // throw error
                MORIS_ERROR( mLibraryHandle, "Library::Library - %s", tError.c_str() );
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
        load_function( const std::string &aFunctionName )
        {
            MORIS_HMR_USER_FUNCTION aUserFunction = reinterpret_cast< MORIS_HMR_USER_FUNCTION >( dlsym( mLibraryHandle, aFunctionName.c_str() ) );

            // make sure that loading succeeded
            MORIS_ERROR( aUserFunction, "Could not find symbol %s within file %s", aFunctionName.c_str(), mPath.c_str() );

            // return function handle
            return aUserFunction;
        }

        // -----------------------------------------------------------------------------
    };
    // -----------------------------------------------------------------------------
}    // namespace moris::hmr

#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_LOAD_USER_LIBRARY_HPP_ */
