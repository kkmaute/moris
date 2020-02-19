/*
 * fn_hmr_load_user_library.hpp
 *
 *  Created on: Okt 13, 2019
 *      Author: schmidt
 */

#ifndef PROJECTS_HMR_SRC_FN_EXEC_LOAD_USER_LIBRARY_HPP_
#define PROJECTS_HMR_SRC_FN_EXEC_LOAD_USER_LIBRARY_HPP_

#include <string>

// dynamic linker function
#include "dlfcn.h"

#include "HMR_Globals.hpp"
#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
    namespace tsa
    {
        class Time_Solver;
    }

    namespace fem
    {
        class Field_Interpolator_Manager;
    }

// -----------------------------------------------------------------------------

        /**
         * Interface for user defined function
         */
        typedef void ( *MORIS_USER_FUNCTION ) ( Cell< moris::real > & aParameter );

        typedef bool ( *MORIS_SOL_CRITERIA_FUNC ) ( moris::tsa::Time_Solver * aTimeSolver );

        typedef void ( *MORIS_PARAMETER_FUNCTION ) ( moris::Cell< moris::Cell< moris::ParameterList > > & aParameterList );

        typedef void ( *MORIS_FEM_FREE_FUNCTION ) ( moris::Matrix< moris::DDRMat >                & aPropMatrix,
                                                    moris::Cell< moris::Matrix< moris::DDRMat > > & aParameters,
                                                    moris::fem::Field_Interpolator_Manager        * aFIManager );

// -----------------------------------------------------------------------------

        /**
         * Wrapper class for shared object file
         */
        class Library_IO
        {
            // path to library file
            std::string mPath;

            // handle to shared object
            void *  mLibraryHandle;
// -----------------------------------------------------------------------------
        public:
// -----------------------------------------------------------------------------

            Library_IO( const std::string & aPath ) : mPath( std::getenv( "PWD" ) )
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

// -----------------------------------------------------------------------------

            ~Library_IO()
            {
                // close handle to library
                dlclose( mLibraryHandle );
            }

// -----------------------------------------------------------------------------

            MORIS_USER_FUNCTION
            load_function( const std::string & aFunctionName )
            {
                MORIS_USER_FUNCTION aUserFunction
                    = reinterpret_cast<MORIS_USER_FUNCTION>
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

            MORIS_PARAMETER_FUNCTION
            load_parameter_file( const std::string & aFunctionName )
            {
                MORIS_PARAMETER_FUNCTION aUserFunction
                    = reinterpret_cast<MORIS_PARAMETER_FUNCTION>
                    ( dlsym( mLibraryHandle, aFunctionName.c_str() ) );

//                // create error message
                std::string tError =  "Could not find symbol " + aFunctionName
                        + "  within file " + mPath;

                // make sure that loading succeeded
                MORIS_ERROR( aUserFunction, tError.c_str() );

                // return function handle
                return aUserFunction;
            }

// -----------------------------------------------------------------------------

            MORIS_FEM_FREE_FUNCTION
            load_fem_free_functions( const std::string & aFunctionName )
            {
                MORIS_FEM_FREE_FUNCTION aUserFunction
                    = reinterpret_cast<MORIS_FEM_FREE_FUNCTION>
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

            MORIS_SOL_CRITERIA_FUNC
            load_sol_criteria_functions( const std::string & aFunctionName )
            {
                MORIS_SOL_CRITERIA_FUNC aUserFunction
                    = reinterpret_cast<MORIS_SOL_CRITERIA_FUNC>
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

    } /* namespace moris */



#endif /* PROJECTS_HMR_SRC_FN_EXEC_LOAD_USER_LIBRARY_HPP_ */
