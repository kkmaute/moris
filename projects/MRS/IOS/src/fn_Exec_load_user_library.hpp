/*
 * fn_Exec_load_user_library.hpp
 *
 *  Created on: Okt 13, 2019
 *      Author: schmidt
 */

#ifndef PROJECTS_HMR_SRC_FN_EXEC_LOAD_USER_LIBRARY_HPP_
#define PROJECTS_HMR_SRC_FN_EXEC_LOAD_USER_LIBRARY_HPP_

#include <string>

// dynamic linker function
#include "dlfcn.h"

//#include "HMR_Globals.hpp"
#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Param_List.hpp"

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

    namespace hmr
    {
        class Element;
    }

    // -----------------------------------------------------------------------------

    /**
     * Interface for user defined function
     */
    typedef Matrix<DDRMat> ( *MORIS_CRITERIA_FUNCTION ) (const moris::Matrix<DDRMat>& );

    typedef Matrix<DDRMat> ( *MORIS_OBJECTIVE_CONSTRAINT_FUNCTION ) (
            const moris::Matrix<DDRMat>&,
            const moris::Matrix<DDRMat>& );

    typedef void ( *MORIS_CRITERIA_INITIALIZE_FUNCTION ) (
            moris::Matrix<DDRMat>& aAdvs,
            moris::Matrix<DDRMat>& aLowerBounds,
            moris::Matrix<DDRMat>& aUpperBounds);

    typedef Matrix<DDSMat> ( *MORIS_CONSTRAINT_TYPES_FUNCTION ) ( );

    typedef bool ( *MORIS_SOL_CRITERIA_FUNC ) ( moris::tsa::Time_Solver * aTimeSolver );

    typedef bool ( *MORIS_POINTER_FUNC ) ( void * aPointer );

    typedef void ( *MORIS_PARAMETER_FUNCTION ) ( moris::Cell< moris::Cell< moris::ParameterList > > & aParameterList );

    typedef real ( *MORIS_GEN_FIELD_FUNCTION ) (
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aParameters );

    typedef void ( *MORIS_GEN_SENSITIVITY_FUNCTION ) (
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< moris::real* >& aParameters,
            moris::Matrix< DDRMat >&           aReturnValue);

    typedef void ( *MORIS_FEM_FREE_FUNCTION ) (
            moris::Matrix< moris::DDRMat >                & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > > & aParameters,
            moris::fem::Field_Interpolator_Manager        * aFIManager );

    typedef sint  ( *MORIS_USER_DEFINED_REFINEMENT_FUNCTION ) (
            hmr::Element                   * aElement,
            const Matrix< DDRMat > & aElementLocalValues);

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

         MORIS_CRITERIA_FUNCTION
         load_criteria_function( const std::string & aFunctionName )
         {
             MORIS_CRITERIA_FUNCTION aUserFunction
             = reinterpret_cast<MORIS_CRITERIA_FUNCTION>
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

         MORIS_OBJECTIVE_CONSTRAINT_FUNCTION
         load_objective_constraint_function( const std::string & aFunctionName )
         {
             MORIS_OBJECTIVE_CONSTRAINT_FUNCTION aUserFunction
             = reinterpret_cast<MORIS_OBJECTIVE_CONSTRAINT_FUNCTION>
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

         MORIS_CRITERIA_INITIALIZE_FUNCTION
         load_criteria_initialize_function( const std::string & aFunctionName )
         {
             MORIS_CRITERIA_INITIALIZE_FUNCTION aUserFunction
             = reinterpret_cast<MORIS_CRITERIA_INITIALIZE_FUNCTION>
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

         MORIS_CONSTRAINT_TYPES_FUNCTION
         load_constraint_types_function(const std::string & aFunctionName )
         {
             MORIS_CONSTRAINT_TYPES_FUNCTION aUserFunction
             = reinterpret_cast<MORIS_CONSTRAINT_TYPES_FUNCTION>
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

         MORIS_GEN_FIELD_FUNCTION
         load_gen_field_function( const std::string & aFunctionName )
         {
             MORIS_GEN_FIELD_FUNCTION aUserFunction
             = reinterpret_cast<MORIS_GEN_FIELD_FUNCTION>
             ( dlsym( mLibraryHandle, aFunctionName.c_str() ) );

             // create error message
             std::string tError =  "Could not find symbol " + aFunctionName
                     + "  within file " + mPath;

             // make sure that loading succeeded
             MORIS_ERROR( aUserFunction, tError.c_str() );

             // return function handle
             return aUserFunction;
         }

         MORIS_GEN_SENSITIVITY_FUNCTION
         load_gen_sensitivity_function( const std::string & aFunctionName )
         {
             MORIS_GEN_SENSITIVITY_FUNCTION aUserFunction
             = reinterpret_cast<MORIS_GEN_SENSITIVITY_FUNCTION>
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

         MORIS_POINTER_FUNC
         load_pointer_functions( const std::string & aFunctionName )
         {
             MORIS_POINTER_FUNC aUserFunction
             = reinterpret_cast<MORIS_POINTER_FUNC>
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

         MORIS_USER_DEFINED_REFINEMENT_FUNCTION
         load_user_defined_refinement_functions( const std::string & aFunctionName )
         {
             MORIS_USER_DEFINED_REFINEMENT_FUNCTION aUserFunction
             = reinterpret_cast<MORIS_USER_DEFINED_REFINEMENT_FUNCTION>
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
