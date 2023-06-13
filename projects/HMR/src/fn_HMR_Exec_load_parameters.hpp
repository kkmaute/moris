/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_HMR_Exec_load_parameters.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_LOAD_PARAMETERS_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_LOAD_PARAMETERS_HPP_

#include <string>

#include "HMR_Globals.hpp"
#include "HMR_Tools.hpp"
#include "typedefs.hpp" //COR/src

// -----------------------------------------------------------------------------

namespace moris::hmr
{
// -----------------------------------------------------------------------------

    void
    load_file_list_from_xml(
            const std::string & aFilePath,
            ParameterList     & aFileList )
    {

        Cell< std::string > tKeys;

        tKeys.push_back( "input_mesh_database" );
        tKeys.push_back( "output_mesh_database" );
        tKeys.push_back( "last_step_exodus" );
        tKeys.push_back( "refined_exodus" );
        tKeys.push_back( "output_exodus" );
        tKeys.push_back( "output_field_database" );
        tKeys.push_back( "coefficient_database" );
        tKeys.push_back( "output_exodus_order" );
        // create parameter list
        for( std::string tKey : tKeys )
        {
            if( tKey == "output_exodus_order" )
            {
                aFileList.insert( tKey, ( sint ) 0 );
            }
            else
            {
                aFileList.insert( tKey, std::string( "" )  );
            }
        }

        // create temporary Parser object
        XML_Parser tParser( aFilePath );
        Cell< std::string > tFirst;
        Cell< std::string > tSecond;

        // fill parser with values
        tParser.get_keys_from_subtree( "moris.hmr", "files", 0, tFirst, tSecond );

        uint tCount = 0;

        // loop over all entries and set parameters
        for( std::string tKey : tFirst )
        {
            if( tKey == "output_exodus_order" )
            {
                aFileList.set( tKey, ( sint ) stoi( tSecond( tCount++ ) ) );
            }
            else
            {
                aFileList.set( tKey, tSecond( tCount++ )  );
            }
        }

    }

// -----------------------------------------------------------------------------

    void
    load_refinement_parameters_from_xml(
            const std::string & aFilePath,
            ParameterList     & aRefParams )
    {
        // create temporary Parser object
        XML_Parser tParser( aFilePath );
        Cell< std::string > tFirst;
        Cell< std::string > tSecond;
        tParser.get_keys_from_subtree( "moris.hmr", "refinement", 0, tFirst, tSecond );

        aRefParams.insert( "library", "" );
        aRefParams.insert( "function", "" );
        uint tCount=0;
        for( std::string tKey : tFirst )
        {

            if( tKey == "library" )
            {
                aRefParams.set( "library", tSecond( tCount++ ) );
            }
            else if( tKey == "function" )
            {
                aRefParams.set( "function", tSecond( tCount++ ) );
            }
            else
            {
                std::string tType  = tKey.substr( 0, tKey.find_first_of("_") );
                std::string tLabel = tKey.substr( tKey.find_first_of("_")+1, tKey.size() );

                if( tType == "real" )
                {
                    aRefParams.insert( tLabel, ( moris::real ) stod( tSecond( tCount++ ) ) );
                }
                else if( tType == "int" )
                {
                    aRefParams.insert( tLabel, ( moris::sint ) stoi( tSecond( tCount++ ) ) );
                }
                else if ( tType == "bool" )
                {
                    aRefParams.insert( tLabel, ( moris::sint ) string_to_bool( tSecond( tCount++ ) ) );
                }
                else if ( tType == "string" )
                {
                    aRefParams.insert( tLabel, tSecond( tCount++ ) );
                }
                else
                {
                    MORIS_ERROR( false, "custom parameters in refinement tag must begin with real, int, bool or string" );
                }
            }
        }
    }

// -----------------------------------------------------------------------------

    void
    load_field_parameters_from_xml(
            const std::string       & aFilePath,
            Cell< ParameterList >   & aFieldParams )
    {
        // clean up output
        aFieldParams.clear();

        // create a parser object for this settings path
        XML_Parser tParser( aFilePath );

        // create empty parameter list
        ParameterList tParams;
        tParams.insert( "label", "untitled" );
        tParams.insert( "lagrange_order", (sint) 0 );
        tParams.insert( "bspline_order", (sint) 0 );
        tParams.insert( "source", "" );
        tParams.insert( "perform_mapping", (sint) 1 );

        tParams.insert( "min_volume_refinement_level", (sint) 0 );
        tParams.insert( "max_volume_refinement_level", (sint) gMaxNumberOfLevels );
        tParams.insert( "min_surface_refinement_level", (sint) 0 );
        tParams.insert( "max_surface_refinement_level", (sint) gMaxNumberOfLevels );

        // get number of fields from parser
        uint tNumberOfFields = tParser.count_keys_in_subtree( "moris.hmr", "field" );

        aFieldParams.resize( tNumberOfFields, tParams );

        // loop over all fields
        for( uint f=0; f<tNumberOfFields; ++f )
        {
            Cell< std::string > tFirst;
            Cell< std::string > tSecond;

            tParser.get_keys_from_subtree( "moris.hmr", "field", f, tFirst, tSecond );

            // copy key to settings struct
            for( uint k=0; k<tFirst.size(); ++k )
            {
                std::string tKey = tFirst( k );

                if ( tKey == "label" )
                {
                    aFieldParams( f ).set( "label", std::string(  tSecond( k ) ) );
                }
                else if ( tKey == "lagrange_order" )
                {
                    aFieldParams( f ).set( "lagrange_order", (sint) stoi( tSecond( k ) ) );

                }
                else if ( tKey == "bspline_order" )
                {
                    aFieldParams( f ).set( "bspline_order", (sint) stoi( tSecond( k ) ) );

                }
                else if ( tKey == "source" )
                {
                    aFieldParams( f ).set( "source",  tSecond( k ));
                }
                else if ( tKey == "min_volume_refinement_level" )
                {
                    aFieldParams( f ).set( "min_volume_refinement_level", (sint) stoi( tSecond( k ) ) );

                }
                else if ( tKey == "max_volume_refinement_level" )
                {
                    aFieldParams( f ).set( "max_volume_refinement_level", (sint) stoi( tSecond( k ) ) );

                }
                else if ( tKey == "min_surface_refinement_level" )
                {
                    aFieldParams( f ).set( "min_surface_refinement_level", (sint) stoi( tSecond( k ) ) );

                }
                else if ( tKey == "max_surface_refinement_level" )
                {
                    aFieldParams( f ).set( "max_surface_refinement_level", (sint) stoi( tSecond( k ) ) );

                }
                else if ( tKey == "perform_mapping" )
                {
                    aFieldParams( f ).set( "perform_mapping", (sint) string_to_bool( tSecond( k ) ) );

                }
            }
        }
    }

// -----------------------------------------------------------------------------
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_LOAD_PARAMETERS_HPP_ */

