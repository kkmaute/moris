/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Parameters.cpp
 *
 */

#include <string>
#include "SDF_Tools.hpp"
#include "cl_SDF_Parameters.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    Parameter_List
    create_sdf_parameter_list()
    {
        Parameter_List aParameterList;

        aParameterList.insert( "input_mesh", std::string( "Mesh.exo" ) );
        aParameterList.insert( "output_mesh", std::string( "Mesh.exo" ) );

        return aParameterList;
    }

    //-------------------------------------------------------------------------------

    Parameter_List
    create_sdf_object_parameter_list()
    {
        Parameter_List aParameterList;
        aParameterList.insert( "label", std::string( "Object" ) );
        aParameterList.insert( "stl_file", std::string( "object.obj" ) );
        aParameterList.insert( "output_values", std::string( "" ) );
        aParameterList.insert( "output_hdf5", std::string( "" ) );
        aParameterList.insert( "candidate_search_depth", 2 );
        aParameterList.insert( "candidate_search_epsilon", 0.01 );

        return aParameterList;
    }

    //-------------------------------------------------------------------------------

    void
    load_sdf_parameter_list_from_xml(
            const std::string        &aFilePath,
            Parameter_List           &aGlobalParameters,
            Vector< Parameter_List > &aObjectParameters )
    {
        // create temporary Parser object
        XML_Parser tParser( aFilePath );

        std::string tValue;

        // global parameters
        aGlobalParameters = create_sdf_parameter_list();

        //            tParser.get( "moris.sdf.parameters.verbose", tValue );
        //            sint tSwitch = (sint) string_to_bool( tValue );
        //            aGlobalParameters.set("verbose",  tSwitch );

        // object parameters
        uint tNumberOfObjects = tParser.count_keys_in_subtree( "moris.sdf", "object" );

        // create cell
        aObjectParameters.clear();
        aObjectParameters.resize( tNumberOfObjects, create_sdf_object_parameter_list() );

        // loop over all objects
        for ( uint b = 0; b < tNumberOfObjects; ++b )
        {
            Vector< std::string > tFirst;
            Vector< std::string > tSecond;

            tParser.get_keys_from_subtree( "moris.sdf", "object", b, tFirst, tSecond );

            // copy key to settings struct
            for ( uint k = 0; k < tFirst.size(); ++k )
            {
                const std::string &tKey = tFirst( k );

                if ( tKey == "label" )
                {
                    aObjectParameters( b ).set( "label", std::string( tSecond( k ) ) );
                }
                else if ( tKey == "stl_file" )
                {
                    aObjectParameters( b ).set( "stl_file", std::string( tSecond( k ) ) );
                }
                else if ( tKey == "output_values" )
                {
                    aObjectParameters( b ).set( "output_values", std::string( tSecond( k ) ) );
                }
                else if ( tKey == "output_hdf5" )
                {
                    aObjectParameters( b ).set( "output_hdf5", std::string( tSecond( k ) ) );
                }
                else if ( tKey == "candidate_search_depth" )
                {
                    aObjectParameters( b ).set( "candidate_search_depth", (sint)std::stoi( tSecond( k ) ) );
                }
                else if ( tKey == "candidate_search_epsilon" )
                {
                    aObjectParameters( b ).set( "candidate_search_epsilon", (real)std::stod( tSecond( k ) ) );
                }
            }
        }
    }

    //-------------------------------------------------------------------------------
    }
