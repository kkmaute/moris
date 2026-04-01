/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_STK_Parameters.hpp
 *
 */

#pragma once

#include "fn_STK_Enums.hpp"
#include "cl_Parameter_List.hpp"

namespace moris::prm
{
    // creates a parameter list with default inputs
    inline Parameter_List
    create_stk_parameter_list()
    {
        Parameter_List tParameterList( "STK" );

        // decomposition and decomposition related parameters
        tParameterList.insert( "input_file", "" );
        tParameterList.insert( "periodic_workspace", false );
        tParameterList.insert( "periodic_side_set_pair", "" );
        tParameterList.insert( "requested_QIs", Vector< std::string >() );

        return tParameterList;
    }

    /**
     * Creates a basic parameter list for SQIs with no specific parameters set.
     *
     * @return SQI parameter list. Cannot necessarily be used directly, needs additional parameters inserted by insert_SQI_parameters()
     */
    static Parameter_List create_SQI_parameter_list()
    {
        Parameter_List tParameterList( "SQI" );

        tParameterList.insert( "mesh_set_names", "" );                             // Name of sidesets to create surface mesh from
        tParameterList.insert_enum( "SQI_type", sqi::SQI_Type_String::values );    // Type of SQI to be computed
        tParameterList.insert( "SQI_name", "" );                                   // Name of the SQI, used for choosing design criteria for OPT or for output

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    static void insert_SQI_parameters( Parameter_List& aSQIParameterList, sqi::SQI_Type aSQIType )
    {
        switch ( aSQIType )
        {
            case sqi::SQI_Type::VOLUME:
                aSQIParameterList.set( "SQI_type", sqi::SQI_Type::VOLUME );
                break;
            case sqi::SQI_Type::RAYCAST_SHAPE_DIAMETER:
                aSQIParameterList.set( "SQI_type", sqi::SQI_Type::RAYCAST_SHAPE_DIAMETER );
                aSQIParameterList.insert( "number_of_polar_rays", 20, 1, 1000 );                               // Number of rays to be cast in the polar direction for the shape diameter function
                aSQIParameterList.insert( "number_of_azimuth_rays", 1, 1, 1000 );                              // Number of rays to be cast in the azimuth direction for the shape diameter function
                aSQIParameterList.insert( "cone_angle", 30.0, 0.0, 179.9999999 );                              // Cone angle in degrees for the shape diameter function
                aSQIParameterList.insert( "agglomeration_exponent", 1.0, 1.0, 1000.0 );                        // Exponent for the agglomeration function. Must be an even integer
                aSQIParameterList.insert( "agglomeration_reference", 1.0, MORIS_REAL_EPS, MORIS_REAL_MAX );    // Reference value for the agglomeration function
                aSQIParameterList.insert( "agglomeration_shift", 0.0, -MORIS_REAL_MAX, MORIS_REAL_MAX );       // Shift value for the agglomeration function
                break;

            case sqi::SQI_Type::INSCRIBED_CIRCLE_SHAPE_DIAMETER:
                aSQIParameterList.set( "SQI_type", sqi::SQI_Type::INSCRIBED_CIRCLE_SHAPE_DIAMETER );
                aSQIParameterList.insert( "number_of_samples", 1, 1, 1000 );                                   // Number of unique minimum inscribed circles to compute for each vertex, taking the average of these values as the shape diameter for the facet. More circles makes the measure less noisy, but makes the hessian more dense
                aSQIParameterList.insert( "cone_angle", 120.0, 0.0, 179.9999999 );                             // Cone angle in degrees for the shape diameter function
                aSQIParameterList.insert( "minimum_relative_chord_length", 0.25, 1e-6, 1.0 );                  // Exclude diameters if the chord between the two vertices is less than this fraction of the diameter, to avoid very high sensitivities for edges that are nearly parallel to the surface
                aSQIParameterList.insert( "agglomeration_exponent", 1.0, 1.0, 1000.0 );                        // Exponent for the agglomeration function. Must be an even integer
                aSQIParameterList.insert( "agglomeration_reference", 1.0, MORIS_REAL_EPS, MORIS_REAL_MAX );    // Reference value for the agglomeration function
                aSQIParameterList.insert( "agglomeration_shift", 0.0, -MORIS_REAL_MAX, MORIS_REAL_MAX );       // Shift value for the agglomeration function
                break;

            case sqi::SQI_Type::SHORTEST_DISTANCE_SHAPE_DIAMETER:
                aSQIParameterList.set( "SQI_type", sqi::SQI_Type::SHORTEST_DISTANCE_SHAPE_DIAMETER );
                aSQIParameterList.insert( "number_of_samples", 1, 1, 1000 );                                   // Number of unique minimal distances to compute for each vertex, taking the average of these values as the shape diameter for the facet. More samples makes the measure less noisy, but makes the hessian more dense
                aSQIParameterList.insert( "cone_angle", 20.0, 0.0, 179.9999999 );                              // Cone angle in degrees for the shape diameter function
                aSQIParameterList.insert( "agglomeration_exponent", 1.0, 1.0, 1000.0 );                        // Exponent for the agglomeration function. Must be an even integer
                aSQIParameterList.insert( "agglomeration_reference", 1.0, MORIS_REAL_EPS, MORIS_REAL_MAX );    // Reference value for the agglomeration function
                aSQIParameterList.insert( "agglomeration_shift", 0.0, -MORIS_REAL_MAX, MORIS_REAL_MAX );       // Shift value for the agglomeration function
                break;
        }
    }

    //------------------------------------------------------------------------------

    // creates a SQI parameter list with default inputs
    inline Parameter_List
    create_SQI_parameter_list( sqi::SQI_Type aSQI )
    {
        Parameter_List tParameterList = create_SQI_parameter_list();
        insert_SQI_parameters( tParameterList, aSQI );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::prm
