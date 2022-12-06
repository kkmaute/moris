/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_MORIS_GENERAL_Parameters.hpp
 *
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_MORIS_GENERAL_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_MORIS_GENERAL_PARAMETERS_HPP_

#include "cl_Param_List.hpp"

namespace moris
{
    namespace prm
    {

        //------------------------------------------------------------------------------

        // creates a parameter list with default inputs
        inline ParameterList
        create_moris_general_parameter_list()
        {
            ParameterList tParameterList;

            return tParameterList;
        }

        //------------------------------------------------------------------------------

        inline void
        create_refinement_parameterlist( ParameterList& aParameterlist )
        {
            aParameterlist.insert( "field_names", "" );
            aParameterlist.insert( "levels_of_refinement", "" );
            aParameterlist.insert( "refinement_pattern", "" );

            aParameterlist.insert( "refinement_function_name", "" );

            aParameterlist.insert( "remeshing_copy_old_pattern_to_pattern", "" );
        }

        //------------------------------------------------------------------------------

        inline void
        create_remeshing_parameterlist( ParameterList& aParameterlist )
        {
            // Remeshing mode. Options are "ab_initio", "former"
            aParameterlist.insert( "mode", "" );
            aParameterlist.insert( "remeshing_refinement_pattern", "" );
            aParameterlist.insert( "refinement_function_name", "" );

            aParameterlist.insert( "remeshing_frequency", MORIS_SINT_MAX );

            // mode "ab_initio"
            aParameterlist.insert( "remeshing_field_names", "" );
            aParameterlist.insert( "remeshing_levels_of_refinement", "" );
            aParameterlist.insert( "remeshing_refinement_pattern", "" );

            // mode "based_on_previous"
            aParameterlist.insert( "remeshing_maximum_refinement_level", "" );
            aParameterlist.insert( "remeshing_minimum_refinement_level", "" );

            aParameterlist.insert( "remeshing_copy_old_pattern_to_pattern", "" );

            // minimum refinement level per pattern and level
            // input: pattern, level, iter, level, iter,level, iter, ... ; pattern, level, iter, level, iter,level, iter, ...
            aParameterlist.insert( "minimum_refinement_level", "" );

            aParameterlist.insert( "output_meshes", false );

            // modde "previous"
        }

        //------------------------------------------------------------------------------

        inline void
        create_mapping_parameterlist( ParameterList& aParameterlist )
        {
            aParameterlist.insert( "adv_field", "" );
            aParameterlist.insert( "dof_type", "" );
            aParameterlist.insert( "reinitialization_frequency", 1 );
            aParameterlist.insert( "output_mesh_file", "" );
            aParameterlist.insert( "time_offset", 0.0 );
        }

        //------------------------------------------------------------------------------

    }    // namespace prm
}    // namespace moris

#endif /* PROJECTS_PRM_SRC_FN_PRM_MORIS_GENERAL_PARAMETERS_HPP_ */
