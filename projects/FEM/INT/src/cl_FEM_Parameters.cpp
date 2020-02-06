/*
 * cl_FEM_Parameters.cpp
 *
 *  Created on: Feb 6, 2020
 *      Author: noel
 */

#include "cl_FEM_Parameters.hpp" //FEM/INT/src

#include "assert.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_unique.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        // creates a property parameter list with default inputs
        ParameterList create_property_parameter_list()
        {
            ParameterList tParameterList;

//            tParameterList.insert( "number_of_elements_per_dimension", std::string( "2, 2" ) );
//            tParameterList.insert( "domain_dimensions", std::string( "1, 1" ) );
//            tParameterList.insert( "domain_offset", std::string( "0, 0 ") );
//            tParameterList.insert( "domain_sidesets", std::string( "" ) );
//            tParameterList.insert( "lagrange_output_meshes", std::string( "" ) );
//
//            tParameterList.insert( "lagrange_input_meshes", std::string( "" ) );
//
//            tParameterList.insert( "refinement_buffer", 0 );
//            tParameterList.insert( "staircase_buffer", 0 );
//
//            tParameterList.insert( "lagrange_orders", std::string( "1" ) );
//            tParameterList.insert( "lagrange_pattern", std::string( "0" ) );
//
//            tParameterList.insert( "bspline_orders", std::string( "1" ) );
//            tParameterList.insert( "bspline_pattern", std::string( "0" ) );
//
//            tParameterList.insert( "union_pattern", 6 );
//            tParameterList.insert( "working_pattern", 7 );
//
//            tParameterList.insert( "lagrange_to_bspline", std::string( "0" ) );
//
//            tParameterList.insert( "severity_level", 1 );
//            tParameterList.insert( "truncate_bsplines", 1 );
//
//            tParameterList.insert( "use_multigrid", 0 );
//            tParameterList.insert( "use_refinement_interrelation", 0 );
//            tParameterList.insert( "renumber_lagrange_nodes", 0 );
//            tParameterList.insert( "use_number_aura", 0 );
//
//            tParameterList.insert( "initial_refinement", 0 );
//            tParameterList.insert( "additional_lagrange_refinement", 0 );
//
//            tParameterList.insert( "max_refinement_level", -1 );

            return tParameterList;
        }

//------------------------------------------------------------------------------
        // creates a constitutive model parameter list with default inputs
        ParameterList create_constitutive_model_parameter_list()
        {
            ParameterList tParameterList;

            return tParameterList;
        }

//------------------------------------------------------------------------------
        // creates a stabilization parameter parameter list with default inputs
        ParameterList create_stabilization_parameter_parameter_list()
        {
            ParameterList tParameterList;

            return tParameterList;
        }

//------------------------------------------------------------------------------
        // creates an IWG parameter list with default inputs
        ParameterList create_IWG_parameter_list()
        {
            ParameterList tParameterList;

            return tParameterList;
        }

//------------------------------------------------------------------------------
        // creates an IQI parameter list with default inputs
        ParameterList create_IQI_parameter_list()
        {
            ParameterList tParameterList;

            return tParameterList;
        }

//------------------------------------------------------------------------------
    }
}


