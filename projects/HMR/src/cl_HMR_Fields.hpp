/*
 * cl_HMR_Fields.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_FIELDS_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_FIELDS_HPP_

#include <string>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "cl_XML_Parser.hpp"
#include "HMR_Globals.hpp"
#include "HMR_Tools.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        ParameterList
        create_fields_parameter_list()
        {
            ParameterList aParams;
            aParams.insert( "label","untitled" );
            aParams.insert( "input_hdf5","" );
            aParams.insert( "input_values","" );
            aParams.insert( "input_coeffs","" );
            aParams.insert( "output_hdf5","" );
            aParams.insert( "output_values","" );
            aParams.insert( "output_coeffs","" );
            aParams.insert( "refine",  (sint) 0 );
            aParams.insert( "lagrange_order",  (sint) 1 );
            aParams.insert( "bspline_order",  (sint) 1 );
            aParams.insert( "min_volume_refinement_level",  ( sint ) 0 );
            aParams.insert( "max_volume_refinement_level",  ( sint ) gMaxNumberOfLevels );
            aParams.insert( "min_surface_refinement_level", ( sint ) 0 );
            aParams.insert( "max_surface_refinement_level", ( sint ) gMaxNumberOfLevels );
            return aParams;
        }

//------------------------------------------------------------------------------
        /**
         * loads the field parameters from an XML file into a Cell.
         * This function is called by the HMR executable
         */
        void
        load_field_parameters_from_xml( const std::string        & aSettingsPath,
                                        Cell< ParameterList >    & aSettings )
        {

            // clean up output
            aSettings.clear();

            // create a parser object for this settings path
            XML_Parser tParser( aSettingsPath );

            // get number of fields from parser
            uint tNumberOfFields = tParser.count_keys_in_subtree( "moris.hmr", "field" );

            // create an empty parameter list
            ParameterList tParams =  create_fields_parameter_list();

            aSettings.resize( tNumberOfFields, tParams );

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
                        aSettings( f ).set( "label", std::string(  tSecond( k ) ) );
                    }
                    else if ( tKey == "input_hdf5" )
                    {
                        aSettings( f ).set( "input_hdf5", tSecond( k ) );

                    }
                    else if ( tKey == "input_values" )
                    {
                        aSettings( f ).set( "input_values", tSecond( k ) );

                    }
                    else if ( tKey == "input_coeffs" )
                    {
                        aSettings( f ).set( "input_coeffs", tSecond( k ) );
                    }
                    else if ( tKey == "output_hdf5" )
                    {
                        aSettings( f ).set( "output_hdf5", tSecond( k ) );
                    }
                    else if ( tKey == "output_values" )
                    {
                        aSettings( f ).set( "output_values", tSecond( k ) );
                    }
                    else if ( tKey == "lagrange_order" )
                    {
                        aSettings( f ).set( "lagrange_order", tSecond( k ) );
                    }
                    else if ( tKey == "bspline_order" )
                    {
                        aSettings( f ).set( "bspline_order", tSecond( k ) );
                    }
                    else if ( tKey == "output_coeffs" )
                    {
                        aSettings( f ).set( "output_coeffs",  tSecond( k ) );
                    }
                    else if ( tKey == "min_volume_refinement_level" )
                    {
                        aSettings( f ).set( "min_volume_refinement_level", (sint) std::stoi( tSecond( k ) ) );
                    }
                    else if ( tKey == "max_volume_refinement_level" )
                    {
                        aSettings( f ).set( "max_volume_refinement_level", (sint) std::stoi( tSecond( k ) ) );
                    }
                    else if ( tKey == "min_surface_refinement_level" )
                    {
                        aSettings( f ).set( "min_surface_refinement_level", (sint) std::stoi( tSecond( k ) ) );
                    }
                    else if ( tKey == "max_surface_refinement_level" )
                    {
                        aSettings( f ).set( "max_surface_refinement_level", (sint) std::stoi( tSecond( k ) ) );
                    }
                    else if ( tKey == "refine" )
                    {
                        aSettings( f ).set( "refine", (sint) string_to_bool( tSecond( k ) ) );
                    }
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */


#endif /* PROJECTS_HMR_SRC_CL_HMR_FIELDS_HPP_ */
