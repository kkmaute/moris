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

            ParameterList tParams;

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
                    // this functionality is not supported yet
                    // else if ( tKey == "interpolation_order" )
                    //{
                    //    aSettings( f ).set( "interpolation_order", (sint) stoi( tSecond( k ) ) );
                    //}
                    else if ( tKey == "input_values" )
                    {
                        aSettings( f ).set( "input_values", parallelize_path( tSecond( k ) ) );

                    }
                    else if ( tKey == "input_coeffs" )
                    {
                        aSettings( f ).set( "input_coeffs",  parallelize_path( tSecond( k ) ) );
                    }
                    else if ( tKey == "output_values" )
                    {
                        aSettings( f ).set( "output_values", parallelize_path( tSecond( k ) ) );
                    }
                    else if ( tKey == "output_coeffs" )
                    {
                        aSettings( f ).set( "output_coeffs", parallelize_path( tSecond( k ) ) );
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
