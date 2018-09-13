/*
 * cl_HMR_Fields.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_FIELDS_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_FIELDS_HPP_

#include <string>
#include "typedefs.hpp" //COR/src
#include "cl_Mat.hpp"
#include "cl_Cell.hpp"
#include "cl_XML_Parser.hpp"
#include "HMR_Tools.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        /**
         * a struct that defines fields that are defined in the xml file
         */
        struct Field_Parameters
        {
            std::string mLabel;
            std::string mInputValuesPath;
            std::string mInputCoeffsPath;
            std::string mOutputValuesPath;
            std::string mOutputCoeffsPath;
            bool        mRefinementFlag;
            uint        mInterpolationOrder;
        };

//------------------------------------------------------------------------------

        void
        load_field_parameters_from_xml( const std::string        & aSettingsPath,
                                        Cell< Field_Parameters > & aSettings )
        {

            // clean up output
            aSettings.clear();

            // create a parser object for this settings path
            XML_Parser tParser( aSettingsPath );

            // get number of fields from parser
            uint tNumberOfFields = tParser.count_keys_in_subtree( "moris.hmr", "field" );


            // loop over all fields
            for( uint f=0; f<tNumberOfFields; ++f )
            {
                Cell< std::string > tFirst;
                Cell< std::string > tSecond;

                tParser.get_keys_from_subtree( "moris.hmr", "field", f, tFirst, tSecond );
                Field_Parameters tParams;

                // copy key to settings struct
                for( uint k=0; k<tFirst.size(); ++k )
                {
                    std::string tKey = tFirst( k );

                    if ( tKey == "label" )
                    {
                        tParams.mLabel = tSecond( k );
                    }
                    else if ( tKey == "interpolation_order" )
                    {
                        tParams.mInterpolationOrder = stoi( tSecond( k ) );
                    }
                    else if ( tKey == "input_values" )
                    {
                        tParams.mInputValuesPath = parallelize_path( tSecond( k ) );
                    }
                    else if ( tKey == "input_coeffs" )
                    {
                        tParams.mInputCoeffsPath  = parallelize_path( tSecond( k ) );
                    }
                    else if ( tKey == "output_values" )
                    {
                        tParams.mOutputValuesPath = parallelize_path( tSecond( k ) );
                    }
                    else if ( tKey == "output_coeffs" )
                    {
                        tParams.mOutputCoeffsPath = parallelize_path( tSecond( k ) );
                    }
                    else if ( tKey == "refine" )
                    {
                        tParams.mRefinementFlag = string_to_bool( tSecond( k ) );
                    }
                }

                aSettings.push_back( tParams );
            }
        }

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */


#endif /* PROJECTS_HMR_SRC_CL_HMR_FIELDS_HPP_ */
