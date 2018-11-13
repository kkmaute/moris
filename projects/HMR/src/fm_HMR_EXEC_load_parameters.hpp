/*
 * fm_HMR_EXEC_load_parameters.hpp
 *
 *  Created on: Nov 12, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_FM_HMR_EXEC_LOAD_PARAMETERS_HPP_
#define PROJECTS_HMR_SRC_FM_HMR_EXEC_LOAD_PARAMETERS_HPP_

#include <string>

#include "typedefs.hpp" //COR/src
#include "HMR_Globals.hpp"

// -----------------------------------------------------------------------------

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        void
        load_file_list_from_xml(
                const std::string & aFilePath,
                ParameterList     & aFileList )
        {
            // list of filenames
            Cell< std::string > tKeys;

            tKeys.push_back("input_mesh_database");
            tKeys.push_back("output_mesh_database");
            tKeys.push_back("input_exodus");
            tKeys.push_back("last_step_exodus");
            tKeys.push_back("refined_step_exodus");
            tKeys.push_back("output_exodus");
            tKeys.push_back("input_field_database");
            tKeys.push_back("output_field_database");
            tKeys.pysh_back("coefficient_database");

            // create parameter list with empty entries
            for( std::string tKey : tKeys )
            {
                aFileList.insert( tKey.c_str(), std::string("") );
            }

            // create temporary Parser object
            XML_Parser tParser( aFilePath );
            Cell< std::string > tFirst;
            Cell< std::string > tSecond;

            // fill parser with values
            tParser.get_keys_from_subtree( "moris.hmr", "files", 0, tFirst, tSecond );

            // loop over all entries
            for( uint k=0; k<tFirst.size(); ++k )
            {
                // compare first with second
                for( std::string tKey : tKeys )
                {
                    if ( tKey == tFirst ( k ) )
                    {
                        // copy value into key list
                        aFileList.set( tKey.c_str(), tSecond( k )  );
                        break;
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        load_refinement_parameters_from_xml(
                const std::string & aFilePath,
                ParameterList     & aRefParams )
        {
            // list of filenames
            Cell< std::string > tKeys;

            tKeys.push_back("library");
            tKeys.push_back("function");

            // create parameter list with empty entries
            for( std::string tKey : tKeys )
            {
                aRefParams.insert( tKey.c_str(), std::string("") );
            }


            // create temporary Parser object
            XML_Parser tParser( aFilePath );
            Cell< std::string > tFirst;
            Cell< std::string > tSecond;
            tParser.get_keys_from_subtree( "moris.hmr", "refinement", 0, tFirst, tSecond );


            for( std::string tKey : tFirst )
            {
                case( "library" ) :
                case( "function") :
                {
                    aRefParams.set( tKey.c_str(), tSecond( k ) );
                    break;
                }
                default :
                {
                    std::string tLabel = tSecond( k );

                }
            }
        }
// -----------------------------------------------------------------------------
    }
}



#endif /* PROJECTS_HMR_SRC_FM_HMR_EXEC_LOAD_PARAMETERS_HPP_ */
