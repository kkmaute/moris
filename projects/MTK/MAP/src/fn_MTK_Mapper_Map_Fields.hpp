/*
 * fn_MTK_Mapper_Map_Fields.hpp
 *
 *  Created on: Nov 5, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_MAP_SRC_FN_MTK_MAPPER_MAP_FIELDS_HPP_
#define PROJECTS_MTK_MAP_SRC_FN_MTK_MAPPER_MAP_FIELDS_HPP_

#include <string>

#include "cl_MTK_Mapper_Arguments.hpp"
#include "cl_MTK_Mapper_Parameters.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_HMR.hpp"

namespace moris
{
    namespace mapper
    {
// -----------------------------------------------------------------------------

        void
        map_fields( const Arguments & aArguments )
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 1: Process input parameters
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // load parameters from XML file
            ParameterList tParameters;
            load_parameters_from_xml( aArguments.get_parameter_path(),
                    tParameters );

            // get input and output path for HMR

            std::string tHmrInputPath = tParameters.get< std::string >( "input_hmr_db" );
            std::string tHmrOutputPath = tParameters.get< std::string >( "output_hmr_db" );

            // pointer to working mesh
            mtk::Mesh * tMesh = nullptr;

            // run in HMR mode
            if( tHmrInputPath.size() > 0 )
            {
                // create HMR object
                hmr::HMR tHMR( tHmrInputPath );

                // load output mesh if pattern is given
                if( tHmrOutputPath.size() == 0 )
                {
                    tHMR.load_output_pattern_from_path( tHmrOutputPath );
                }

                // finalize mesh
                tHMR.finalize();
            }
            else
            {
                MORIS_ERROR( false, "The mapper needs an HMR Database to run" );
            }





        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */



#endif /* PROJECTS_MTK_MAP_SRC_FN_MTK_MAPPER_MAP_FIELDS_HPP_ */
