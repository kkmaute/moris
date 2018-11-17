/*
 * fn_HMR_Exec_dump_meshes.hpp
 *
 *  Created on: Nov 13, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_MESHES_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_MESHES_HPP_

#include <string>

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Arguments.hpp"

namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------

        void
        dump_meshes( const Arguments & aArguments, HMR * aHMR )
        {

            // container for file list
            ParameterList tFileList;

            // load the passed files from the parameter list
            load_file_list_from_xml(
                    aArguments.get_parameter_path(), tFileList );

            // get default order for output
            uint tOutputMeshOrder = tFileList.get< sint >( "output_exodus_order" );

            // grab output path from file list
            std::string tOutputDB = tFileList.get< std::string >("output_mesh_database");

            // test if an output database path is given
            if( tOutputDB.size() > 0 )
            {
                aHMR->save_to_hdf5( tOutputDB );
            }

            // grab coefficient path from file list
            std::string tCoeffPath = tFileList.get< std::string >("coefficient_database");

            // test if coefficient path is given
            if( tCoeffPath.size() > 0 )
            {
                aHMR->save_coeffs_to_hdf5_file( tCoeffPath );
            }

            // check state
            switch ( aArguments.get_state() )
            {
                case( State::INITIALIZE_MESH ) :
                {
                    // grab output path from file list
                    std::string tOutputExo = tFileList.get< std::string > ("output_exodus");

                    // test if path is given
                    if( tOutputExo.size() > 0 )
                    {
                        // write file and copy timestep from arguments
                        aHMR->save_to_exodus( tOutputExo, aArguments.get_timestep() , tOutputMeshOrder );
                    }
                    break;
                }
                case( State::REFINE_MESH ) :
                {
                    // grab path to last step exodus
                    std::string tLastStep =  tFileList.get< std::string > ("last_step_exodus");

                    // test if path is given
                    if( tLastStep.size() > 0 )
                    {
                        aHMR->save_last_step_to_exodus( tLastStep, aArguments.get_timestep() , tOutputMeshOrder );
                    }

                    // test if refined mesh is given
                    std::string tRefinedMesh = tFileList.get< std::string > ("refined_exodus");
                    {
                        aHMR->save_to_exodus( tRefinedMesh, aArguments.get_timestep() , tOutputMeshOrder );
                    }

                    break;
                }
                case( State::MAP_FIELDS ) :
                {
                    // grab output path from file list
                    std::string tOutputExo = tFileList.get< std::string > ("output_exodus");

                    // test if path is given
                    if( tOutputExo.size() > 0 )
                    {
                        // write file and copy timestep from arguments
                        aHMR->save_to_exodus( tOutputExo, aArguments.get_timestep() );
                    }
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Invalid state in dump_meshes.");
                }
            }
        }

//--------------------------------------------------------------------------------
    }
}

#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_MESHES_HPP_ */
