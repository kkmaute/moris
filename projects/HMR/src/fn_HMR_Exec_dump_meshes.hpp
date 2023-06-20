/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_HMR_Exec_dump_meshes.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_MESHES_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_MESHES_HPP_

#include <string>

#include "cl_HMR.hpp"
#include "cl_HMR_Arguments.hpp"
#include "assert.hpp"
#include "typedefs.hpp"

namespace moris::hmr
{
//--------------------------------------------------------------------------------

    void dump_meshes( const Arguments & aArguments,
                      const Paramfile & aParamfile,
                            HMR       * aHMR )
    {
        // test if an output database path is given
        if( aParamfile.get_output_db_path().size() > 0 )
        {
            aHMR->save_to_hdf5( aParamfile.get_output_db_path(),0 ); //FIXME
        }

        // test if coefficient path is given
        if( aParamfile.get_coefficient_db_path().size() > 0 )
        {
            aHMR->save_coeffs_to_hdf5_file( aParamfile.get_coefficient_db_path(), 0 );
        }

        // loop over all output meshes
        for( uint m=0; m<aParamfile.get_number_of_meshes(); ++m )
        {
            // test if path is given
            if ( aParamfile.get_mesh_path( m ).size() > 0 )
            {
                // get orde rof mesh
                uint tOrder = aParamfile.get_mesh_order( m );

                uint tIndex = 0;

                MORIS_ERROR(false, "HMR::get_mesh_index() this function is not udated yet ");
                // get index of mesh order
                if( tOrder <= 2 )
                {
//                        tIndex = aHMR->get_mesh_index( tOrder, aHMR->get_parameters()->get_lagrange_output_pattern() );
                }
                else
                {
//                        tIndex = aHMR->get_mesh_index( tOrder, aHMR->get_parameters()->get_refined_output_pattern() );
                }

                // dump mesh
                aHMR->save_to_exodus( tIndex,
                                      aParamfile.get_mesh_path( m ),
                                      aArguments.get_timestep() );

                // also save last step
                if( aArguments.get_state() == State::REFINE_MESH && tOrder < 3 )
                {
                    // get path
                    std::string tOrgPath = aParamfile.get_mesh_path( m );

                    // add suffix
                    std::string tPath =
                            tOrgPath.substr(0,tOrgPath.find_last_of(".")) // base path
                            + "_last_step" +
                            tOrgPath.substr( tOrgPath.find_last_of("."), tOrgPath.length() ); // file extension

                    MORIS_ERROR(false, "HMR::get_mesh_index() this function is not udated yet ");
//                        tIndex = aHMR->get_mesh_index( tOrder, aHMR->get_parameters()->get_lagrange_input_pattern() );

                    // dump mesh
                    aHMR->save_to_exodus( tIndex,
                                          tPath,
                                          aArguments.get_timestep() );
                }
            }
        }
    }

//--------------------------------------------------------------------------------
}

#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_MESHES_HPP_ */

