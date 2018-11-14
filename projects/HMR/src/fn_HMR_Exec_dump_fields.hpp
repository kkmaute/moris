/*
 * fn_HMR_Exec_dump_fields.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_FIELDS_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_FIELDS_HPP_

#include <string>

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Arguments.hpp"
#include "HDF5_Tools.hpp"

namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------

        void
        dump_fields(
                ParameterList aFileList,
                Cell< std::shared_ptr< Field > > aFields )
        {

            // get path to output
            std::string tPath =  aFileList.get< std::string > ( "output_field_database" );

            if( tPath.size() > 0 )
            {
                // create output file
                hid_t tFileID = create_hdf5_file( tPath );

                // create hdf5 error
                herr_t tStatus = 0;

                for( std::shared_ptr< Field > tField : aFields )
                {
                    save_matrix_to_hdf5_file(
                            tFileID,
                            tField->get_label(),
                            tField->get_coefficients(),
                            tStatus );
                }

                // close file
                close_hdf5_file( tFileID );
            }


        }
//--------------------------------------------------------------------------------
    }
}


#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_FIELDS_HPP_ */
