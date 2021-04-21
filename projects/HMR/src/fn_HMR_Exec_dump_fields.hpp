/*
 * fn_HMR_Exec_dump_fields.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_FIELDS_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_FIELDS_HPP_

#include <string>

#include "cl_HMR.hpp"
#include "cl_HMR_Arguments.hpp"
#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "HDF5_Tools.hpp"

namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------

        void dump_fields( const Paramfile & aParams, Cell< std::shared_ptr< Field > > & aFields )
        {
            // before we dump the fields, we must find out if they are all
            // written into the same output database or not
            // we create a cell of strings

            Cell< std::string > tOutputFiles;

            uint tNumberOfFields = aParams.get_number_of_fields();

            // next, we map the target paths with the field IDs
            map< moris_id, std::string > tFieldIdToTarget;
            for( uint f=0; f<tNumberOfFields; ++f )
            {
                tFieldIdToTarget[ aParams.get_field_params( f ).mID ] =  aParams.get_field_params( f ).mTarget;

                tOutputFiles.push_back( aParams.get_field_params( f ).mTarget );
            }

            // next, we make the files unique
            unique( tOutputFiles );

            // remember number of files
            uint tNumberOfFiles = tOutputFiles.size();

            // now, we crate a map with the indices
            map< std::string, uint > tFileMap;

            for( uint k=0; k<tNumberOfFiles; ++k )
            {
                tFileMap[ tOutputFiles( k ) ] = k;
            }

            // create a moris mat with flags
            Matrix< DDUMat > tFileFlags( tNumberOfFiles, 1, 0 );

            // Cell of file IDs
            Cell< hid_t > tFileIDs( tNumberOfFiles, 0 );

            // loop over all fields
            for(  std::shared_ptr< Field > tField : aFields )
            {
                // get path
                std::string tFilePath = tFieldIdToTarget.find( tField->get_id() );

                // error handler
                herr_t tStatus = 0;

                // get index of output file
                uint tIndex = tFileMap.find( tFilePath );

                // test if file has been created already
                if ( tFileFlags( tIndex ) == 0 )
                {
                    // create file
                    tFileIDs( tIndex ) = create_hdf5_file( tFilePath );
                    tFileFlags( tIndex ) = 1;
                }

                save_matrix_to_hdf5_file( tFileIDs( tIndex ),
                                          tField->get_label(),
                                          tField->get_coefficients(),
                                          tStatus );
            }

            // close all open files
            for( hid_t tFileID : tFileIDs )
            {
                close_hdf5_file( tFileID );
            }
        }
//--------------------------------------------------------------------------------
    }
}

#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_DUMP_FIELDS_HPP_ */
