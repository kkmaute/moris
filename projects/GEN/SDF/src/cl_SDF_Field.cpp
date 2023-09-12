/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Field.cpp
 *
 */

#include "HDF5_Tools.hpp"
#include "cl_SDF_Field.hpp"
#include "SDF_Tools.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        Field::Field(
                const std::string & aLabel,
                Mesh              & aMesh,
                Matrix< DDRMat >  & aData ) :
                        mLabel( aLabel ),
                        mMesh( aMesh ),
                        mData( aData )
        {

        }

//-------------------------------------------------------------------------------

        void
        Field::save_field_to_hdf5( const std::string & aFilePath )
        {
            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            // Create a new file using default properties
            hid_t tFileID = H5Fcreate(
                    tFilePath.c_str(),
                    H5F_ACC_TRUNC,
                    H5P_DEFAULT,
                    H5P_DEFAULT);

            // error handler
            herr_t tStatus;

            // save label of this field
            save_string_to_hdf5_file(
                    tFileID,
                    "Label",
                    mLabel,
                    tStatus );

            // save par rank
            save_scalar_to_hdf5_file(
                    tFileID,
                    "ParRank",
                    par_rank(),
                    tStatus );

            // save node IDs
            save_matrix_to_hdf5_file(
                    tFileID,
                    "NodeIDs",
                    mMesh.get_node_ids(),
                    tStatus );

            // save values
            save_matrix_to_hdf5_file(
                    tFileID,
                    "NodeValues",
                    mData,
                    tStatus );

            // save order to file
            /*save_scalar_to_hdf5_file(
                                tFileID,
                                "LagrangeOrder",
                                mMesh.get_order(),
                                tStatus ); */

            // close file
            tStatus = H5Fclose( tFileID );
        }

//-------------------------------------------------------------------------------
    }
}

