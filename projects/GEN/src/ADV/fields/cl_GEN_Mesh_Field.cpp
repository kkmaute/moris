/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Mesh_Field.cpp
 *
 */

#include "cl_GEN_Mesh_Field.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "HDF5_Tools.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Mesh_Field::Mesh_Field(
                mtk::Mesh*  aMesh,
                std::string aFieldName,
                mtk::EntityRank  aEntityRank )
                 : Field_Discrete_Integration( Matrix< DDRMat >{{}}, aMesh->get_num_nodes())
                 , mMesh(aMesh)
                 , mFieldName(aFieldName)
                 , mEntityRank(aEntityRank)
                 , mUseOwnData(false)
         {
         }

        //--------------------------------------------------------------------------------------------------------------
         Mesh_Field::Mesh_Field(
                mtk::Mesh*  aMesh,
                std::string aFileName,
                std::string aFieldName,
                std::string aFileFormat,
                real        aOffset,
                mtk::EntityRank  aEntityRank )
                : Field_Discrete_Integration( Matrix< DDRMat>{{}}, aMesh == nullptr ? 0 : aMesh->get_num_nodes())
                , mMesh(aMesh)
                , mFieldName(aFieldName)
                , mOffset(aOffset)
                , mEntityRank(aEntityRank)
                , mUseOwnData(true)
        {
            // initialize error flag
            bool tError=true;

            // check for empty file name - use trivial field values
            if ( aFileName.empty() )
            {
                mFieldData.set_size(1,1,mOffset);

                return;
            }

            // read data from exodus file
            if ( aFileFormat == "exodus" )
            {
                // open and query exodus output file (set verbose to true to get basic mesh information)
                moris::mtk::Exodus_IO_Helper tExoIO(aFileName.c_str(),0,false,false);

                // get number of nodes from mesh on file
                uint tNumNodes = tExoIO.get_number_of_nodes();

                // get field index from field name
                uint tFieldIndex = tExoIO.get_field_index_by_name( aFieldName );

                // check that node numbers are consistent
                if ( aMesh != nullptr )
                {
                    MORIS_ERROR( tNumNodes == aMesh->get_num_nodes(),
                            "Mesh_Field::Mesh_Field - inconsistent number of nodes.");
                }

                mFieldData=tExoIO.get_nodal_field_vector(tFieldIndex);

                // set error flag to false
                tError=false;
            }

            // read data from hdf5 file
            if ( aFileFormat == "hdf5" )
            {
                // open hdf5 file
                hid_t tFileID = open_hdf5_file( aFileName );
                herr_t tStatus = 0;

                // load data from file
                load_matrix_from_hdf5_file( tFileID, aFieldName, mFieldData, tStatus);

                // close file
                close_hdf5_file( tFileID );

                // set error flag to false
                tError=false;
            }

            // check the data was read
            MORIS_ERROR( tError == false,
                    "Mesh_Field::Mesh_Field - incorrect file format.");
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Mesh_Field::get_field_value(uint aNodeIndex)
        {
            if ( mUseOwnData )
            {
                if ( aNodeIndex < mFieldData.numel() )
                {
                    return mFieldData(aNodeIndex) + mOffset;
                }
                else
                {
                    return mOffset;
                }
            }
            else
            {
                return mMesh->get_entity_field_value_real_scalar({{moris_index(aNodeIndex)}}, mFieldName, mEntityRank)(0);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>&
        Mesh_Field::get_dfield_dadvs(uint aNodeIndex)
        {
            MORIS_ERROR(false, "get_dfield_dadvs function is not implemented for a mesh field.");
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

