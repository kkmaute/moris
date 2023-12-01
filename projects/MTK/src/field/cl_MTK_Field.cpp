/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Field.cpp
 *
 */

#include <iostream>
#include <cstdio>

// HD5 c-interface
#include "hdf5.h"

#include "fn_save_matrix_to_binary_file.hpp"

#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Field.hpp"
#include "linalg_typedefs.hpp"

#include "HDF5_Tools.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "fn_dot.hpp"

namespace moris
{
    namespace mtk
    {
        Field::Field(
                Mesh_Pair              aMeshPair,
                uint const &           aNumberOfFields,
                enum Field_Entity_Type aFieldEntityType )
                : mMeshPair( aMeshPair )
                , mNumberOfFields( aNumberOfFields )
                , mFieldEntityType( aFieldEntityType )
        {
            // get interpolation mesh
            mtk::Mesh* tIPmesh = mMeshPair.get_interpolation_mesh();

            if ( mFieldEntityType == Field_Entity_Type::NODAL )
            {
                // size matrix of nodal values
                mValues.set_size( tIPmesh->get_num_nodes(), mNumberOfFields, MORIS_REAL_MIN );
            }
            else if ( mFieldEntityType == Field_Entity_Type::ELEMENTAL )
            {
                mValues.set_size( tIPmesh->get_num_elems(), mNumberOfFields, MORIS_REAL_MIN );
            }
        }

        //------------------------------------------------------------------------------

        Field::Field( enum Field_Entity_Type aFieldEntityType )
                : mMeshPair( nullptr, nullptr )
                , mNumberOfFields( 1 )
                , mFieldEntityType( aFieldEntityType )
        {
        }

        //------------------------------------------------------------------------------

        Field::~Field()
        {
        }

        //------------------------------------------------------------------------------

        enum Field_Entity_Type
        Field::get_field_entity_type()
        {
            return mFieldEntityType;
        }

        //------------------------------------------------------------------------------

        Mesh_Pair
        Field::get_mesh_pair()
        {
            return mMeshPair;
        }

        //------------------------------------------------------------------------------

        const std::string&
        Field::get_mesh_pair_label()
        {
            MORIS_ERROR( not mMeshPairLabel.empty(), "Field::get_mesh_pair_label(), Mesh pair label not set" );
            return mMeshPairLabel;
        }

        //------------------------------------------------------------------------------

        void
        Field::set_mesh_pair( Mesh_Pair aMeshPair )
        {
            // check whether field is unlocked
            this->error_if_locked();

            // set new mesh pair
            mMeshPair = aMeshPair;

            // set update flag to false; nodal values do not need to be updated
            mUpdateNodalValues = true;

            // update coefficient data as underlying discretization may have changed
            this->update_coefficient_data();

            // lock mesh
            mFieldIsLocked = true;
        }

        //------------------------------------------------------------------------------

        void
        Field::update_field()
        {
            // set update flag to false; nodal values do not need to be updated
            mUpdateNodalValues = true;

            // update coefficient data as underlying discretization may have changed
            this->update_coefficient_data();

            // lock mesh
            mFieldIsLocked = true;
        }

        //------------------------------------------------------------------------------

        moris::real
        Field::get_max_value()
        {
            return mValues.max();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field::get_values()
        {
            // check whether nodal values are updated; if not compute them first
            if ( this->nodal_values_need_update() )
            {
                // compute nodal values for current coefficients
                this->compute_nodal_values();

                // set update flag to false
                mUpdateNodalValues = false;
            }

            // return nodal value vector
            return mValues;
        }

        //------------------------------------------------------------------------------

        moris::real
        Field::get_value(
                const uint& aIndex,
                const uint& aFieldIndex )
        {
            // check whether nodal values are updated; if not compute them first
            if ( this->nodal_values_need_update() )
            {
                this->compute_nodal_values();

                // set update flag to false
                mUpdateNodalValues = false;
            }

            return mValues( aIndex, aFieldIndex );
        }

        //------------------------------------------------------------------------------

        void
        Field::get_value(
                Matrix< IndexMat > const & aIndex,
                Matrix< DDRMat >&          aValues,
                Matrix< IndexMat > const & aFieldIndex )
        {
            // check whether nodal values are updated; if not compute them first
            if ( this->nodal_values_need_update() )
            {
                this->compute_nodal_values();

                // set update flag to false
                mUpdateNodalValues = false;
            }

            // get number of fields and nodal values
            uint tNumFields  = aFieldIndex.numel();
            uint tNumIndices = aIndex.numel();

            // set size for requested values
            aValues.resize( tNumIndices, tNumFields );

            // assemble requested values into matrix
            for ( uint Ik = 0; Ik < tNumFields; Ik++ )
            {
                uint tFieldIndex = aFieldIndex( Ik );

                for ( uint Ii = 0; Ii < tNumIndices; Ii++ )
                {
                    aValues( Ii, Ik ) = mValues( aIndex( Ii ), tFieldIndex );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Field::set_values( const Matrix< DDRMat >& aValues )
        {
            // check whether field is unlocked
            this->error_if_locked();

            // check that input vector has proper size
            MORIS_ERROR(
                    aValues.n_rows() == mValues.n_rows() && aValues.n_cols() == mValues.n_cols(),
                    "mtk::Field::set_value_vector - input vector has incorrect size: %zu vs %zu vs %d.\n",
                    aValues.n_rows(),
                    mValues.n_rows(),
                    mMeshPair.get_interpolation_mesh()->get_num_nodes() );

            // set nodal value vector using child implementation
            this->set_nodal_value_vector( aValues );

            // set update flag to false, i.e., nodal values should not be updated
            mUpdateNodalValues = false;

            // lock field
            mFieldIsLocked = true;
        }

        //------------------------------------------------------------------------------

        void
        Field::set_nodal_value_vector( const Matrix< DDRMat >& aNodalValues )
        {
            // FIXME - loop should not be needed if input vector has correct size
            // copy values
            for ( uint tNodeIndex = 0; tNodeIndex < mValues.n_rows(); ++tNodeIndex )
            {
                for ( uint tFieldIndex = 0; tFieldIndex < mValues.n_cols(); ++tFieldIndex )
                {
                    mValues( tNodeIndex, tFieldIndex ) = aNodalValues( tNodeIndex, tFieldIndex );
                }
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field::get_coefficients()
        {
            // call child implementation to populate coefficient vector
            this->get_coefficient_vector();

            // return coefficient vector
            return mCoefficients;
        }

        //------------------------------------------------------------------------------

        void
        Field::set_coefficients( const Matrix< DDRMat >& aCoefficients )
        {
            // check whether field is unlocked
            this->error_if_locked();

            // check whether vector has the correct size
            MORIS_ERROR( (sint)aCoefficients.n_rows() == mNumberOfCoefficients,
                    "mtk::Field::set_coefficients - coefficient vector has incorrect length: %zu vs %d \n",
                    aCoefficients.n_rows(),
                    mNumberOfCoefficients );

            // set coefficient vector using child implementation
            this->set_coefficient_vector( aCoefficients );

            // set update flag to false; nodal values do not need to be updated
            mUpdateNodalValues = true;

            // lock field
            mFieldIsLocked = true;
        }

        //------------------------------------------------------------------------------

        void
        Field::set_coefficient_vector( const Matrix< DDRMat >& aCoefficients )
        {
            mCoefficients = aCoefficients;
        }

        //------------------------------------------------------------------------------

        const Matrix< IdMat >&
        Field::get_coefficient_ids_and_owners()
        {
            // call child implementation to populate coefficient vector
            return this->get_coefficient_id_and_owner_vector();
        }

        //------------------------------------------------------------------------------

        bool
        Field::nodal_values_need_update()
        {
            if ( mFieldEntityType == Field_Entity_Type::NODAL )
            {
                // get interpolation mesh
                mtk::Mesh* tIPmesh = mMeshPair.get_interpolation_mesh();

                // update is required if either update flag is set or number of nodes have changed
                // FIXME: check on changing nodes should not be needed; if mesh is changed this
                //        should be handled explicitly by setting a new mesh pair
                if ( mUpdateNodalValues || tIPmesh->get_num_nodes() != mValues.n_rows() )
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        //------------------------------------------------------------------------------

        void
        Field::save_coefficients_to_hdf5(
                const std::string& aFilePath,
                const bool         aCreateNewFile )
        {
            MORIS_ERROR( mNumberOfCoefficients > -1,
                    "mtk::Field::save_coefficients_to_hdf5 - coefficient vector not set.\n" );

            // test if file exists
            std::string tFilePath = make_path_parallel( aFilePath );

            // test if file exists
            std::ifstream tFile( tFilePath );
            bool          tFileExists = tFile ? true : false;

            tFile.close();

            // delete file if it exists and user does not want to keep it
            if ( aCreateNewFile && tFileExists )
            {
                std::remove( tFilePath.c_str() );
                tFileExists = false;
            }

            hid_t tFileID;

            if ( tFileExists )
            {
                tFileID = open_hdf5_file( aFilePath );
            }
            else
            {
                tFileID = create_hdf5_file( aFilePath );
            }

            herr_t tStatus;

            save_matrix_to_hdf5_file( tFileID,
                    this->get_label(),
                    this->get_coefficients(),
                    tStatus );

            // close file
            tStatus = close_hdf5_file( tFileID );
        }

        //------------------------------------------------------------------------------

        void
        Field::save_nodal_values_to_hdf5(
                const std::string& aFilePath,
                const bool         aCreateNewFile )
        {
            // test if file exists
            std::string tFilePath = make_path_parallel( aFilePath );

            // test if file exists
            std::ifstream tFile( tFilePath );
            bool          tFileExists;
            if ( tFile )
            {
                tFileExists = true;
            }
            else
            {
                tFileExists = false;
            }

            tFile.close();

            // delete file if it exists and user does not want to keep it
            if ( aCreateNewFile && tFileExists )
            {
                std::remove( tFilePath.c_str() );
                tFileExists = false;
            }

            hid_t tFileID;

            if ( tFileExists )
            {
                tFileID = open_hdf5_file( aFilePath );
            }
            else
            {
                tFileID = create_hdf5_file( aFilePath );
            }

            herr_t tStatus;

            save_matrix_to_hdf5_file( tFileID,
                    this->get_label(),
                    this->get_values(),
                    tStatus );

            // close file
            tStatus = close_hdf5_file( tFileID );
        }

        //------------------------------------------------------------------------------

        void
        Field::load_coefficients_from_hdf5(
                const std::string& aFilePath )
        {
            hid_t            tFile   = open_hdf5_file( aFilePath );
            herr_t           tStatus = 0;
            Matrix< DDRMat > tMat;
            load_matrix_from_hdf5_file( tFile,
                    this->get_label(),
                    tMat,
                    tStatus );

            // check that coefficient vector is column vector
            MORIS_ERROR( tMat.n_cols() == 1,
                    "mtk::Field::load_coefficients_from_hdf5 - vector needs to be column vector.\n" );

            // update number of coefficients
            mNumberOfCoefficients = tMat.n_rows();

            this->unlock_field();
            this->set_coefficients( tMat );

            tStatus = close_hdf5_file( tFile );
        }

        //------------------------------------------------------------------------------

        void
        Field::load_nodal_values_from_hdf5(
                const std::string& aFilePath )
        {
            hid_t            tFile   = open_hdf5_file( aFilePath );
            herr_t           tStatus = 0;
            Matrix< DDRMat > tMat;
            load_matrix_from_hdf5_file( tFile,
                    this->get_label(),
                    tMat,
                    tStatus );

            this->unlock_field();
            this->set_values( tMat );

            tStatus = close_hdf5_file( tFile );
        }

        //------------------------------------------------------------------------------

        void
        Field::save_nodal_values_to_binary( const std::string& aFilePath )
        {
            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            save_matrix_to_binary_file( this->get_values(), tFilePath );
        }

        //------------------------------------------------------------------------------

        void
        Field::save_coefficients_to_binary( const std::string& aFilePath )
        {
            MORIS_ERROR( mNumberOfCoefficients > -1,
                    "mtk::Field::save_coefficients_to_hdf5 - coefficient vector not set.\n" );

            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            save_matrix_to_binary_file( this->get_coefficients(), tFilePath );
        }

        //------------------------------------------------------------------------------

        void
        Field::save_field_to_exodus( const std::string& aFileName )
        {
            mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();

            // set mesh
            moris::mtk::Writer_Exodus tExodusWriter( tMesh );

            // set file names
            tExodusWriter.write_mesh(
                    "./",
                    aFileName,
                    "./",
                    "field_temp" );

            // set standard field names
            moris::Cell< std::string > tNodalFieldNames( mNumberOfFields );

            for ( uint tIndex = 0; tIndex < mNumberOfFields; ++tIndex )
            {
                tNodalFieldNames( tIndex ) = "Field-" + std::to_string( tIndex );
            }

            // pass nodal field names to writer
            tExodusWriter.set_nodal_fields( tNodalFieldNames );

            // set time
            tExodusWriter.set_time( 0.0 );

            // write all fields
            for ( uint tIndex = 0; tIndex < mNumberOfFields; ++tIndex )
            {
                tExodusWriter.write_nodal_field(
                        tNodalFieldNames( tIndex ),
                        this->get_values().get_column( tIndex ) );
            }

            // finalize and write mesh
            tExodusWriter.save_mesh();
        }

        //------------------------------------------------------------------------------

        void
        Field::load_field_from_exodus(
                const std::string&      aFileName,
                const moris_index       aTimeIndex,
                const Matrix< DDUMat >& aFiledIndices )
        {
            // open and query exodus file
            Exodus_IO_Helper tExoIO( aFileName, 0, false, false );

            // check that number of requested field indices is equal to number of fields in nodal data
            MORIS_ERROR( aFiledIndices.numel() == mNumberOfFields,
                    "Field::load_field_from_exodus - number of requested fields incorrect.\n" );

            // get number of nodes
            uint tNumNodes = tExoIO.get_number_of_nodes();

            // check that number of nodes match with the field nodes
            MORIS_ERROR( tNumNodes == mValues.n_rows(),
                    "Field::load_field_from_exodus - number of nodes in exodus file incorrect.\n" );

            // allocate temporary memory to store nodal values
            Matrix< DDRMat > tNodalValues( tNumNodes, mNumberOfFields );

            // get IP Mesh to determine nodal value mapping
            mtk::Mesh* tIPMesh = mMeshPair.get_interpolation_mesh();

            // loop over all requested field indices
            for ( uint tFieldIndex = 0; tFieldIndex < mNumberOfFields; ++tFieldIndex )
            {
                // get exodus field index
                moris_index tExoFieldIndex = aFiledIndices( tFieldIndex );

                // get exodus field for current field index
                const Matrix< DDRMat >& tExodusData = tExoIO.get_nodal_field_vector( tExoFieldIndex, aTimeIndex );

                // assign nodal values based on Field index to node id mapping and from
                for ( uint tNodeIndex = 0; tNodeIndex < tNumNodes; ++tNodeIndex )
                {
                    // getting mtk mesh id for this vertex
                    moris_id tIPVertexID = tIPMesh->get_mtk_vertex( tNodeIndex ).get_id();

                    // assign nodal value based on node IDs
                    tNodalValues( tNodeIndex, tFieldIndex ) =
                            tExodusData( tExoIO.get_node_index_by_Id( tIPVertexID ) );
                }
            }

            this->unlock_field();
            this->set_values( tNodalValues );
        }

        //------------------------------------------------------------------------------

        void
        Field::error_if_locked() const
        {
            MORIS_ERROR( !mFieldIsLocked,
                    "Field is locked. You are not allowed to change the mesh as well as the field or coefficient vector" );
        }

        //------------------------------------------------------------------------------
    }    // namespace mtk
} /* namespace moris */
