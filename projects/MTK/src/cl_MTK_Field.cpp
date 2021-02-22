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

        //------------------------------------------------------------------------------

        Field::~Field()
        {

        }

        //------------------------------------------------------------------------------

        Mesh_Pair * Field::get_mesh_pair()
        {
            MORIS_ERROR( mMeshPair != nullptr, " Field::get_mesh_pair()(), Mesh_Manager not set" );

            return mMeshPair;
        }

        //------------------------------------------------------------------------------

        void Field::set_mesh_pair( Mesh_Pair * aMeshPair)
        {
            this->error_if_locked();

            mMeshPair = aMeshPair;

            this->compute_nodal_values();

            mFieldIsLocked = true;
        }

        //------------------------------------------------------------------------------

        void Field::save_field_to_hdf5(
                const std::string & aFilePath,
                const bool          aCreateNewFile )
        {
            // test if file exists
            std::string tFilePath = make_path_parallel( aFilePath );

            // test if file exists
            std::ifstream tFile( tFilePath );
            bool tFileExists;
            if( tFile )
            {
                tFileExists = true;
            }
            else
            {
                tFileExists = false;
            }

            tFile.close();

            // delete file if it exists and user does not want to keep it
            if( aCreateNewFile && tFileExists )
            {
                std::remove( tFilePath.c_str() );
                tFileExists = false;
            }

            hid_t tFileID;

            if( tFileExists )
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

        void Field::save_node_values_to_hdf5(
                const std::string & aFilePath,
                const bool          aCreateNewFile )
        {
            // test if file exists
            std::string tFilePath = make_path_parallel( aFilePath );

            // test if file exists
            std::ifstream tFile( tFilePath );
            bool tFileExists;
            if( tFile )
            {
                tFileExists = true;
            }
            else
            {
                tFileExists = false;
            }

            tFile.close();

            // delete file if it exists and user does not want to keep it
            if( aCreateNewFile && tFileExists )
            {
                std::remove( tFilePath.c_str() );
                tFileExists = false;
            }

            hid_t tFileID;

            if( tFileExists )
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
                    this->get_nodal_values(),
                    tStatus );

            // close file
            tStatus = close_hdf5_file( tFileID );
        }

        //------------------------------------------------------------------------------

        void Field::load_field_from_hdf5(
                const std::string & aFilePath,
                const uint          aBSplineOrder )
        {
            hid_t tFile    = open_hdf5_file( aFilePath );
            herr_t tStatus = 0;
            Matrix<DDRMat> tMat;
            load_matrix_from_hdf5_file( tFile,
                    this->get_label(),
                    tMat,
                    tStatus );

            this->set_coefficients( tMat );

            tStatus = close_hdf5_file( tFile );
        }

        //------------------------------------------------------------------------------

        void Field::save_node_values_to_binary( const std::string & aFilePath )
        {
            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            save_matrix_to_binary_file( this->get_nodal_values(), tFilePath );
        }

        //------------------------------------------------------------------------------

        void Field::save_bspline_coeffs_to_binary( const std::string & aFilePath )
        {
            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            save_matrix_to_binary_file( this->get_coefficients(), tFilePath );
        }

        //------------------------------------------------------------------------------

        void Field::save_field_to_exodus( const std::string & aFileName )
        {
            mtk::Mesh * tMesh = mMeshPair->mInterpolationMesh;

            moris::mtk::Writer_Exodus tExodusWriter( tMesh );

            tExodusWriter.write_mesh(
                    "./",
                    aFileName,
                    "./",
                    "field_temp");

            // set standard field names
            moris::Cell<std::string> tNodalFieldNames( 1 );

            tNodalFieldNames( 0 ) = "Field";

            // pass nodal field names to writer
            tExodusWriter.set_nodal_fields( tNodalFieldNames );

            tExodusWriter.set_time( 0.0 );

            tExodusWriter.write_nodal_field( tNodalFieldNames( 0 ), mNodalValues );

            tExodusWriter.save_mesh( );
        }

        //------------------------------------------------------------------------------

        void Field::error_if_locked() const
        {
            MORIS_ERROR( !mFieldIsLocked,
                    "Field is locked. You are not allowed to change the mesh as well as the field or coefficient vector");
        }


        //------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
