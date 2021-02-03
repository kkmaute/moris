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
                    this->get_node_values(),
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
            load_matrix_from_hdf5_file( tFile,
                    this->get_label(),
                    this->get_coefficients(),
                    tStatus );

            tStatus = close_hdf5_file( tFile );
        }

        //------------------------------------------------------------------------------

        void Field::save_node_values_to_binary( const std::string & aFilePath )
        {
            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            save_matrix_to_binary_file( this->get_node_values(), tFilePath );
        }

        //------------------------------------------------------------------------------

        void Field::save_bspline_coeffs_to_binary( const std::string & aFilePath )
        {
            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            save_matrix_to_binary_file( this->get_coefficients(), tFilePath );
        }

        //------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
