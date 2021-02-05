
#include "cl_MTK_Field_Proxy.hpp"

#include "fn_dot.hpp"

#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {

        //--------------------------------------------------------------------------------------------------------------

        Field_Proxy::Field_Proxy(
                std::shared_ptr<mtk::Mesh_Manager> aMeshManager,
                uint                               aMeshIndex,
                uint                               aDiscretizationMeshIndex )
                : Field(aDiscretizationMeshIndex)
                , mMeshManager(aMeshManager)
                , mMeshIndex(aMeshIndex)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Field_Proxy::~Field_Proxy()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field_Proxy::get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return mNodalValues(aNodeIndex);
        }

        //------------------------------------------------------------------------------

        std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > Field_Proxy::get_mesh_pair()
        {
            return std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> >( mMeshIndex, mMeshManager );
        }

        //--------------------------------------------------------------------------------------------------------------

        //        void Field::save_field_to_hdf5(
        //                const std::string & aFilePath,
        //                const bool          aCreateNewFile )
        //        {
        //            // test if file exists
        //            std::string tFilePath = make_path_parallel( aFilePath );
        //
        //            // test if file exists
        //            std::ifstream tFile( tFilePath );
        //            bool tFileExists;
        //            if( tFile )
        //            {
        //                tFileExists = true;
        //            }
        //            else
        //            {
        //                tFileExists = false;
        //            }
        //
        //            tFile.close();
        //
        //            // delete file if it exists and user does not want to keep it
        //            if( aCreateNewFile && tFileExists )
        //            {
        //                std::remove( tFilePath.c_str() );
        //                tFileExists = false;
        //            }
        //
        //            hid_t tFileID;
        //
        //            if( tFileExists )
        //            {
        //                tFileID = open_hdf5_file( aFilePath );
        //            }
        //            else
        //            {
        //                tFileID = create_hdf5_file( aFilePath );
        //            }
        //
        //            herr_t tStatus;
        //
        //            save_matrix_to_hdf5_file( tFileID,
        //                    this->get_name(),
        //                    this->get_coefficients(),
        //                    tStatus );
        //
        //            // close file
        //            tStatus = close_hdf5_file( tFileID );
        //        }
        //
        //        //------------------------------------------------------------------------------
        //
        //        void Field::save_node_values_to_hdf5(
        //                const std::string & aFilePath,
        //                const bool          aCreateNewFile )
        //        {
        //            // test if file exists
        //            std::string tFilePath = make_path_parallel( aFilePath );
        //
        //            // test if file exists
        //            std::ifstream tFile( tFilePath );
        //            bool tFileExists;
        //            if( tFile )
        //            {
        //                tFileExists = true;
        //            }
        //            else
        //            {
        //                tFileExists = false;
        //            }
        //
        //            tFile.close();
        //
        //            // delete file if it exists and user does not want to keep it
        //            if( aCreateNewFile && tFileExists )
        //            {
        //                std::remove( tFilePath.c_str() );
        //                tFileExists = false;
        //            }
        //
        //            hid_t tFileID;
        //
        //            if( tFileExists )
        //            {
        //                tFileID = open_hdf5_file( aFilePath );
        //            }
        //            else
        //            {
        //                tFileID = create_hdf5_file( aFilePath );
        //            }
        //
        //            herr_t tStatus;
        //
        //            // FIXME
        //            save_matrix_to_hdf5_file( tFileID,
        //                    this->get_name(),
        //                    this->get_node_values(),
        //                    tStatus );
        //
        //            // close file
        //            tStatus = close_hdf5_file( tFileID );
        //        }
        //
        //        //------------------------------------------------------------------------------
        //
        //        void Field::load_field_from_hdf5(
        //                const std::string & aFilePath,
        //                const uint          aBSplineOrder )
        //        {
        //            hid_t tFile    = open_hdf5_file( aFilePath );
        //            herr_t tStatus = 0;
        ////            load_matrix_from_hdf5_file( tFile,
        ////                    this->get_name(),
        ////                    this->get_coefficients(),
        ////                    tStatus );
        //
        //            tStatus = close_hdf5_file( tFile );
        //        }
        //
        //        //------------------------------------------------------------------------------
        //
        //        void Field::save_node_values_to_binary( const std::string & aFilePath )
        //        {
        //            // make path parallel
        //            std::string tFilePath = parallelize_path( aFilePath );
        //
        //            //FIXME
        ////            save_matrix_to_binary_file( this->get_node_values(), tFilePath );
        //        }
        //
        //        //------------------------------------------------------------------------------
        //
        //        void Field::save_bspline_coeffs_to_binary( const std::string & aFilePath )
        //        {
        //            // make path parallel
        //            std::string tFilePath = parallelize_path( aFilePath );
        //
        ////            save_matrix_to_binary_file( this->get_coefficients(), tFilePath );
        //        }
        //
        //        //------------------------------------------------------------------------------

    }
}

