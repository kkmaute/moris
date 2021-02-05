#include "cl_MTK_BSpline_Field.hpp"
#include "fn_dot.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::BSpline_Field(
                std::shared_ptr<mtk::Mesh_Manager> aMeshManager,
                uint                               aMeshIndex,
                uint                               aDiscretizationMeshIndex )
                : Field(aDiscretizationMeshIndex)
                , mMeshManager(aMeshManager)
                , mMeshIndex(aMeshIndex)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::~BSpline_Field()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void BSpline_Field::set_coefficients(Matrix<DDRMat> aCoefficients)
        {
            mCoefficients = aCoefficients;
        }

        //--------------------------------------------------------------------------------------------------------------

        void BSpline_Field::transfer_coefficients(const BSpline_Field& aField)
        {
            mCoefficients = std::move(aField.mCoefficients);
        }

        //--------------------------------------------------------------------------------------------------------------

        real BSpline_Field::get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
//            // Get mesh
//            Interpolation_Mesh* tInterpolationMesh =
//                    mMeshManager->get_interpolation_mesh( mMeshIndex );
//
//            // Get T-matrix
//            Matrix<IndexMat> tBSplineIndices = mMesh->get_bspline_inds_of_node_loc_ind(aNodeIndex, this->get_discretization_mesh_index());
//            Matrix<DDRMat> tMatrix = mMesh->get_t_matrix_of_node_loc_ind(tNodeIndex, this->get_discretization_mesh_index());
//
//            // Compute field value
//            real tValue = 0.0;
//            for (uint tBSpline = 0; tBSpline < tBSplineIndices.length(); tBSpline++)
//            {
//                tValue += tMatrix(tBSpline) * mCoefficients(tBSplineIndices(tBSpline));
//            }
            return mNodalValues(aNodeIndex);
        }

        //------------------------------------------------------------------------------

        std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > BSpline_Field::get_mesh_pair()
        {
            return std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> >( mMeshIndex, mMeshManager );
        }

        //--------------------------------------------------------------------------------------------------------------

        void BSpline_Field::evaluate_nodal_values()
        {
            Interpolation_Mesh* tInterpolationMesh =
                    mMeshManager->get_interpolation_mesh( mMeshIndex );

            mNodalValues.set_size(tInterpolationMesh->get_num_nodes(), 1);

            for( uint Ik = 0; Ik < tInterpolationMesh->get_num_nodes(); ++Ik )
            {
                // get pointer to node
                auto tNode = &tInterpolationMesh->get_mtk_vertex( Ik );

                // get PDOFs from node
                auto tBSplines = tNode->
                        get_interpolation( this->get_discretization_mesh_index() )->
                        get_coefficients();

                // get T-Matrix
                const Matrix< DDRMat > & tTMatrix = *tNode->
                        get_interpolation( this->get_discretization_mesh_index() )->
                        get_weights();

                // get number of coefficients
                uint tNumberOfCoeffs = tTMatrix.length();

                MORIS_ASSERT( tNumberOfCoeffs > 0, "No coefficients defined for node" ) ;

                // fill coeffs vector
                Matrix< DDRMat > tCoeffs( tNumberOfCoeffs, 1 );
                for( uint Ii = 0; Ii < tNumberOfCoeffs; ++Ii )
                {
                    tCoeffs( Ii ) = mCoefficients( tBSplines( Ii )->get_index() );
                }

                // write value into solution
                mNodalValues( Ik ) = moris::dot( tTMatrix, tCoeffs );
            }
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

