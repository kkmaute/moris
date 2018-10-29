
// HD5 c-interface
#include "hdf5.h"

#include "fn_save_matrix_to_binary_file.hpp"

#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "HMR_Tools.hpp"
#include "HDF5_Tools.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        Field::Field(
                const std::string             & aLabel,
                std::shared_ptr< mtk::Mesh >    aMesh,
                const uint                    & aBSplineOrder,
                std::shared_ptr< Database >     aDatabase,
                Lagrange_Mesh_Base *            aLagrangeMesh ) :
                        mtk::Field( aLabel, aMesh ),
                        mBSplineOrder( aBSplineOrder ),
                        mDatabase( aDatabase ),
                        mLagrangeMesh( aLagrangeMesh ),
                        mFieldIndex( aLagrangeMesh->create_field_data( aLabel ) )

        {

        }

//------------------------------------------------------------------------------

        Field::~Field()
        {
        }

//------------------------------------------------------------------------------

        // parameter copied from input settings
        void
        Field::set_min_surface_level( const uint & aLevel )
        {
            mMinSurfaceLevel = aLevel;
        }

//------------------------------------------------------------------------------

        void
        Field::set_min_volume_level( const uint & aLevel )
        {
            mMinVolumeLevel = aLevel;
        }

//------------------------------------------------------------------------------

        void
        Field::set_max_surface_level( const uint & aLevel )
        {
            mMaxSurfaceLevel = aLevel;
        }

//------------------------------------------------------------------------------

        void
        Field::set_max_volume_level( const uint & aLevel )
        {
            mMaxVolumeLevel = aLevel;
        }

//------------------------------------------------------------------------------

        uint
        Field::get_min_surface_level() const
        {
            return mMinSurfaceLevel;
        }

//------------------------------------------------------------------------------

        // parameter copied from input settings
        uint
        Field::get_min_volume_level() const
        {
            return mMinVolumeLevel;
        }

//------------------------------------------------------------------------------

        // parameter copied from input settings
        uint
        Field::get_max_surface_level() const
        {
            return mMaxSurfaceLevel;
        }

//------------------------------------------------------------------------------

        // parameter copied from input settings
        uint
        Field::get_max_volume_level() const
        {
            return mMaxVolumeLevel;
        }

//------------------------------------------------------------------------------

        const std::string &
        Field::get_label() const
        {
            return mLagrangeMesh->get_field_label( mFieldIndex );

        }

//------------------------------------------------------------------------------

        void
        Field::set_label( const std::string & aLabel )
        {
            mLagrangeMesh->set_field_label( mFieldIndex, aLabel );
        }

//------------------------------------------------------------------------------


        Matrix< DDRMat > &
        Field::get_node_values()
        {
            return mLagrangeMesh->get_field_data( mFieldIndex );
        }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        Field::get_node_values() const
        {
            return mLagrangeMesh->get_field_data( mFieldIndex );
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > &
        Field::get_coefficients()
        {
            return mLagrangeMesh->get_field_coeffs( mFieldIndex );
        }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        Field::get_coefficients() const
        {
            return mLagrangeMesh->get_field_coeffs( mFieldIndex );
        }

//------------------------------------------------------------------------------

        void
        Field::change_mesh( Lagrange_Mesh_Base * aMesh, const uint aFieldIndex )
        {
            mLagrangeMesh = aMesh;
            mFieldIndex   = aFieldIndex;
        }

//------------------------------------------------------------------------------

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
                    this->get_label(),
                    tStatus );

            // save par rank
            save_scalar_to_hdf5_file(
                    tFileID,
                    "ParRank",
                    par_rank(),
                    tStatus );

            // save interpolation order of Lagrange mesh
            save_scalar_to_hdf5_file(
                    tFileID,
                    "LagrangeOrder",
                    mLagrangeMesh->get_order(),
                    tStatus );

            // save interpolation order of B-Spline mesh
            save_scalar_to_hdf5_file(
                    tFileID,
                    "BSplineOrder",
                    mBSplineOrder,
                    tStatus );

           // save node values
           save_matrix_to_hdf5_file(
                    tFileID,
                    "NodeValues",
                    this->get_node_values(),
                    tStatus );

           // create ID matrix
           // get number of nodes
           uint tNumberOfNodes =
                   this->get_mesh()->get_number_of_nodes_on_proc();

           Matrix< IdMat > tIDs( tNumberOfNodes, 1 );
           for( uint k=0; k<tNumberOfNodes; ++k )
           {
               tIDs( k ) = this->get_mesh()->get_node_by_index( k )->get_id();
           }

           // save IDs
           save_matrix_to_hdf5_file(
                               tFileID,
                               "NodeIDs",
                               tIDs,
                               tStatus );

           // test if B-Spline coefficients exist
           if ( this->get_coefficients().length() > 0 )
           {
               // save coefficients
               save_matrix_to_hdf5_file(
                        tFileID,
                        "BSplineCoefficients",
                        this->get_coefficients(),
                        tStatus );

               // get mesh
               BSpline_Mesh_Base * tBMesh = mLagrangeMesh->get_bspline_mesh( mBSplineOrder );

               uint tNumberOfBasis = tBMesh->get_number_of_active_basis_on_proc();

               MORIS_ERROR( tNumberOfBasis == this->get_coefficients().length(),
                       "Something went wrong while saving B-Spline coeffs into HDF5 file" );

               // save IDs
               Matrix< IdMat > tBSplines( tNumberOfBasis, 1 );

               for( uint k=0; k<tNumberOfBasis; ++k )
               {
                   tBSplines( k ) = tBMesh->get_active_basis( k )->get_id();
               }

               // save coefficients
               save_matrix_to_hdf5_file(
                       tFileID,
                       "BSplineIDs",
                       tBSplines,
                       tStatus );
           }

           // close file
           tStatus = H5Fclose( tFileID );
        }

//------------------------------------------------------------------------------

        void
        Field::load_field_from_hdf5( const std::string & aFilePath )
        {
            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            // opens an existing file with read and write access
            hid_t tFileID = H5Fopen(
                    tFilePath.c_str(),
                    H5F_ACC_RDWR,
                    H5P_DEFAULT);

            // error handler
            hid_t tStatus;

            // par rank to read
            moris_id tRank;
            load_scalar_from_hdf5_file(
                    tFileID,
                    "ParRank",
                    tRank,
                    tStatus );

            // make sure that par rank is correct
            MORIS_ERROR(
                    tRank == par_rank(),
                    "Tried to read from HDF5 file with wrong par rank" );

            // test if label is not set
            if( this->get_label().size() == 0 )
            {
                std::string tLabel;
                load_string_from_hdf5_file(
                        tFileID,
                        "Label",
                        tLabel,
                        tStatus);

                // write label of this field
                this->set_label( tLabel );
            }

            // Matrix with IDs
            Matrix< IdMat > tNodeIDs;
            load_matrix_from_hdf5_file(
                    tFileID,
                    "NodeIDs",
                    tNodeIDs,
                    tStatus );

            // Orde of B-Splines
            load_scalar_from_hdf5_file(
                    tFileID,
                    "BSplineOrder",
                    mBSplineOrder,
                    tStatus );

            /// fixme: why is this uncommented?
            // test if B-Spline coefficients exist
            /*if( H5Lexists( tFileID, "BSplineCoefficients", H5P_DEFAULT ) )
            {

                // container for coefficients
                Matrix< DDRMat > & tCoeffs = this->get_coefficients();

                // load B-Spline coefficients
                load_matrix_from_hdf5_file(
                        tFileID,
                        "BSplineCoefficients",
                        tCoeffs,
                        tStatus );

                // load IDs
                Matrix< IdMat > tBSplineIDs;
                load_matrix_from_hdf5_file(
                        tFileID,
                        "BSplineCoefficients",
                        tBSplineIDs,
                        tStatus );

                // create map
                map< moris_id, real > tMap;

                // get number of basis
                uint tNumberOfBasis = tBSplineIDs.length();

                for( uint k=0; k<tNumberOfBasis; ++k )
                {
                    tMap[ tBSplineIDs( k ) ] = tCoeffs( k );
                }

                // get mesh
                BSpline_Mesh_Base * tBMesh = mLagrangeMesh->get_bspline_mesh();

                // copy data into coeffs
                for( uint k=0; k<tNumberOfBasis; ++k )
                {
                    tCoeffs( k ) = tMap.find( tBMesh->get_active_basis( k )->get_id() );
                }
            } */


            // get pointer to field data
            Matrix< DDRMat > & tNodeValues = this->get_node_values();

            load_matrix_from_hdf5_file(
                    tFileID,
                    "NodeValues",
                    tNodeValues,
                    tStatus );

            // close file
            tStatus = H5Fclose( tFileID );

            // create map
            map< moris_id, real > tMap;

            // get number of nodes
            uint tNumberOfNodes = tNodeIDs.length();

            // fill map with data
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                tMap[ tNodeIDs( k ) ] = tNodeValues( k );
            }

            // copy data into field
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                tNodeValues( k ) = tMap.find(
                        mMesh->get_glb_entity_id_from_entity_loc_index(
                                k, EntityRank::NODE ) );
            }
        }

//------------------------------------------------------------------------------


        void
        Field::save_node_values_to_binary( const std::string & aFilePath )
        {
            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            save_matrix_to_binary_file( this->get_node_values(), tFilePath );
        }

//------------------------------------------------------------------------------

        void
        Field::save_bspline_coeffs_to_binary( const std::string & aFilePath )
        {
            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            save_matrix_to_binary_file( this->get_coefficients(), tFilePath );
        }

//------------------------------------------------------------------------------

        void
        Field::project_coefficients()
        {

        }

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
