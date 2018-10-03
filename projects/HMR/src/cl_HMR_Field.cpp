
// HD5 c-interface
#include "hdf5.h"

#include "fn_save_matrix_to_binary_file.hpp"

#include "HMR_Tools.hpp"
#include "HMR_HDF5_Tools.hpp"

#include "cl_HMR_Field.hpp"
#include "cl_HMR_Block.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        Field::Field(
                const std::string             & aLabel,
                std::shared_ptr< mtk::Mesh >    aMesh,
                std::shared_ptr< Database >     aDatabase,
                Lagrange_Mesh_Base *            aLagrangeMesh ) :
                        mtk::Field( aLabel, aMesh ),
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

            //! error handler
            herr_t tStatus;

            // save label of this field
            save_string_to_hdf5_file(
                    tFileID,
                    "Label",
                    this->get_label(),
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
                    mLagrangeMesh->get_bspline_order(),
                    tStatus );

            // save node values
           save_matrix_to_hdf5_file(
                    tFileID,
                    "NodeValues",
                    this->get_node_values(),
                    tStatus );

            // save coefficients
            save_matrix_to_hdf5_file(
                    tFileID,
                    "BSplineCoefficients",
                    this->get_coefficients(),
                    tStatus );

            // close file
            tStatus = H5Fclose( tFileID );
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
    } /* namespace hmr */
} /* namespace moris */
