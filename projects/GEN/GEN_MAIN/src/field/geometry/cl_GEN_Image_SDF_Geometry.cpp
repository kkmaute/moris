#include "cl_GEN_Image_SDF_Geometry.hpp"

#include "HDF5_Tools.hpp"
#include "fn_Parsing_Tools.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        real
        Image_SDF_Geometry::get_field_value( const Matrix< DDRMat >& aCoordinates )
        {
            // set pointer to correct field value function
            if ( mDomainDimensions.numel() == 3 )
            {
                return get_field_value_3d( aCoordinates );
            }
            else
            {
                return get_field_value_2d( aCoordinates );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Image_SDF_Geometry::get_field_value_3d( const Matrix< DDRMat >& aCoordinates )
        {
            real tEpsilon = 1E-12;

            // check whether coordinates are outside image domain
            if ( aCoordinates( 0 ) - mDomainOffset( 0 ) + tEpsilon < 0.0 ) return mSdfDefault;
            if ( aCoordinates( 1 ) - mDomainOffset( 1 ) + tEpsilon < 0.0 ) return mSdfDefault;
            if ( aCoordinates( 2 ) - mDomainOffset( 2 ) + tEpsilon < 0.0 ) return mSdfDefault;

            if ( mDomainOffset( 0 ) + mDomainDimensions( 0 ) - aCoordinates( 0 ) + tEpsilon < 0 ) return mSdfDefault;
            if ( mDomainOffset( 1 ) + mDomainDimensions( 1 ) - aCoordinates( 1 ) + tEpsilon < 0 ) return mSdfDefault;
            if ( mDomainOffset( 2 ) + mDomainDimensions( 2 ) - aCoordinates( 2 ) + tEpsilon < 0 ) return mSdfDefault;

            // compute i,j,k position
            sint tI = std::floor( ( aCoordinates( 0 ) - mDomainOffset( 0 ) ) / mVoxelSizeX );
            sint tJ = std::floor( ( aCoordinates( 1 ) - mDomainOffset( 1 ) ) / mVoxelSizeY );
            sint tK = std::floor( ( aCoordinates( 2 ) - mDomainOffset( 2 ) ) / mVoxelSizeZ );

            // clip i, j, k
            tI = std::max( 0, std::min( tI, (sint)mVoxelsInX - 1 ) );
            tJ = std::max( 0, std::min( tJ, (sint)mVoxelsInY - 1 ) );
            tK = std::max( 0, std::min( tK, (sint)mVoxelsInZ - 1 ) );

            uint tRow = tK * mVoxelsInX * mVoxelsInY + tJ * mVoxelsInX + tI;

            return mSdfField( tRow ) + mSdfShift;
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Image_SDF_Geometry::get_field_value_2d( const Matrix< DDRMat >& aCoordinates )
        {
            real tEpsilon = 1E-12;

            // check whether coordinates are outside image domain
            if ( aCoordinates( 0 ) - mDomainOffset( 0 ) + tEpsilon < 0.0 ) return mSdfDefault;
            if ( aCoordinates( 1 ) - mDomainOffset( 1 ) + tEpsilon < 0.0 ) return mSdfDefault;

            if ( mDomainOffset( 0 ) + mDomainDimensions( 0 ) - aCoordinates( 0 ) + tEpsilon < 0 ) return mSdfDefault;
            if ( mDomainOffset( 1 ) + mDomainDimensions( 1 ) - aCoordinates( 1 ) + tEpsilon < 0 ) return mSdfDefault;

            // compute i,j position
            sint tI = std::floor( ( aCoordinates( 0 ) - mDomainOffset( 0 ) ) / mVoxelSizeX );
            sint tJ = std::floor( ( aCoordinates( 1 ) - mDomainOffset( 1 ) ) / mVoxelSizeY );

            // clip i, j
            tI = std::max( 0, std::min( tI, (sint)mVoxelsInX - 1 ) );
            tJ = std::max( 0, std::min( tJ, (sint)mVoxelsInY - 1 ) );

            uint tRow = tJ * mVoxelsInX + tI;

            return mSdfField( tRow ) + mSdfShift;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Image_SDF_Geometry::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
        {
            MORIS_ERROR( false,
                    "Image_SDF_Geometry::get_dfield_dadvs(), Sensitivities cannot be calculated for Voxel field." );
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Image_SDF_Geometry::get_dfield_dcoordinates(
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities )
        {
            MORIS_ERROR( false, "get_dfield_dcoordinates not implemented for voxel input geometry." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Image_SDF_Geometry::read_image_sdf_data( std::string aImageFiledName )
        {
            // open hdf5 file
            hid_t  tFileID = open_hdf5_file( aImageFiledName );
            herr_t tStatus = 0;

            // load image dimensions (number of pixels/voxels)
            Matrix< DDRMat > tDimensions;
            load_matrix_from_hdf5_file( tFileID, "Dimensions", tDimensions, tStatus );

            // load sdf from file
            load_matrix_from_hdf5_file( tFileID, "SDF", mSdfField, tStatus );

            // close file
            close_hdf5_file( tFileID );

            // check for correct dimensions
            MORIS_ERROR( mDomainDimensions.numel() == tDimensions.numel(),
                    "Image_SDF_Geometry::read_image_sdf_data - mismatch of dimensions in SDF file and domain dimensions." );

            // extract dimensions and compute size of voxel
            mVoxelsInX = tDimensions( 0 );
            mVoxelsInY = tDimensions( 1 );

            mVoxelSizeX = mDomainDimensions( 0 ) / mVoxelsInX;
            mVoxelSizeY = mDomainDimensions( 1 ) / mVoxelsInY;

            if ( tDimensions.numel() > 2 )
            {
                mVoxelsInZ = tDimensions( 2 );

                mVoxelSizeY = mDomainDimensions( 2 ) / mVoxelsInZ;
            }

            // check for proper size of sdf file (stored in single vector)
            MORIS_ERROR( mSdfField.numel() == mVoxelsInX * mVoxelsInY * mVoxelsInZ,
                    "Image_SDF_Geometry::read_image_sdf_data - mismatch of data size in SDF file." );

            // scale sdf file

            // check whether automatic scaling factor needs to be computed
            if ( std::abs( mSDFScaling ) < MORIS_REAL_EPS )
            {
                // compute range of sdf values
                real tSdfDifference = mSdfField.max() - mSdfField.min();

                // compute size measure of physical domain of image
                real tPhysicalDim;
                if ( tDimensions.numel() > 2 )
                {
                    tPhysicalDim = std::pow( mDomainDimensions( 0 ) * mDomainDimensions( 1 ) * mDomainDimensions( 2 ), 1.0 / 3.0 );
                }
                else
                {
                    tPhysicalDim = std::pow( mDomainDimensions( 0 ) * mDomainDimensions( 1 ), 1.0 / 2.0 );
                }

                mSDFScaling = tPhysicalDim / tSdfDifference;
            }

            mSdfField = mSDFScaling * mSdfField;

            // compute default value
            if ( mSdfDefault < 0 )
            {
                mSdfDefault = mSdfField.min();
            }
            else
            {
                mSdfDefault = mSdfField.max();
            }
        }

        //--------------------------------------------------------------------------------------------------------------
    }    // namespace ge
}    // namespace moris
