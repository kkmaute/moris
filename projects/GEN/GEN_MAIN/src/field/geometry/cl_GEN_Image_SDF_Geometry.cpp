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

            // compute relative position in image
            const real tIpos = ( aCoordinates( 0 ) - mDomainOffset( 0 ) ) / mVoxelSizeX;
            const real tJpos = ( aCoordinates( 1 ) - mDomainOffset( 1 ) ) / mVoxelSizeY;
            const real tKpos = ( aCoordinates( 2 ) - mDomainOffset( 2 ) ) / mVoxelSizeZ;

            real tInterpolatedValue;

            // linear interpolation
            if ( mDoInterpolate )
            {
                // compute i,j,k position of lower left voxel
                sint tI = std::floor( tIpos );
                sint tJ = std::floor( tJpos );
                sint tK = std::floor( tKpos );

                // clip i,j,k position (enforce i,j,k to be within domain)
                tI = std::max( 0, std::min( tI, (sint)mVoxelsInX - 2 ) );
                tJ = std::max( 0, std::min( tJ, (sint)mVoxelsInY - 2 ) );
                tK = std::max( 0, std::min( tK, (sint)mVoxelsInZ - 2 ) );

                const uint tRow = tK * mVoxelsInX * mVoxelsInY + tJ * mVoxelsInX + tI;

                // compute relative position in between voxels
                real tIloc = tIpos - (real)tI;
                real tJloc = tJpos - (real)tJ;
                real tKloc = tKpos - (real)tK;

                // clip relative position in between voxels and shift-scale relative positions
                tIloc = -1.0 + 2.0 * std::max( 0.0, std::min( tIloc, 1.0 ) );
                tJloc = -1.0 + 2.0 * std::max( 0.0, std::min( tJloc, 1.0 ) );
                tKloc = -1.0 + 2.0 * std::max( 0.0, std::min( tKloc, 1.0 ) );

                tInterpolatedValue                                                                                                                      //
                        = 0.125 * ( 1.0 - tIloc ) * ( 1.0 - tJloc ) * ( 1.0 - tKloc ) * mSdfField( tRow )                                               //
                        + 0.125 * ( 1.0 + tIloc ) * ( 1.0 - tJloc ) * ( 1.0 - tKloc ) * mSdfField( tRow + 1 )                                           //
                        + 0.125 * ( 1.0 + tIloc ) * ( 1.0 + tJloc ) * ( 1.0 - tKloc ) * mSdfField( tRow + mVoxelsInX + 1 )                              //
                        + 0.125 * ( 1.0 - tIloc ) * ( 1.0 + tJloc ) * ( 1.0 - tKloc ) * mSdfField( tRow + mVoxelsInX )                                  //
                        + 0.125 * ( 1.0 - tIloc ) * ( 1.0 - tJloc ) * ( 1.0 + tKloc ) * mSdfField( tRow + mVoxelsInX * mVoxelsInY )                     //
                        + 0.125 * ( 1.0 + tIloc ) * ( 1.0 - tJloc ) * ( 1.0 + tKloc ) * mSdfField( tRow + mVoxelsInX * mVoxelsInY + 1 )                 //
                        + 0.125 * ( 1.0 + tIloc ) * ( 1.0 + tJloc ) * ( 1.0 + tKloc ) * mSdfField( tRow + mVoxelsInX * mVoxelsInY + mVoxelsInX + 1 )    //
                        + 0.125 * ( 1.0 - tIloc ) * ( 1.0 + tJloc ) * ( 1.0 + tKloc ) * mSdfField( tRow + mVoxelsInX * mVoxelsInY + mVoxelsInX );
            }
            else
            {
                // compute i,j,k position of closest voxel
                sint tI = std::round( tIpos );
                sint tJ = std::round( tJpos );
                sint tK = std::round( tKpos );

                // clip i,j,k position (allow i,j,k  on boundaries)
                tI = std::max( 0, std::min( tI, (sint)mVoxelsInX - 1 ) );
                tJ = std::max( 0, std::min( tJ, (sint)mVoxelsInY - 1 ) );
                tK = std::max( 0, std::min( tK, (sint)mVoxelsInZ - 1 ) );

                const uint tRow = tK * mVoxelsInX * mVoxelsInY + tJ * mVoxelsInX + tI;

                tInterpolatedValue = mSdfField( tRow );
            }

            return ( tInterpolatedValue + mSdfShift );
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

            // compute relative position in image
            const real tIpos = ( aCoordinates( 0 ) - mDomainOffset( 0 ) ) / mVoxelSizeX;
            const real tJpos = ( aCoordinates( 1 ) - mDomainOffset( 1 ) ) / mVoxelSizeY;

            real tInterpolatedValue;

            if ( mDoInterpolate )
            {
                // compute i,j position of lower left voxel
                sint tI = std::floor( tIpos );
                sint tJ = std::floor( tJpos );

                // clip i,j position (enforce i,j to be within domain)
                tI = std::max( 0, std::min( tI, (sint)mVoxelsInX - 2 ) );
                tJ = std::max( 0, std::min( tJ, (sint)mVoxelsInY - 2 ) );

                // compute relative position in between voxels
                real tIloc = tIpos - (real)tI;
                real tJloc = tJpos - (real)tJ;

                // clip relative position in between voxels and shift-scale relative positions
                tIloc = -1.0 + 2.0 * std::max( 0.0, std::min( tIloc, 1.0 ) );
                tJloc = -1.0 + 2.0 * std::max( 0.0, std::min( tJloc, 1.0 ) );

                // linear interpolation
                const uint tRow = tJ * mVoxelsInX + tI;

                tInterpolatedValue                                                                         //
                        = 0.25 * ( 1.0 - tIloc ) * ( 1.0 - tJloc ) * mSdfField( tRow )                     //
                        + 0.25 * ( 1.0 + tIloc ) * ( 1.0 - tJloc ) * mSdfField( tRow + 1 )                 //
                        + 0.25 * ( 1.0 + tIloc ) * ( 1.0 + tJloc ) * mSdfField( tRow + mVoxelsInX + 1 )    //
                        + 0.25 * ( 1.0 - tIloc ) * ( 1.0 + tJloc ) * mSdfField( tRow + mVoxelsInX );
            }
            else
            {
                // compute i,j,k position of closest voxel
                sint tI = std::round( tIpos );
                sint tJ = std::round( tJpos );

                // clip i,j position (allow i,j,k  on boundaries)
                tI = std::max( 0, std::min( tI, (sint)mVoxelsInX - 1 ) );
                tJ = std::max( 0, std::min( tJ, (sint)mVoxelsInY - 1 ) );

                // linear interpolation
                const uint tRow = tJ * mVoxelsInX + tI;

                tInterpolatedValue = mSdfField( tRow );
            }

            return tInterpolatedValue + mSdfShift;
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

            mVoxelSizeX = mDomainDimensions( 0 ) / ( mVoxelsInX - 1 );
            mVoxelSizeY = mDomainDimensions( 1 ) / ( mVoxelsInY - 1 );

            if ( tDimensions.numel() > 2 )
            {
                mVoxelsInZ = tDimensions( 2 );

                mVoxelSizeZ = mDomainDimensions( 2 ) / ( mVoxelsInZ - 1 );
            }

            // check for proper size of sdf file (stored in single vector)
            MORIS_ERROR( mSdfField.numel() == mVoxelsInX * mVoxelsInY * mVoxelsInZ,
                    "Image_SDF_Geometry::read_image_sdf_data - mismatch of data size in SDF file." );

            // scale sdf file

            // check whether automatic scaling factor needs to be computed
            if ( std::abs( mSdfScaling ) < MORIS_REAL_EPS )
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

                mSdfScaling = tPhysicalDim / tSdfDifference;
            }

            mSdfField = mSdfScaling * mSdfField;

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
