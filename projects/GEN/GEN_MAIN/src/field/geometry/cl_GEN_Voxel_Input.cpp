#include "cl_GEN_Voxel_Input.hpp"

#include "cl_Ascii.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Voxel_Input::Voxel_Input(
                Matrix<DDRMat>            aConstants,
                std::string               aVoxelFieldName,
                Matrix<DDRMat>            aDomainDimensions,
                Matrix<DDRMat>            aDomainOffset,
                Geometry_Field_Parameters aParameters)
                : Field(aConstants, aParameters)
                , Geometry(aParameters)
                , mDomainDimensions( aDomainDimensions )
                , mDomainOffset( aDomainOffset )
        {
            this->read_voxel_data( aVoxelFieldName );
        }

        //--------------------------------------------------------------------------------------------------------------

        real Voxel_Input::get_field_value( const Matrix<DDRMat>& aCoordinates )
        {
            // set pointer to correct field value function
            if ( mDomainDimensions.numel() == 3 )
            {
                return get_field_value_3d(aCoordinates );
            }
            else
            {
                return get_field_value_2d(aCoordinates );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        real Voxel_Input::get_field_value_3d( const Matrix<DDRMat>& aCoordinates )
        {
            moris::real tVoxelSizeX = mDomainDimensions( 0 ) / mVoxelsInX;
            moris::real tVoxelSizeY = mDomainDimensions( 1 ) / mVoxelsInY;
            moris::real tVoxelSizeZ = mDomainDimensions( 2 ) / mVoxelsInZ;

            MORIS_ASSERT(
                    aCoordinates( 0 ) - mDomainOffset( 0 ) >= 0.0 &&
                    aCoordinates( 1 ) - mDomainOffset( 1 ) >= 0.0 &&
                    aCoordinates( 1 ) - mDomainOffset( 2 ) >= 0.0,
                    "Voxel_Input::get_field_value_2d() - invalid domain dimensions; check offset.\n");

            moris::uint tI = std::floor( ( aCoordinates( 0 ) - mDomainOffset( 0 ) ) / tVoxelSizeX );
            moris::uint tJ = std::floor( ( aCoordinates( 1 ) - mDomainOffset( 1 ) ) / tVoxelSizeY );
            moris::uint tK = std::floor( ( aCoordinates( 2 ) - mDomainOffset( 2 ) ) / tVoxelSizeZ );

            moris::real tEpsilon = 1E-12;

            if( aCoordinates( 0 ) >=  mDomainDimensions( 0 ) + mDomainOffset( 0  ) -tEpsilon )
            {
                tI = mVoxelsInX-1;
            }
            if( aCoordinates( 1 ) >=  mDomainDimensions( 1 ) + mDomainOffset( 1 ) -tEpsilon )
            {
                tJ = mVoxelsInY-1;
            }
            if( aCoordinates( 2 ) >=  mDomainDimensions( 2 ) + mDomainOffset( 2 ) -tEpsilon )
            {
                tK = mVoxelsInZ-1;
            }

            // mVoxelField columns are ordered - VoxelIndex - GainsId - I - J - K
            moris::uint tRow =
                    tI * mVoxelsInY * mVoxelsInZ +
                    tJ * mVoxelsInZ +
                    tK;

            return ( moris::real )mVoxelField( tRow, 1 );
        }
        //--------------------------------------------------------------------------------------------------------------

        real Voxel_Input::get_field_value_2d( const Matrix<DDRMat>& aCoordinates )
        {
            moris::real tVoxelSizeX = mDomainDimensions( 0 ) / mVoxelsInX;
            moris::real tVoxelSizeY = mDomainDimensions( 1 ) / mVoxelsInY;

            MORIS_ASSERT(
                    aCoordinates( 0 ) - mDomainOffset( 0 ) >= 0.0 &&
                    aCoordinates( 1 ) - mDomainOffset( 1 ) >= 0.0,
                    "Voxel_Input::get_field_value_2d() - invalid domain dimensions; check offset.\n");

            moris::uint tI = std::floor( ( aCoordinates( 0 ) - mDomainOffset( 0 ) ) / tVoxelSizeX );
            moris::uint tJ = std::floor( ( aCoordinates( 1 ) - mDomainOffset( 1 ) ) / tVoxelSizeY );

            moris::real tEpsilon = 1E-12;

            if( aCoordinates( 0 ) >=  mDomainDimensions( 0 ) + mDomainOffset( 0  ) -tEpsilon )
            {
                tI = mVoxelsInX-1;
            }
            if( aCoordinates( 1 ) >=  mDomainDimensions( 1 ) + mDomainOffset( 1 ) -tEpsilon )
            {
                tJ = mVoxelsInY-1;
            }

            // mVoxelField columns are ordered - VoxelIndex - GainsId - I - J
            moris::uint tRow = tI * mVoxelsInY + tJ;

            return ( moris::real )mVoxelField( tRow, 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Voxel_Input::get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates)
        {
            MORIS_ERROR( false,
                    "Voxel_Input::get_dfield_dadvs(), Sensitivities cannot be calculated for Voxel field.");
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Voxel_Input::get_dfield_dcoordinates(
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities)
        {
            MORIS_ERROR(false, "get_dfield_dcoordinates not implemented for voxel input geometry.");
        }

        //--------------------------------------------------------------------------------------------------------------

        void Voxel_Input::read_voxel_data( std::string aVoxelFieldName)
        {
            // build Ascii reader
            moris::Ascii tAsciiReader( aVoxelFieldName, moris::FileMode::OPEN_RDONLY );

            // get number of lines in asci file
            moris::uint tNumLines = tAsciiReader.length();

            // get number of spatial dimensions
            uint tNumSpaceDim = mDomainDimensions.numel();

            // set matrix for voxel field
            mVoxelField.set_size( tNumLines, tNumSpaceDim+2, MORIS_UINT_MAX );

            // initialize matrix for temporarily storing data of single voxel
            Matrix<DDUMat> tMatrix;

            for( uint Ik = 0; Ik < tNumLines; Ik++ )
            {
                // get line from ascii file
                std::string tFileLine = tAsciiReader.line( Ik );

                // convert line into numerical values
                delimited_string_to_mat(tFileLine," ",tMatrix);

                MORIS_ASSERT( tMatrix.numel() == mVoxelField.n_cols(),
                        "Voxel_Input::read_voxel_data - Incorrect number of columns in voxel file.\n");

                // store data
                mVoxelField.get_row( Ik ) = trans( tMatrix );
            }

            // check that all voxels have been initialized
            MORIS_ERROR( mVoxelField.max() != MORIS_UINT_MAX,
                    "Voxel_Input::Voxel_Input() - Matrix not correctly initialized.\n");

            // check that smallest grain id equals 1
            MORIS_ERROR( mVoxelField.get_column( 1 ).min() == 1,
                    "Voxel_Input::Voxel_Input() - Voxel index needs to be 1.\n");

            // get number of grains assuming that grains are numbered consecutively
            mNumGrainInd = mVoxelField.get_column( 1 ).max();

            // Voxel indices start with 1
            mVoxelsInX   = mVoxelField.get_column( 2 ).max() + 1;
            mVoxelsInY   = mVoxelField.get_column( 3 ).max() + 1;

            if (tNumSpaceDim == 3)
            {
                mVoxelsInZ   = mVoxelField.get_column( 4 ).max() + 1;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
