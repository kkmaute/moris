/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Voxel_Input.cpp
 *
 */

#include "cl_GEN_Voxel_Input.hpp"

#include "cl_Ascii.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_trans.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Voxel_Input::Voxel_Input(
            std::string      aVoxelFileName,
            Matrix< DDRMat > aDomainDimensions,
            Matrix< DDRMat > aDomainOffset,
            Matrix< DDRMat > aGrainIdToValueMap,
            Node_Manager&    aNodeManager )
            : mVoxelFileName( aVoxelFileName )
            , mNodeManager( aNodeManager )
            , mDomainDimensions( aDomainDimensions )
            , mDomainOffset( aDomainOffset )
            , mGrainIdToValueMap( aGrainIdToValueMap )
    {
        this->read_voxel_data( aVoxelFileName );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Voxel_Input::get_file_name()
    {
        return mVoxelFileName;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Node_Manager& Voxel_Input::get_node_manager()
    {
        return mNodeManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Voxel_Input::get_voxel_ID( const Matrix< DDRMat >& aCoordinates )
    {
        // set pointer to correct field value function
        if ( mDomainDimensions.numel() == 3 )
        {
            return get_voxel_ID_3d( aCoordinates );
        }
        else
        {
            return get_voxel_ID_2d( aCoordinates );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Voxel_Input::get_voxel_ID_3d( const Matrix< DDRMat >& aCoordinates )
    {
        moris::real tVoxelSizeX = mDomainDimensions( 0 ) / mVoxelsInX;
        moris::real tVoxelSizeY = mDomainDimensions( 1 ) / mVoxelsInY;
        moris::real tVoxelSizeZ = mDomainDimensions( 2 ) / mVoxelsInZ;

        MORIS_ASSERT(
                aCoordinates( 0 ) - mDomainOffset( 0 ) >= 0.0               //
                        && aCoordinates( 1 ) - mDomainOffset( 1 ) >= 0.0    //
                        && aCoordinates( 1 ) - mDomainOffset( 2 ) >= 0.0,
                "Voxel_Input::get_voxel_ID_3d() - invalid domain dimensions; check offset.\n" );

        moris::uint tI = std::floor( ( aCoordinates( 0 ) - mDomainOffset( 0 ) ) / tVoxelSizeX );    // K
        moris::uint tJ = std::floor( ( aCoordinates( 1 ) - mDomainOffset( 1 ) ) / tVoxelSizeY );
        moris::uint tK = std::floor( ( aCoordinates( 2 ) - mDomainOffset( 2 ) ) / tVoxelSizeZ );    // I

        moris::real tEpsilon = 1E-12;

        if ( aCoordinates( 0 ) >= mDomainDimensions( 0 ) + mDomainOffset( 0 ) - tEpsilon )
        {
            tI = mVoxelsInX - 1;    // K
        }
        if ( aCoordinates( 1 ) >= mDomainDimensions( 1 ) + mDomainOffset( 1 ) - tEpsilon )
        {
            tJ = mVoxelsInY - 1;
        }
        if ( aCoordinates( 2 ) >= mDomainDimensions( 2 ) + mDomainOffset( 2 ) - tEpsilon )
        {
            tK = mVoxelsInZ - 1;    // I
        }

        // mVoxelInput columns are ordered - VoxelIndex - GainsId - I - J - K
        moris::uint tRow =
                tI * mVoxelsInY * mVoxelsInZ + tJ * mVoxelsInZ + tK;

        return mVoxelField( tRow, 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Voxel_Input::get_voxel_ID_2d( const Matrix< DDRMat >& aCoordinates )
    {
        int tGrainID;

        // check whether coordinates are within voxel domain
        if ( aCoordinates( 0 ) - mDomainOffset( 0 ) >= 0.0                                 //
                && mDomainOffset( 0 ) + mDomainDimensions( 0 ) - aCoordinates( 0 ) >= 0    //
                && aCoordinates( 1 ) - mDomainOffset( 1 ) >= 0.0                           //
                && mDomainOffset( 1 ) + mDomainDimensions( 1 ) - aCoordinates( 1 ) >= 0 )
        {
            moris::real tVoxelSizeX = mDomainDimensions( 0 ) / mVoxelsInX;
            moris::real tVoxelSizeY = mDomainDimensions( 1 ) / mVoxelsInY;

            moris::uint tI = std::floor( ( aCoordinates( 0 ) - mDomainOffset( 0 ) ) / tVoxelSizeX );
            moris::uint tJ = std::floor( ( aCoordinates( 1 ) - mDomainOffset( 1 ) ) / tVoxelSizeY );

            moris::real tEpsilon = 1E-12;

            if ( aCoordinates( 0 ) >= mDomainDimensions( 0 ) + mDomainOffset( 0 ) - tEpsilon )
            {
                tI = mVoxelsInX - 1;
            }
            if ( aCoordinates( 1 ) >= mDomainDimensions( 1 ) + mDomainOffset( 1 ) - tEpsilon )
            {
                tJ = mVoxelsInY - 1;
            }

            // mVoxelInput columns are ordered - VoxelIndex - GainsId - I - J
            moris::uint tRow = tI * mVoxelsInY + tJ;

            tGrainID = mVoxelField( tRow, 1 );
        }
        else
        {
            tGrainID = mNumGrainInd + 1;
        }

        // return mapped value if map exists
        if ( mGrainIdToValueMap.numel() > 0 )
        {
            return mGrainIdToValueMap( tGrainID - 1 );
        }

        return tGrainID;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Voxel_Input::read_voxel_data( std::string aVoxelFieldName )
    {
        // build Ascii reader
        moris::Ascii tAsciiReader( aVoxelFieldName, moris::FileMode::OPEN_RDONLY );

        // get number of lines in asci file
        moris::uint tNumLines = tAsciiReader.length();

        // get number of spatial dimensions
        uint tNumSpaceDim = mDomainDimensions.numel();

        // set matrix for voxel field
        mVoxelField.set_size( tNumLines, tNumSpaceDim + 2, MORIS_UINT_MAX );

        // initialize matrix for temporarily storing data of single voxel
        Matrix< DDUMat > tMatrix;

        for ( uint Ik = 0; Ik < tNumLines; Ik++ )
        {
            // get line from ascii file
            std::string tFileLine = tAsciiReader.line( Ik );

            // convert line into numerical values
            delimited_string_to_mat( tFileLine, " ", tMatrix );

            MORIS_ASSERT( tMatrix.numel() == mVoxelField.n_cols(),
                    "Voxel_Input::read_voxel_data - Incorrect number of columns in voxel file.\n" );

            // store data
            mVoxelField.get_row( Ik ) = trans( tMatrix );
        }

        // check that all voxels have been initialized
        MORIS_ERROR( mVoxelField.max() != MORIS_UINT_MAX,
                "Voxel_Input::Voxel_Input() - Matrix not correctly initialized.\n" );

        // check that smallest grain id equals 1
        MORIS_ERROR( mVoxelField.get_column( 1 ).min() == 1,
                "Voxel_Input::Voxel_Input() - Voxel index needs to be 1.\n" );

        // get number of grains assuming that grains are numbered consecutively
        mNumGrainInd = mVoxelField.get_column( 1 ).max();

        // Voxel indices start with 1

        if ( tNumSpaceDim == 3 )
        {
            mVoxelsInZ = mVoxelField.get_column( 2 ).max() + 1;
            mVoxelsInY = mVoxelField.get_column( 3 ).max() + 1;
            mVoxelsInX = mVoxelField.get_column( 4 ).max() + 1;
        }
        else
        {
            mVoxelsInY = mVoxelField.get_column( 3 ).max() + 1;
            mVoxelsInX = mVoxelField.get_column( 2 ).max() + 1;
        }
    }

    //--------------------------------------------------------------------------------------------------------------
}