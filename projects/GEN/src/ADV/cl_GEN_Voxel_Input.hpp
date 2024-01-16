/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Voxel_Input.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Analytic.hpp"
#include "cl_Library_IO.hpp"

namespace moris::ge
{
    class Voxel_Input
    {

      private:
        moris::Matrix< DDRMat > mDomainDimensions;
        moris::Matrix< DDRMat > mDomainOffset;

        moris::Matrix< DDRMat > mGrainIdToValueMap;

        moris::Matrix< DDUMat > mVoxelField;
        moris::uint             mVoxelsInX;
        moris::uint             mVoxelsInY;
        moris::uint             mVoxelsInZ;

        moris::uint mNumGrainInd;

      public:
        /**
         * Constructor
         */
        Voxel_Input(
                std::string               aVoxelFieldName,
                Matrix< DDRMat >          aDomainDimensions,
                Matrix< DDRMat >          aDomainOffset,
                Matrix< DDRMat >          aGrainIdToValueMap );

        /**
         * Given a node coordinate, returns the field value.
         *
         * @param aCoordinates Coordinate values
         * @return Field value
         */
        uint get_voxel_ID( const Matrix< DDRMat >& aCoordinates );

        /**
         * Return number of voxel of different color.
         */
        moris::uint
        get_num_voxel_IDs()
        {
            return mNumGrainInd;
        };

      private:
        void read_voxel_data( std::string aVoxelFieldName );

        uint get_voxel_ID_2d( const Matrix< DDRMat >& aCoordinates );
        uint get_voxel_ID_3d( const Matrix< DDRMat >& aCoordinates );
    };
}
