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

namespace moris::gen
{
    class Voxel_Input
    {

      private:
        std::string             mVoxelFileName;
        Node_Manager&           mNodeManager;

        Vector< real > mDomainDimensions;
        Vector< real > mDomainOffset;
        Matrix< DDUMat > mVoxelField;

        uint             mVoxelsInX;
        uint             mVoxelsInY;
        uint             mVoxelsInZ;

        uint mNumGrainInd;

      public:
        /**
         * Constructor
         */
        Voxel_Input(
                const std::string&    aVoxelFileName,
                const Vector< real >& aDomainDimensions,
                const Vector< real >& aDomainOffset,
                Node_Manager&         aNodeManager );

        /**
         * Gets the file name of the voxel field.
         *
         * @return File name
         */
        std::string get_file_name();

        /**
         * Gets the node manager being stored, so that each voxel geometry can use it.
         *
         * @return Node manager
         */
        const Node_Manager& get_node_manager();

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
        uint
        get_num_voxel_IDs()
        {
            return mNumGrainInd;
        };

      private:
        void read_voxel_data( const std::string& aVoxelFieldName );

        uint get_voxel_ID_2d( const Matrix< DDRMat >& aCoordinates );
        uint get_voxel_ID_3d( const Matrix< DDRMat >& aCoordinates );
    };
}
