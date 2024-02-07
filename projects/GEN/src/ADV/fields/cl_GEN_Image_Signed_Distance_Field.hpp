/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Image_Signed_Distance_Field.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Analytic.hpp"
#include "cl_Library_IO.hpp"

namespace moris::gen
{
    class Image_Signed_Distance_Field : public Field_Analytic< 0 >
    {
        /*   Image-based signed distance field
         *
         *   This field can be generated, for example via Matlab's bwdist function
         *
         *   pixel/voxel data stored in hdf5 file in single column format:
         *
         *   vector component of i,j pixel:    j*mVoxelInX+i
         *   vector component of i,j,k voxel:  k*(mVoxelInX+mVoxelInY)+j*mVoxelInX+i
         *
         *   note: i,j,k is zero based
         */

      private:
        Matrix< DDRMat > mDomainDimensions;    /// physical dimension of image
        Matrix< DDRMat > mDomainOffset;        /// offset of image

        Matrix< DDRMat > mSdfField;    /// SDF field (single column)

        uint mVoxelsInX = 1;    /// number of pixels/voxels in x-direction
        uint mVoxelsInY = 1;    /// number of pixels/voxels in y-direction
        uint mVoxelsInZ = 1;    /// number of pixels/voxels in z-direction

        real mVoxelSizeX = MORIS_REAL_MAX;    /// size of voxel in x-direction
        real mVoxelSizeY = MORIS_REAL_MAX;    /// size of voxel in y-direction
        real mVoxelSizeZ = MORIS_REAL_MAX;    /// size of voxel in z-direction

        real mSdfScaling = -1.0;              /// scaling of sdf field (0: automatic)
        real mSdfShift   = 0.0;               /// shift of sdf field
        real mSdfDefault = MORIS_REAL_MAX;    /// sdf value at points outside image

        bool mDoInterpolate = false;    /// flag to perform interpolation

      public:
        /**
         * Constructor, sets the pointers to advs and constant parameters for evaluations.
         *
         * @param aImageFileName Name of the image file to read
         * @param aDomainDimensions Dimensions of the domain
         * @param aDomainDimensions Amount to offset the domain by
         * @param aSdfScaling Amount to scale the SDF field by
         * @param aSdfShift Amount to shift the SDF field by
         * @param aSdfDefault Default SDF field outside of the domain
         */
        Image_Signed_Distance_Field(
                std::string               aImageFileName,
                Matrix< DDRMat >          aDomainDimensions,
                Matrix< DDRMat >          aDomainOffset,
                real                      aSdfScaling,
                real                      aSdfShift,
                real                      aSdfDefault,
                bool                      aSdfInterpolate )
                : Field_Analytic< 0 >( {} )
                , mDomainDimensions( aDomainDimensions )
                , mDomainOffset( aDomainOffset )
                , mSdfScaling( aSdfScaling )
                , mSdfShift( aSdfShift )
                , mSdfDefault( aSdfDefault )
                , mDoInterpolate( aSdfInterpolate )
        {
            this->read_image_sdf_data( aImageFileName );
        }

        /**
         * Given a node coordinate, returns the field value.
         *
         * @param aCoordinates Coordinate values
         * @return Distance to the interface
         */
        real get_field_value( const Matrix< DDRMat >& aCoordinates );

        /**
         * Given a node coordinate, evaluates the sensitivity of the field with respect to all of the
         * field variables.
         *
         * @param aCoordinates Coordinate values
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >& get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates );

        /**
         * Given nodal coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities );

      private:
        void read_image_sdf_data( std::string aImageFileName );

        real get_field_value_2d( const Matrix< DDRMat >& aCoordinates );
        real get_field_value_3d( const Matrix< DDRMat >& aCoordinates );
    };
}
