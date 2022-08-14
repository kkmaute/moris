#ifndef MORIS_CL_GEN_IMAGE_SDF_GEOMETRY_HPP
#define MORIS_CL_GEN_IMAGE_SDF_GEOMETRY_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"
#include "cl_Library_IO.hpp"

namespace moris
{
    namespace ge
    {
        class Image_SDF_Geometry : public Geometry
                , public Field_Analytic
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
             * @tparam Vector_Type Type of vector where ADVs are stored
             * @param aADVs ADV vector
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aParameters Additional parameters
             */
            template< typename Vector_Type >
            Image_SDF_Geometry(
                    Vector_Type&              aADVs,
                    Matrix< DDUMat >          aGeometryVariableIndices,
                    Matrix< DDUMat >          aADVIndices,
                    Matrix< DDRMat >          aConstants,
                    std::string               aImageFileName,
                    Matrix< DDRMat >          aDomainDimensions,
                    Matrix< DDRMat >          aDomainOffset,
                    real                      aSdfScaling,
                    real                      aSdfShift,
                    real                      aSdfDefault,
                    bool                      aSdfInterpolate,
                    Geometry_Field_Parameters aParameters = {} )
                    : Field( aADVs, aGeometryVariableIndices, aADVIndices, aConstants, aParameters )
                    , Geometry( aParameters )
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
             * @return Distance to this geometry
             */
            real get_field_value( const Matrix< DDRMat >& aCoordinates );

            /**
             * Given a node coordinate, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables.
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
    }    // namespace ge
}    // namespace moris

#endif    // MORIS_CL_GEN_IMAGE_SDF_GEOMETRY_HPP
