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

#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"
#include "cl_Library_IO.hpp"

namespace moris::ge
{
    class Voxel_Input : public Field_Analytic
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
        Voxel_Input(
                Vector_Type&              aADVs,
                Matrix< DDUMat >          aGeometryVariableIndices,
                Matrix< DDUMat >          aADVIndices,
                Matrix< DDRMat >          aConstants,
                std::string               aVoxelFieldName,
                Matrix< DDRMat >          aDomainDimensions,
                Matrix< DDRMat >          aDomainOffset,
                Matrix< DDRMat >          aGrainIdToValueMap,
                Level_Set_Parameters aParameters = Level_Set_Parameters() )
                : Field_Analytic( aADVs, aGeometryVariableIndices, aADVIndices, aConstants )
                , mDomainDimensions( aDomainDimensions )
                , mDomainOffset( aDomainOffset )
                , mGrainIdToValueMap( aGrainIdToValueMap )
        {
            this->read_voxel_data( aVoxelFieldName );
        }

        /**
         * Constructor with only constant parameters
         *
         * @param aConstants The constant field variables not filled by ADVs
         * @param aParameters Additional parameters
         */
        Voxel_Input(
                Matrix< DDRMat >          aConstants,
                std::string               aVoxelFieldName,
                Matrix< DDRMat >          aDomainDimensions,
                Matrix< DDRMat >          aDomainOffset,
                Matrix< DDRMat >          aGrainIdToValueMap,
                Level_Set_Parameters aParameters = Level_Set_Parameters() );

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

        /**
         * Return number of voxel of different color.
         */
        moris::uint
        get_num_voxel_Ids()
        {
            return mNumGrainInd;
        };

      private:
        void read_voxel_data( std::string aVoxelFieldName );

        real get_field_value_2d( const Matrix< DDRMat >& aCoordinates );
        real get_field_value_3d( const Matrix< DDRMat >& aCoordinates );
    };
}
