/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Single_Grain.cpp
 *
 */

#include "cl_GEN_Single_Grain.hpp"
#include "cl_GEN_Voxel_Input.hpp"

#include "cl_Ascii.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Single_Grain::Single_Grain(
            std::shared_ptr< Field > aVoxelGeometry,
            uint                      aIndex )
            : Field_Analytic( {} )
            , mVoxelGeometry(aVoxelGeometry)
            , mIndex(aIndex)
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    real Single_Grain::get_field_value( const Matrix<DDRMat>& aCoordinates )
    {
        moris::real tGeoValue = reinterpret_cast< Voxel_Input* >(mVoxelGeometry.get())->get_field_value(aCoordinates);

        moris::real tValue = -1.0;

        if( tGeoValue == ( real ) mIndex )
        {
            tValue =  1.0;
        }

        return tValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>& Single_Grain::get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates)
    {
        MORIS_ERROR( false, "Voxel_Input::get_dfield_dadvs(), Sensitivities cannot be calculated for Voxel field.");
        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Single_Grain::get_dfield_dcoordinates(
            const Matrix<DDRMat>& aCoordinates,
            Matrix<DDRMat>&       aSensitivities)
    {
        MORIS_ERROR(false, "get_dfield_dcoordinates not implemented for single grain field.");
    }

    //--------------------------------------------------------------------------------------------------------------

}
