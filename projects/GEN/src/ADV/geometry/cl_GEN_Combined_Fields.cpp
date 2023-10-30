/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Combined_Fields.cpp
 *
 */

#include "cl_GEN_Combined_Fields.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Combined_Fields::Combined_Fields(
            Cell< std::shared_ptr< Field > > aFields,
            bool                             aUseMinimum )
            : Field( {{}} )
            , mFields( aFields )
            , mScale( 2 * aUseMinimum - 1 )
    {
        MORIS_ERROR( mFields.size() > 0, "A GEN Combined_Fields must be created with at least one field." );
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Combined_Fields::get_field_value(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates)
    {
        real tResult = mScale * mFields( 0 )->get_field_value( aNodeIndex, aCoordinates );
        for (uint tFieldIndex = 1; tFieldIndex < mFields.size(); tFieldIndex++)
        {
            tResult = std::min( tResult, mScale * mFields( tFieldIndex )->get_field_value( aNodeIndex, aCoordinates ) );
        }
        return mScale * tResult;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>&
    Combined_Fields::get_dfield_dadvs(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates)
    {
        // Find which field is the minimum
        real tMin = mScale * mFields( 0 )->get_field_value( aNodeIndex, aCoordinates );
        uint tMinFieldIndex = 0;
        for ( uint tFieldIndex = 1; tFieldIndex < mFields.size(); tFieldIndex++ )
        {
            real tResult = mScale * mFields( tFieldIndex )->get_field_value( aNodeIndex, aCoordinates );
            if ( tResult < tMin )
            {
                tMin = tResult;
                tMinFieldIndex = tFieldIndex;
            }
        }

        // Return relevant sensitivity
        return mFields( tMinFieldIndex )->get_dfield_dadvs( aNodeIndex, aCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Combined_Fields::get_dfield_dcoordinates(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates,
            Matrix<DDRMat>&       aSensitivities)
    {
        // Find which field is the minimum
        real tMin = mScale * mFields( 0 )->get_field_value( aNodeIndex, aCoordinates );
        uint tMinFieldIndex = 0;
        for ( uint tFieldIndex = 1; tFieldIndex < mFields.size(); tFieldIndex++ )
        {
            real tResult = mScale * mFields( tFieldIndex )->get_field_value( aNodeIndex, aCoordinates );
            if ( tResult < tMin )
            {
                tMin = tResult;
                tMinFieldIndex = tFieldIndex;
            }
        }

        // Get relevant sensitivity
        mFields( tMinFieldIndex )->get_dfield_dcoordinates( aNodeIndex, aCoordinates, aSensitivities );
    }

    //--------------------------------------------------------------------------------------------------------------

}
