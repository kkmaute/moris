/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Plane.cpp
 *
 */

#include "cl_GEN_Plane.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Plane::Plane(
            const ADV&  aXCenter,
            const ADV&  aYCenter,
            const ADV&  aZCenter,
            const ADV&  aXNormal,
            const ADV&  aYNormal,
            const ADV&  aZNormal,
            std::string aName )
            : Field_Analytic< 3 >( { aXCenter, aYCenter, aZCenter, aXNormal, aYNormal, aZNormal }, std::move( aName ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Plane::get_field_value( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        real tXCenter = mADVHandler.get_variable( 0 );
        real tYCenter = mADVHandler.get_variable( 1 );
        real tZCenter = mADVHandler.get_variable( 2 );
        real tXNormal = mADVHandler.get_variable( 3 );
        real tYNormal = mADVHandler.get_variable( 4 );
        real tZNormal = mADVHandler.get_variable( 5 );

        // Evaluate field value
        return tXNormal * ( aCoordinates( 0 ) - tXCenter ) + tYNormal * ( aCoordinates( 1 ) - tYCenter ) + tZNormal * ( aCoordinates( 2 ) - tZCenter );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Plane::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        real tXCenter = mADVHandler.get_variable( 0 );
        real tYCenter = mADVHandler.get_variable( 1 );
        real tZCenter = mADVHandler.get_variable( 2 );
        real tXNormal = mADVHandler.get_variable( 3 );
        real tYNormal = mADVHandler.get_variable( 4 );
        real tZNormal = mADVHandler.get_variable( 5 );

        // Evaluate sensitivities
        mSensitivities( 0 ) = -tXNormal;
        mSensitivities( 1 ) = -tYNormal;
        mSensitivities( 2 ) = -tZNormal;
        mSensitivities( 3 ) = aCoordinates( 0 ) - tXCenter;
        mSensitivities( 4 ) = aCoordinates( 1 ) - tYCenter;
        mSensitivities( 5 ) = aCoordinates( 2 ) - tZCenter;

        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Plane::get_dfield_dcoordinates(
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        aSensitivities( 0 ) = mADVHandler.get_variable( 3 );
        aSensitivities( 1 ) = mADVHandler.get_variable( 4 );
        aSensitivities( 2 ) = mADVHandler.get_variable( 5 );
    }

    //--------------------------------------------------------------------------------------------------------------
}
